/* Native ISA transcript evaluator. Sources: Intel SDM ch.8 and instruction
 * definitions; original implementation, MIT. No libm or long-double arithmetic.
 * Every temporary is ten bytes. A tbyte round trip preserves a register value;
 * only an explicitly selected qword store rounds to binary64. */
#include <stdint.h>
#include <string.h>

typedef struct { unsigned char b[10]; } Ext;
typedef struct { int32_t op, a, b, site, phase; uint64_t bits; } Node;
enum { INPUT, CONST, ADD, SUB, MUL, DIV, NEG, LOG2, LOG2P1, EXPM1_2,
       RNDINT, SCALE, LN2, L2E, SQRT };

static Ext load64(uint64_t bits) {
    Ext out;
    __asm__ volatile("fldl %1; fstpt %0" : "=m"(out) : "m"(bits) : "st");
    return out;
}
static uint64_t store64(Ext x) {
    uint64_t out;
    __asm__ volatile("fldt %1; fstpl %0" : "=m"(out) : "m"(x) : "st");
    return out;
}
#define BIN_CASE(code, instruction) case code: \
    __asm__ volatile("fldt %1; fldt %2; " instruction "; fstpt %0" \
        : "=m"(z) : "m"(a), "m"(b) : "st", "st(1)"); break
#define UN_CASE(code, instruction) case code: \
    __asm__ volatile("fldt %1; " instruction "; fstpt %0" \
        : "=m"(z) : "m"(a) : "st"); break
static Ext instruction(int op, Ext a, Ext b) {
    Ext z;
    switch(op) {
        BIN_CASE(ADD, "faddp");
        BIN_CASE(SUB, "fsubrp");
        BIN_CASE(MUL, "fmulp");
        BIN_CASE(DIV, "fdivrp");
        /* FYL2X/YL2XP1: st1=a (multiplier), st0=b (argument). */
        BIN_CASE(LOG2, "fyl2x");
        BIN_CASE(LOG2P1, "fyl2xp1");
        UN_CASE(NEG, "fchs");
        UN_CASE(EXPM1_2, "f2xm1");
        UN_CASE(RNDINT, "frndint");
        UN_CASE(SQRT, "fsqrt");
        case SCALE:
            __asm__ volatile("fldt %2; fldt %1; fscale; fstpt %0; fstp %%st(0)"
                : "=m"(z) : "m"(a), "m"(b) : "st", "st(1)");
            break;
        case LN2:
            __asm__ volatile("fldln2; fstpt %0" : "=m"(z) : : "st"); break;
        case L2E:
            __asm__ volatile("fldl2e; fstpt %0" : "=m"(z) : : "st"); break;
        default: memset(&z, 0, sizeof z); break;
    }
    return z;
}

/* raw/values may be NULL. Outputs always contain all ten bytes. Return -1 for
 * malformed graph, never a successful or silently omitted candidate. */
int x87_eval(const Node *nodes, int count, const uint64_t *inputs, int input_count,
             const Ext *extended, int extended_count, uint32_t words, uint64_t mask,
             int out_a, int out_b, Ext *outputs, Ext *raw, Ext *values) {
    if (count <= 0 || count > 128 || out_a < 0 || out_a >= count || out_b < 0 || out_b >= count)
        return -1;
    Ext v[128];
    uint16_t cw = (uint16_t)words, inner = (uint16_t)(words >> 16);
    if ((cw != 0x123f && cw != 0x133f) || (inner != 0x123f && inner != 0x133f)) return -1;
    uint16_t saved;
    __asm__ volatile("fnstcw %0; fldcw %1" : "=m"(saved) : "m"(cw));
    int rc = 0;
    for (int i=0; i<count; ++i) {
        Node n=nodes[i]; Ext z;
        if (n.phase < 0 || n.phase > 1) { rc=-1; break; }
        uint16_t next = n.phase ? inner : (uint16_t)words;
        if (cw != next) { cw=next; __asm__ volatile("fldcw %0" : : "m"(cw)); }
        if (n.op == INPUT) {
            if (n.a < 0 || n.a >= input_count+extended_count) { rc=-1; break; }
            z = n.a < input_count ? load64(inputs[n.a]) : extended[n.a-input_count];
        } else if (n.op == CONST) z=load64(n.bits);
        else {
            if (n.op < ADD || n.op > SQRT || n.a < 0 || n.a >= i || n.b < 0 || n.b >= i) { rc=-1; break; }
            z=instruction(n.op,v[n.a],v[n.b]);
        }
        if (raw) raw[i]=z;
        if (n.site >= 0) {
            if (n.site >= 63) { rc=-1; break; }
            if (mask & (UINT64_C(1)<<n.site)) z=load64(store64(z));
        }
        v[i]=z;
        if (values) values[i]=z;
    }
    if (!rc) { outputs[0]=v[out_a]; outputs[1]=v[out_b]; }
    __asm__ volatile("fldcw %0" : : "m"(saved));
    return rc;
}

/* Batch interfaces remove Python overhead without changing the graph. */
int x87_masks(const Node *nodes, int count, const uint64_t *inputs, int ni,
              const Ext *extended, int ne, uint32_t words, uint64_t first,
              uint64_t number, int oa, int ob, Ext *outputs) {
    for (uint64_t j=0; j<number; ++j)
        if (x87_eval(nodes,count,inputs,ni,extended,ne,words,first+j,oa,ob,outputs+2*j,0,0)) return -1;
    return 0;
}

/* Same interpreter and publication, over a row-major binary64 input grid. */
int x87_grid(const Node *nodes, int count, const uint64_t *inputs, int ni,
             uint64_t rows, uint32_t words, uint64_t mask, int oa, int ob,
             uint64_t *outputs) {
    if (ni < 1) return -1;
    for (uint64_t j=0; j<rows; ++j) {
        Ext pair[2];
        if (x87_eval(nodes,count,inputs+j*ni,ni,0,0,words,mask,oa,ob,pair,0,0)) return -1;
        uint16_t saved, cw=(uint16_t)words;
        __asm__ volatile("fnstcw %0; fldcw %1" : "=m"(saved) : "m"(cw));
        outputs[2*j]=store64(pair[0]); outputs[2*j+1]=store64(pair[1]);
        __asm__ volatile("fldcw %0" : : "m"(saved));
    }
    return 0;
}

/* Independent callers can replay a single native instruction, avoiding the
 * interpreter's indexing and store-mask logic. */
void x87_leaf(int op, const Ext *a, const Ext *b, uint16_t cw, Ext *out) {
    uint16_t saved;
    __asm__ volatile("fnstcw %0; fldcw %1" : "=m"(saved) : "m"(cw));
    *out=instruction(op,*a,*b);
    __asm__ volatile("fldcw %0" : : "m"(saved));
}
uint64_t x87_publish(const Ext *a, uint16_t cw) {
    uint16_t saved;
    __asm__ volatile("fnstcw %0; fldcw %1" : "=m"(saved) : "m"(cw));
    uint64_t result=store64(*a);
    __asm__ volatile("fldcw %0" : : "m"(saved));
    return result;
}
