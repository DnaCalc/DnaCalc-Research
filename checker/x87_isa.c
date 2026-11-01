/* Independent native microcode replay for the rational certificate checker.
 * Original implementation from public ISA instruction definitions, MIT.
 * No graph interpreter, candidate generator or math library is linked. */
#include <stdint.h>
typedef struct { unsigned char bytes[10]; } Register;

int replay(int opcode, const Register *a, const Register *b, uint16_t word, Register *result) {
    uint16_t original;
    __asm__ volatile("fnstcw %0; fldcw %1" : "=m"(original) : "m"(word));
    int supported=1;
    switch (opcode) {
    case 7:
        __asm__ volatile("fldt %1; fldt %2; fyl2x; fstpt %0"
            : "=m"(*result) : "m"(*a), "m"(*b) : "st", "st(1)"); break;
    case 8:
        __asm__ volatile("fldt %1; fldt %2; fyl2xp1; fstpt %0"
            : "=m"(*result) : "m"(*a), "m"(*b) : "st", "st(1)"); break;
    case 9:
        __asm__ volatile("fldt %1; f2xm1; fstpt %0"
            : "=m"(*result) : "m"(*a) : "st"); break;
    case 12:
        __asm__ volatile("fldln2; fstpt %0" : "=m"(*result) : : "st"); break;
    case 13:
        __asm__ volatile("fldl2e; fstpt %0" : "=m"(*result) : : "st"); break;
    default: supported=0; break;
    }
    __asm__ volatile("fldcw %0" : : "m"(original));
    return supported;
}
