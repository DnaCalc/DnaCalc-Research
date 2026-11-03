"""Independent exact checker for R13 x87 instruction/store box certificates.

Only integer/rational rounding code is shared with earlier certificate checkers.
No candidate generator, interpreter, numeric library or product module is imported.
Transcendental leaves are replayed as instructions on the recorded native CPU;
their values are never assumed to be correctly rounded mathematical functions.
"""
import argparse
import base64
from collections import Counter
import ctypes as C
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

from exact import (BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P,
                   binary64, exponent_of, two)

ROOT = Path(__file__).resolve().parent.parent
WORDS = (0x123f, 0x133f)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def records(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines()]


@lru_cache(maxsize=262144)
def decode(raw):
    require(len(raw) == 10, 'not an 80-bit value')
    sig, se = int.from_bytes(raw[:8], 'little'), int.from_bytes(raw[8:], 'little')
    exp = se & 0x7fff
    require(exp < 0x7fff, 'nonfinite value outside finite certificate scope')
    require(bool(sig >> 63) == bool(exp), 'noncanonical extended encoding')
    value = sig * two((exp or 1)-16383-63)
    return -value if se & 0x8000 else value


def negative(raw):
    return bool(raw[9] & 128)


@lru_cache(maxsize=262144)
def encode(value, negative_zero=False):
    sign = value < 0 or (value == 0 and negative_zero)
    value = abs(value)
    if value == 0:
        return bytes(8)+(0x8000 if sign else 0).to_bytes(2, 'little')
    e = max(-16382, exponent_of(value))
    m = value/two(e-63)
    require(m.denominator == 1 and m < 1 << 64, 'not an exact extended value')
    return int(m).to_bytes(8, 'little')+((e+16383 if m >= 1 << 63 else 0)
            | (0x8000 if sign else 0)).to_bytes(2, 'little')


def bits64(value, negative_zero=False):
    sign = value < 0 or (value == 0 and negative_zero)
    value = abs(value)
    if not value:
        return (1 << 63) if sign else 0
    e = max(-1022, exponent_of(value))
    m = value/two(e-52)
    require(m.denominator == 1 and m < 1 << 53, 'not a binary64 result')
    return ((1 << 63) if sign else 0) | ((e+1023) << 52 if m >= 1 << 52 else 0) | (int(m) & ((1 << 52)-1))


@lru_cache(maxsize=4096)
def load64(bits):
    return encode(binary64(bits),bool(bits>>63))


@lru_cache(maxsize=65536)
def loaded_arguments(inputs,extended):
    require(len(extended)%10==0,'incomplete extended argument')
    return tuple(load64(x)for x in inputs)+tuple(extended[i:i+10]for i in range(0,len(extended),10))


@lru_cache(maxsize=262144)
def published_bits(raw):
    value=B.round(decode(raw))
    require(value is not None,'publication overflow')
    return bits64(value,negative(raw))


class Register(C.Structure):
    _fields_ = [('bytes', C.c_ubyte*10)]


class Instructions:
    def __init__(self, library):
        self.library = Path(library)
        self.lib = C.CDLL(str(self.library.resolve()))
        ptr = C.POINTER(Register)
        self.lib.replay.argtypes = [C.c_int, ptr, ptr, C.c_uint16, ptr]
        self.counts = Counter()
        self.transcript = {}

    @lru_cache(maxsize=262144)
    def run(self, kind, a, b, cw):
        x, y, z = Register.from_buffer_copy(a), Register.from_buffer_copy(b), Register()
        opcode = {'fyl2x':7, 'fyl2xp1':8, 'f2xm1':9, 'fldln2':12, 'fldl2e':13}[kind]
        require(self.lib.replay(opcode, C.byref(x), C.byref(y), cw, C.byref(z)) == 1,
                'unsupported hardware leaf')
        self.counts[kind] += 1
        key, result = (kind, a, b, cw), bytes(z)
        prior = self.transcript.get(key)
        require(prior is None or prior == result, 'native leaves differ within the replay')
        self.transcript[key] = result
        return result


def canonical_templates():
    """Derive the declared ISA grammar independently of the search generator."""
    def context(n):
        return [dict(op='input', a=i) for i in range(n)]+[dict(op='const', bits='3ff0000000000000')]

    def append(ns, op, a=0, b=0, stored=True, phase=0):
        node = dict(op=op, a=a, b=b)
        if phase:node['phase']=phase
        if stored:
            node['site'] = sum('site' in n for n in ns)
        ns.append(node)
        return len(ns)-1

    def finish(ns, a, b):
        return dict(nodes=ns, outputs=[a,b], sites=sum('site' in n for n in ns))

    helpers, consumers = [], []
    cases = [('growth','plain'),('discount','plain'),('discount','reciprocal'),
             ('growth','expm1'),('discount','expm1')]
    for sign in ('growth','discount'):
        for unit in ('base2','natural'):
            for association in range(3):
                cases.append((sign,f'correction-{unit}-{association}'))
    for route in range(6):
        for mode, method in cases:
            ns = context(5); one = 5; phase=0
            def op(kind, a=0, b=0, stored=True):
                return append(ns,kind,a,b,stored,phase)
            exponent = op('neg',1,stored=False) if mode == 'discount' and method != 'reciprocal' else 1
            argument = op('add',one,0) if route < 3 else 0
            logop = 'fyl2x' if route < 3 else 'fyl2xp1'
            if route in (0,5):
                x = op(logop,exponent,argument)
            else:
                coefficient = op('fldln2',stored=False) if route in (2,4) else one
                logarithm = op(logop,coefficient,argument)
                x = op('mul',exponent,logarithm)
                if route in (2,4):
                    rom = op('fldl2e',stored=False)
                    x = op('mul',x,rom)
            phase=1
            integer = op('frndint',x,stored=False)
            fraction = op('sub',x,integer)
            fraction_power = op('f2xm1',fraction)
            if method == 'expm1':
                fractional = op('fscale',fraction_power,integer)
                integral = op('fscale',one,integer,stored=False)
                difference = op('sub',integral,one)
                em = op('add',fractional,difference)
                g = op('add',one,em)
            else:
                m = op('add',one,fraction_power)
                g = op('fscale',m,integer)
                if method == 'reciprocal':
                    g = op('div',one,g)
                em = op('sub',g,one)
                if method.startswith('correction-'):
                    phase=0
                    fields=method.split('-');unit,association=fields[1:]
                    coefficient=one if unit=='base2'else op('fldln2',stored=False)
                    denominator=op('fyl2x',coefficient,g)
                    numerator=x if unit=='base2'else op('mul',x,coefficient)
                    if association=='0':
                        a=op('mul',em,numerator);em=op('div',a,denominator)
                    elif association=='1':
                        a=op('div',numerator,denominator);em=op('mul',em,a)
                    else:
                        a=op('div',em,denominator);em=op('mul',a,numerator)
            helpers.append(dict(id=len(helpers),path=route,mode=mode,reconstruction=method,
                                graph=finish(ns,g,em)))
    # Derive the final consumer family independently from leaf partitions.
    from x87_consumer_grammar import templates as consumer_grammar
    consumers=consumer_grammar(context,append,finish)
    return helpers,consumers


def declared_profiles(templates):
    result=[]
    for t in templates:
        for cw in WORDS:
            inner_words=WORDS if any(n.get('phase')for n in t['graph']['nodes'])else(None,)
            for inner in inner_words:
                for mask in range(2**t['graph']['sites']):
                    row=dict(id=len(result),template=t['id'],cw=cw,mask=mask)
                    if inner is not None:row['inner_cw']=inner
                    result.append(row)
    return result


class Exact:
    def __init__(self,instructions):
        self.instructions=instructions
        self.counts=Counter()

    @lru_cache(maxsize=262144)
    def operation(self,kind,a,b,cw,store):
        x,y=decode(a),decode(b)
        if kind in ('fyl2x','fyl2xp1','f2xm1','fldln2','fldl2e'):
            if kind=='fyl2x':require(y>0,'FYL2X domain')
            if kind=='fyl2xp1':require(abs(y)<=Q(1,4),'FYL2XP1 domain')
            if kind=='f2xm1':require(abs(x)<=Q(1,2),'F2XM1 reduced domain')
            raw=self.instructions.run(kind,a,b,cw)
            if kind=='fldln2':require(decode(raw)==Q(0xb17217f7d1cf79ac,1<<64),'ROM ln2 differs')
            if kind=='fldl2e':require(decode(raw)==Q(0xb8aa3b295c17f0bc,1<<63),'ROM log2(e) differs')
        elif kind=='neg':
            raw=a[:9]+bytes([a[9]^128])
        elif kind=='frndint':
            lo=x.numerator//x.denominator; f=x-lo
            integer=lo+int(f>Q(1,2)or(f==Q(1,2)and lo%2))
            raw=encode(Q(integer),negative(a))
        elif kind=='fscale':
            k=abs(y.numerator)//y.denominator
            if y<0:k=-k
            require(-16382<k<16383,'scale outside certificate range')
            value=E.round(x*two(k))
            require(value is not None,'nonfinite FSCALE')
            raw=encode(value,negative(a))
        else:
            require(kind in ('add','sub','mul','div'),'unknown instruction')
            if kind=='add':
                z=x+y; sign=negative(a)and negative(b)
            elif kind=='sub':
                z=x-y; sign=negative(a)and not negative(b)
            elif kind=='mul':
                z=x*y; sign=negative(a)!=negative(b)
            else:
                require(y!=0,'zero denominator outside finite certificate scope')
                z=x/y; sign=negative(a)!=negative(b)
            fmt=P if cw==0x123f else E
            rounded=fmt.round(z)
            require(rounded is not None,'nonfinite arithmetic')
            raw=encode(rounded,sign if z==0 else z<0)
        self.counts[kind]+=1
        if store:
            value=B.round(decode(raw));require(value is not None,'nonfinite store')
            raw=encode(value,negative(raw))
            self.counts['binary64_stores']+=1
        return raw

    def graph(self,graph,inputs,cw,mask,extended=b'',inner_cw=None,trace=False):
        require(cw in WORDS and 0<=mask<1<<graph['sites'],'invalid control word/store mask')
        require(inner_cw is None or inner_cw in WORDS,'invalid inner control word')
        args=loaded_arguments(tuple(inputs),extended)
        values=[];raw=[]
        for i,n in enumerate(graph['nodes']):
            kind=n['op']
            if kind=='input':v=args[n['a']]
            elif kind=='const':v=encode(binary64(int(n['bits'],16)))
            else:
                require(n['a']<i and n['b']<i,'forward graph edge')
                word=inner_cw if n.get('phase',0)==1 else cw
                require(word in WORDS,'missing phase word')
                if trace:
                    raw.append(self.operation(kind,values[n['a']],values[n['b']],word,False))
                v=self.operation(kind,values[n['a']],values[n['b']],word,
                    'site'in n and bool(mask&(1<<n['site'])))
            if trace and kind in ('input','const'):raw.append(v)
            values.append(v)
        result=b''.join(values[o]for o in graph['outputs'])
        return (result,b''.join(raw),b''.join(values))if trace else result


def subtract_box(alive,helper_values,hp,hs,match_sets,box,touched):
    hsids=box['helpers'];failed=int(box['consumers'],16)
    require(hsids and len(hsids)==len(set(hsids))and failed>0,'empty/duplicate first kill box')
    helper=bytes.fromhex(box['helper'])
    removed=0
    for hid in hsids:
        require(0<=hid<len(hp)and hid not in touched,'unknown/duplicate helper at first kill')
        touched.add(hid)
        require(helper_values[hid]==helper,'shared helper value equality fails')
        key=hs[hp[hid]['template']]['mode'],helper
        require(key in match_sets,'missing consumer calculations')
        require(failed==alive[hid]&~match_sets[key],'first kill is incomplete or includes a match')
        require((failed&~alive[hid])==0,'assignment already killed')
        alive[hid]&=~failed
        removed+=failed.bit_count()
    return removed


def replay_row_tables(path, row, hs, cs, consumer_offsets, exact):
    """Verify one independent observation's complete arithmetic tables."""
    checked = Counter()
    before_operations = exact.counts.copy()
    before_native = exact.instructions.counts.copy()
    data = read(path)
    ins = [int(x, 16) for x in row['args']]
    require(data['witness_id'] == row['id'], 'witness ID differs')
    helper_values = []
    require(len(data['helper_tables']) == 4*len(hs), 'missing helper table')
    for ti, t in enumerate(data['helper_tables']):
        template = hs[ti//4]
        cw, inner = WORDS[(ti//2) % 2], WORDS[ti % 2]
        require((t['template'], t['cw'], t['inner_cw']) == (template['id'], cw, inner), 'helper table order differs')
        raw = base64.b64decode(t['values'], validate=True)
        require(len(raw) == 20*(1 << template['graph']['sites']), 'incomplete helper output table')
        for mask in range(1 << template['graph']['sites']):
            value = exact.graph(template['graph'], ins, cw, mask, inner_cw=inner)
            require(value == raw[20*mask:20*(mask+1)], 'helper instruction/store output differs')
            helper_values.append(value)
            checked['helper_executions'] += 1
    match_sets = {}
    for ev in data['consumer_evaluations']:
        key = ev['mode'], bytes.fromhex(ev['helper'])
        require(key not in match_sets and key[0] in ('growth', 'discount'), 'duplicate/unknown consumer input')
        relevant = [c for c in cs if c['mode'] == key[0]]
        require(len(ev['tables']) == 4*len(relevant), 'missing consumer table')
        matches = 0
        for ti, t in enumerate(ev['tables']):
            c = relevant[ti//4]
            cw, inner = WORDS[(ti//2) % 2], WORDS[ti % 2]
            require((t['template'], t['cw'], t['inner_cw']) == (c['id'], cw, inner), 'consumer table order differs')
            raw = base64.b64decode(t['values'], validate=True)
            require(len(raw) == 10*(1 << c['graph']['sites']), 'incomplete consumer output table')
            offset = consumer_offsets[c['id'], cw, inner]
            for mask in range(1 << c['graph']['sites']):
                value = exact.graph(c['graph'], ins, cw, mask, key[1], inner_cw=inner)[:10]
                require(value == raw[10*mask:10*(mask+1)], 'consumer arithmetic/store output differs')
                if published_bits(value) == int(row['expected'], 16):
                    matches |= 1 << (offset+mask)
                checked['consumer_executions'] += 1
        match_sets[key] = matches
    return dict(helper_values=helper_values, match_sets=match_sets, checked=dict(checked),
        operation_misses=dict(exact.counts-before_operations),
        native_replays=dict(exact.instructions.counts-before_native),
        leaves=list(exact.instructions.transcript.items()))


def initialize_row_worker(hs, cs, consumer_offsets, library):
    global ROW_CONTEXT
    ROW_CONTEXT = hs, cs, consumer_offsets, Exact(Instructions(library))


def replay_row_worker(task):
    path, row = task
    return replay_row_tables(path, row, *ROW_CONTEXT)


def row_replays(files, rows, hs, cs, consumer_offsets, library, exact, workers):
    """Only table arithmetic is independent; first-kill coverage stays ordered."""
    tasks = [(path, rows[ri]) for ri, path in enumerate(files)]
    if workers == 1:
        for path, row in tasks:
            yield replay_row_tables(path, row, hs, cs, consumer_offsets, exact)
    else:
        import multiprocessing
        # Spawn avoids inheriting the parent's large coverage bitsets/caches.
        with multiprocessing.get_context('spawn').Pool(workers, initializer=initialize_row_worker,
                initargs=(hs, cs, consumer_offsets, str(library))) as pool:
            yield from pool.imap(replay_row_worker, tasks, chunksize=1)


def check(directory,trace,library,output,leaves=None,workers=1):
    require(type(workers) is int and 1 <= workers <= 8, 'workers must be between1 and8')
    start=time.monotonic();directory=Path(directory)
    source_paths=[Path(__file__),ROOT/'checker/exact.py',ROOT/'checker/x87_isa.c',ROOT/'checker/x87_consumer_grammar.py']
    source_sha={str(p):digest(p)for p in source_paths}
    inputs=list(directory.glob('*.json'))+list((directory/'rows').glob('*.json'))+[directory/'witnesses.jsonl',directory/'first-kills.jsonl',directory/'first-kill-traces.jsonl',Path(trace),Path(library)]
    frozen={str(p):digest(p)for p in inputs}
    hs,cs=canonical_templates()
    require(read(directory/'helpers.json')==hs,'helper grammar differs')
    require(read(directory/'consumers.json')==cs,'consumer grammar differs')
    hp,cp=declared_profiles(hs),declared_profiles(cs)
    require(read(directory/'helper-profiles.json')==hp,'helper profiles incomplete or changed')
    require(read(directory/'consumer-profiles.json')==cp,'consumer profiles incomplete or changed')
    original=records(trace);rows=[];excluded=[]
    for row in original:
        require(row['function']=='PMT'and len(row['args'])==5,'not an original PMT identity')
        require(all(isinstance(x,str)and x.startswith('0x')and len(x)==18 for x in row['args']),'inexact argument encoding')
        bits=[int(x,16)for x in row['args']]
        for x in bits:binary64(x)
        require(row.get('result_index')is None,'selected result outside scope')
        reason=None
        if not 0<bits[0]<=0x3fd0000000000000:reason='outside positive FYL2XP1 domain'
        elif bits[4]not in (0,0x3ff0000000000000):reason='noncanonical timing'
        elif not row['expected'].startswith('0x'):reason='nonnumeric target'
        if reason:excluded.append(dict(reason=reason,witness=row))
        else:
            require(binary64(bits[1])>0,'positive periods required')
            require(binary64(row['expected'])!=0,'nonzero target required')
            rows.append(row)
    require(records(directory/'witnesses.jsonl')==rows,'original witness content/order differs')
    require(read(directory/'excluded.json')==excluded,'scope exclusions differ')
    cbits={mode:sum(1<<c['id']for c in cp if cs[c['template']]['mode']==mode)for mode in ('growth','discount')}
    alive=[cbits[hs[h['template']]['mode']]for h in hp]
    helper_modes=Counter(hs[h['template']]['mode']for h in hp)
    total=sum(helper_modes[mode]*bits.bit_count()for mode,bits in cbits.items())
    consumer_offsets={(p['template'],p['cw'],p['inner_cw']):p['id']for p in cp if p['mask']==0}
    kills=records(directory/'first-kills.jsonl');byrow={}
    for k in kills:byrow.setdefault(k['row'],[]).append(k)
    transcripts=records(directory/'first-kill-traces.jsonl')
    require(len(transcripts)==len(kills),'missing representative first-kill trace')
    trace_index=0
    exact=Exact(Instructions(library));checked=Counter();used=0;excluded_count=0
    files=sorted((directory/'rows').glob('*.json'))
    require([p.name for p in files]==[f'{i:03d}.json'for i in range(len(files))],'missing/unexpected row certificate')
    require(len(files) <= len(rows), 'unknown witness row')
    tables = row_replays(files, rows, hs, cs, consumer_offsets, library, exact, workers)
    worker_operations, worker_native = Counter(), Counter()
    try:
        for ri, replayed in enumerate(tables):
            row = rows[ri]
            ins = [int(x, 16) for x in row['args']]
            helper_values, match_sets = replayed['helper_values'], replayed['match_sets']
            checked.update(replayed['checked'])
            worker_operations.update(replayed['operation_misses'])
            worker_native.update(replayed['native_replays'])
            for key, value in replayed['leaves']:
                prior = exact.instructions.transcript.get(key)
                require(prior is None or prior == value, 'native leaves differ between row workers')
                exact.instructions.transcript[key] = value
            touched=set()
            for box in byrow.get(ri,[]):
                excluded_count+=subtract_box(alive,helper_values,hp,hs,match_sets,box,touched)
                tr=transcripts[trace_index];hid=box['helpers'][0]
                failed=int(box['consumers'],16);cid=(failed&-failed).bit_length()-1
                require((tr['row'],tr['box'],tr['helper_profile'],tr['consumer_profile'])==(ri,trace_index,hid,cid),'first-kill trace identity differs')
                h,c=hp[hid],cp[cid]
                ho,hr,hv=exact.graph(hs[h['template']]['graph'],ins,h['cw'],h['mask'],inner_cw=h['inner_cw'],trace=True)
                co,cr,cv=exact.graph(cs[c['template']]['graph'],ins,c['cw'],c['mask'],ho,inner_cw=c['inner_cw'],trace=True)
                require((ho.hex(),hr,hv,co.hex(),cr,cv)==(tr['helper_output'],base64.b64decode(tr['helper_raw']),base64.b64decode(tr['helper_values']),tr['consumer_output'],base64.b64decode(tr['consumer_raw']),base64.b64decode(tr['consumer_values'])),'full-width first-kill node trace differs')
                trace_index+=1
                checked['first_kill_boxes']+=1
            # Every original first contradiction is still applied in order.
            for hid,bits in enumerate(alive):
                if bits:
                    key=hs[hp[hid]['template']]['mode'],helper_values[hid]
                    require(key in match_sets and bits&~match_sets[key]==0,'a first contradiction was omitted')
            checked['rows']+=1;used+=1
            print(json.dumps(dict(row=ri,id=row['id'],checked=dict(checked),seconds=time.monotonic()-start)),flush=True)
    finally:
        tables.close()
    declared_survivors=[dict(helper=i,consumers=hex(bits))for i,bits in enumerate(alive)if bits]
    require(read(directory/'survivors.json')==declared_survivors,'final survivor set differs')
    remaining=sum(x.bit_count()for x in alive)
    require(total==excluded_count+remaining,'Cartesian coverage count differs')
    require(trace_index==len(kills),'kill record for unknown row')
    require(remaining==0 or used==len(rows),'surviving assignments did not see all scoped witnesses')
    summary=read(directory/'summary.json')
    for field,value in dict(complete=True,helper_profiles=len(hp),consumer_profiles=len(cp),
            declared_assignments=total,excluded_assignments=excluded_count,surviving_assignments=remaining,
            first_kill_boxes=checked['first_kill_boxes'],rows_used=used,scoped_rows=len(rows),excluded_rows=len(excluded),trace_sha256=digest(trace)).items():
        require(summary[field]==value,'summary field differs: '+field)
    for path,expected_digest in summary['source_sha256'].items():require(digest(ROOT/path)==expected_digest,'candidate source differs')  # ubs:ignore - Public file integrity hashes, never secrets or authentication.
    require(source_sha=={str(p):digest(p)for p in source_paths},'checker sources changed')
    require(frozen=={str(p):digest(p)for p in inputs},'proof inputs changed')
    receipt=dict(complete=True,declared_assignments=total,excluded_assignments=excluded_count,
        surviving_assignments=remaining,checked=dict(checked),workers=workers,
        unique_exact_operations=dict(exact.counts) if workers == 1 else None,
        worker_operation_cache_misses=dict(worker_operations),
        native_microcode_replays=dict(exact.instructions.counts if workers == 1 else worker_native+exact.instructions.counts),seconds=time.monotonic()-start,
        inputs_sha256=frozen,checker_sha256=source_sha,
        claim='Complete coverage of the declared four-phase PC53/PC64 store-mask family. Every first kill is reproduced by independent rational arithmetic and direct native instruction leaves.',
        limits='Native AMD microcode replay is conditional on this CPU; not mathematical correctly rounded transcendentals or a software x87 emulation. Scope/input provenance and dead/redundant labels are retained. No fitted helper, identified Excel implementation or corpus gain.')
    if leaves:
        leaves=Path(leaves);leaves.mkdir(parents=True,exist_ok=False)
        with (leaves/'native-leaf-requests.tsv').open('x')as requests,(leaves/'native-leaf-reference.tsv').open('x')as reference:
            requests.write('id\tkind\tcw\ta\tb\n');reference.write('id\tkind\tcw\ta\tb\tresult\n')
            for i,((kind,a,b,cw),result)in enumerate(sorted(exact.instructions.transcript.items())):
                row=f'leaf-{i:05d}\t{kind}\t{cw:04x}\t{a.hex()}\t{b.hex()}'
                requests.write(row+'\n');reference.write(row+'\t'+result.hex()+'\n')
        receipt['native_leaf_requests']=len(exact.instructions.transcript)
        receipt['native_leaf_sha256']={p.name:digest(p)for p in leaves.iterdir()}
    Path(output).write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for a in ('directory','trace','library','output'):p.add_argument(a)
    p.add_argument('--leaves')
    p.add_argument('--workers', type=int, default=1)
    a=p.parse_args();check(a.directory,a.trace,a.library,a.output,a.leaves,a.workers)


if __name__=='__main__':main()
