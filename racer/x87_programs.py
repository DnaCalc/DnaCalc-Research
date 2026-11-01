"""Finite x87 instruction programs, exact native evaluation, and box certificates.

Original ISA/algebra implementation, MIT. No libm, coefficients, fitted helpers,
or selected answers. Program stores are explicit mask bits, never Python floats.
"""
import argparse
import base64
from collections import Counter, defaultdict
import ctypes as C
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

OPS = ('input', 'const', 'add', 'sub', 'mul', 'div', 'neg', 'fyl2x',
       'fyl2xp1', 'f2xm1', 'frndint', 'fscale', 'fldln2', 'fldl2e')
WORDS = (0x123f, 0x133f)
ROOT = Path(__file__).resolve().parent.parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2) + '\n')


class Ext(C.Structure):
    _fields_ = [('b', C.c_ubyte * 10)]


class Node(C.Structure):
    _fields_ = [('op', C.c_int32), ('a', C.c_int32), ('b', C.c_int32),
                ('site', C.c_int32), ('phase', C.c_int32), ('bits', C.c_uint64)]


class Native:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.lib = C.CDLL(str(self.path))
        extp, nodep = C.POINTER(Ext), C.POINTER(Node)
        self.lib.x87_eval.argtypes = [nodep, C.c_int, C.POINTER(C.c_uint64), C.c_int,
            extp, C.c_int, C.c_uint32, C.c_uint64, C.c_int, C.c_int, extp, extp, extp]
        self.lib.x87_masks.argtypes = self.lib.x87_eval.argtypes[:8] + [
            C.c_uint64, C.c_int, C.c_int, extp]
        self.lib.x87_leaf.argtypes = [C.c_int, extp, extp, C.c_uint16, extp]
        self.lib.x87_leaf.restype = None
        self.lib.x87_publish.argtypes = [extp, C.c_uint16]
        self.lib.x87_publish.restype = C.c_uint64

    @staticmethod
    @lru_cache(maxsize=1024)
    def nodes(encoded):
        g = json.loads(encoded)
        return (Node * len(g['nodes']))(*(Node(OPS.index(n['op']), n.get('a', 0),
            n.get('b', 0), n.get('site', -1), n.get('phase',0), int(n.get('bits', '0'), 16))
            for n in g['nodes']))

    def evaluate(self, graph, inputs, cw, mask=0, extended=b'', trace=False, inner_cw=None):
        if len(extended)%10:raise ValueError('incomplete extended operand')
        nodes = self.nodes(json.dumps(graph, sort_keys=True))
        ins = (C.c_uint64 * len(inputs))(*inputs)
        exts = (Ext * (len(extended)//10)).from_buffer_copy(extended)
        outputs = (Ext * 2)()
        raw, values = (Ext * len(nodes))(), (Ext * len(nodes))()
        words=cw|((cw if inner_cw is None else inner_cw)<<16)
        rc = self.lib.x87_eval(nodes, len(nodes), ins, len(ins), exts, len(exts), words,
            mask, *graph['outputs'], outputs, raw if trace else None,
            values if trace else None)
        if rc:
            raise ValueError('malformed instruction graph')
        return (bytes(outputs), bytes(raw), bytes(values)) if trace else bytes(outputs)

    def masks(self, graph, inputs, cw, extended=b'', inner_cw=None):
        if len(extended)%10:raise ValueError('incomplete extended operand')
        nodes = self.nodes(json.dumps(graph, sort_keys=True))
        ins = (C.c_uint64 * len(inputs))(*inputs)
        exts = (Ext * (len(extended)//10)).from_buffer_copy(extended)
        count = 1 << graph['sites']
        outputs = (Ext * (2*count))()
        words=cw|((cw if inner_cw is None else inner_cw)<<16)
        rc = self.lib.x87_masks(nodes, len(nodes), ins, len(ins), exts, len(exts), words,
            0, count, *graph['outputs'], outputs)
        if rc:
            raise ValueError('malformed instruction graph')
        return bytes(outputs)

    def leaf(self, op, a, b, cw):
        aa, bb = Ext.from_buffer_copy(a), Ext.from_buffer_copy(b)
        out = Ext()
        self.lib.x87_leaf(OPS.index(op), C.byref(aa), C.byref(bb), cw, C.byref(out))
        return bytes(out)

    def publish(self, x, cw=0x133f):
        a = Ext.from_buffer_copy(x)
        return self.lib.x87_publish(C.byref(a), cw)


class Builder:
    def __init__(self, inputs):
        self.nodes = [dict(op='input', a=i) for i in range(inputs)]
        self.sites = 0
        self.phase = 0
        self.one = self.constant('3ff0000000000000')

    def constant(self, bits):
        self.nodes.append(dict(op='const', bits=bits))
        return len(self.nodes)-1

    def op(self, op, a=0, b=0, store=True):
        n = dict(op=op, a=a, b=b)
        if self.phase:n['phase']=self.phase
        if store:
            n['site'] = self.sites
            self.sites += 1
        self.nodes.append(n)
        return len(self.nodes)-1

    def graph(self, *outputs):
        return dict(nodes=self.nodes, outputs=list(outputs), sites=self.sites)


def helper_templates():
    templates = []
    # Every branch here is a fixed program description, independent of outputs.
    for path in range(6):
        cases = [('growth', 'plain'), ('discount', 'plain'),
                 ('discount', 'reciprocal'), ('growth', 'expm1'), ('discount', 'expm1')]
        cases += [(mode, f'correction-{units}-{association}')
                  for mode in ('growth', 'discount') for units in ('base2', 'natural')
                  for association in range(3)]
        for mode, reconstruction in cases:
            g = Builder(5)
            r, n = 0, 1
            signed_n = g.op('neg', n, store=False) if mode == 'discount' and reconstruction != 'reciprocal' else n
            if path <= 2:
                base = g.op('add', g.one, r)
                argument, instruction = base, 'fyl2x'
            else:
                argument, instruction = r, 'fyl2xp1'
            if path in (0, 5):
                x = g.op(instruction, signed_n, argument)
            else:
                natural = path in (2, 4)
                coefficient = g.op('fldln2', store=False) if natural else g.one
                log = g.op(instruction, coefficient, argument)
                x = g.op('mul', signed_n, log)
                if natural:
                    rom = g.op('fldl2e', store=False)
                    x = g.op('mul', x, rom)
            g.phase = 1
            whole = g.op('frndint', x, store=False)
            fraction = g.op('sub', x, whole)
            small = g.op('f2xm1', fraction)
            if reconstruction == 'expm1':
                scaled = g.op('fscale', small, whole)
                power = g.op('fscale', g.one, whole, store=False)
                integer = g.op('sub', power, g.one)
                em = g.op('add', scaled, integer)
                growth = g.op('add', g.one, em)
            else:
                mantissa = g.op('add', g.one, small)
                growth = g.op('fscale', mantissa, whole)
                if reconstruction == 'reciprocal':
                    growth = g.op('div', g.one, growth)
                em = g.op('sub', growth, g.one)
                if reconstruction.startswith('correction-'):
                    g.phase = 0
                    _, units, association = reconstruction.split('-')
                    coefficient = g.op('fldln2', store=False) if units == 'natural' else g.one
                    logarithm = g.op('fyl2x', coefficient, growth)
                    t = g.op('mul', x, coefficient) if units == 'natural' else x
                    if association == '0':
                        em = g.op('div', g.op('mul', em, t), logarithm)
                    elif association == '1':
                        em = g.op('mul', em, g.op('div', t, logarithm))
                    else:
                        em = g.op('mul', g.op('div', em, logarithm), t)
            templates.append(dict(id=len(templates), path=path, mode=mode,
                reconstruction=reconstruction, graph=g.graph(growth, em)))
    return templates


def consumer_templates():
    from x87_consumers import templates
    return templates(Builder)


def profiles(templates):
    result=[]
    for t in templates:
        mixed=any(n.get('phase')for n in t['graph']['nodes'])
        for cw in WORDS:
            for inner in (WORDS if mixed else (None,)):
                for mask in range(1<<t['graph']['sites']):
                    p=dict(id=len(result),template=t['id'],cw=cw,mask=mask)
                    if inner is not None:p['inner_cw']=inner
                    result.append(p)
    return result


def build_library(output):
    output = Path(output).resolve()
    command = ['cc', '-std=c11', '-O2', '-Wall', '-Wextra', '-Werror', '-shared',
               '-fPIC', str(ROOT/'racer/x87_native.c'), '-o', str(output)]
    subprocess.run(command, check=True, timeout=120)
    return command


def b64(data):
    return base64.b64encode(data).decode('ascii')


def run(trace, output, library):
    start = time.monotonic()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    (output/'rows').mkdir()
    native = Native(library)
    sources = {str(p.relative_to(ROOT)): digest(p) for p in
        (ROOT/'racer/x87_programs.py', ROOT/'racer/x87_native.c', ROOT/'racer/x87_consumers.py', ROOT/'x87-r13-unit.json')}
    original = [json.loads(l) for l in Path(trace).read_text().splitlines()]
    rows, excluded = [], []
    for row in original:
        args = row['args']
        bits = [int(x, 16) for x in args]
        # Positive finite rate <=1/4 is safely inside FYL2XP1's domain.
        reason = None
        if len(bits) != 5 or any(not x.startswith('0x') or len(x) != 18 for x in args):
            raise ValueError('original exact numeric arguments required')
        if not 0 < bits[0] <= 0x3fd0000000000000:
            reason = 'outside positive FYL2XP1 domain'
        elif bits[4] not in (0, 0x3ff0000000000000):
            reason = 'noncanonical timing'
        elif not row['expected'].startswith('0x'):
            reason = 'nonnumeric target'
        if reason:
            excluded.append(dict(reason=reason, witness=row))
        else:
            rows.append(row)
    write(output/'excluded.json', excluded)
    (output/'witnesses.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    hs, cs = helper_templates(), consumer_templates()
    hp, cp = profiles(hs), profiles(cs)
    write(output/'helpers.json', hs)
    write(output/'consumers.json', cs)
    write(output/'helper-profiles.json', hp)
    write(output/'consumer-profiles.json', cp)
    # Growth helpers compose only with growth consumers; likewise discount.
    cbits = {mode: sum(1 << p['id'] for p in cp if cs[p['template']]['mode'] == mode)
             for mode in ('growth', 'discount')}
    survivors = [cbits[hs[h['template']]['mode']] for h in hp]
    helper_modes=Counter(hs[h['template']]['mode']for h in hp)
    total = sum(helper_modes[mode]*bits.bit_count()for mode,bits in cbits.items())
    killed_total, boxes, used = 0, 0, 0
    consumer_offsets, offset = {}, 0
    for c in cs:
        for cw in WORDS:
            for inner in WORDS:
                consumer_offsets[c['id'], cw, inner] = offset
                offset += 1 << c['graph']['sites']
    kills = (output/'first-kills.jsonl').open('x')
    traces = (output/'first-kill-traces.jsonl').open('x')
    for row_index, row in enumerate(rows):
        inputs = [int(x, 16) for x in row['args']]
        helper_tables, helper_values = [], []
        for h in hs:
            for cw in WORDS:
                for inner in WORDS:
                    values = native.masks(h['graph'], inputs, cw, inner_cw=inner)
                    helper_tables.append(dict(template=h['id'], cw=cw, inner_cw=inner, values=b64(values)))
                    helper_values.extend(values[i:i+20] for i in range(0, len(values), 20))
        groups = defaultdict(list)
        for hid, alive in enumerate(survivors):
            if alive:
                mode = hs[hp[hid]['template']]['mode']
                groups[mode, helper_values[hid], alive].append(hid)
        evaluations, computed, kills_here, remaining = [], {}, 0, 0
        for (mode, value, alive), helpers in groups.items():
            key = mode, value
            if key not in computed:
                matches, tables = 0, []
                for c in cs:
                    if c['mode'] != mode:
                        continue
                    for cw in WORDS:
                        for inner in WORDS:
                            values = native.masks(c['graph'], inputs, cw, value,inner_cw=inner)
                            raw_outputs = b''.join(values[i:i+10] for i in range(0, len(values), 20))
                            tables.append(dict(template=c['id'], cw=cw, inner_cw=inner,values=b64(raw_outputs)))
                            offset = consumer_offsets[c['id'], cw, inner]
                            for m in range(len(raw_outputs)//10):
                                published = native.publish(raw_outputs[10*m:10*(m+1)], cw)
                                if published == int(row['expected'], 16):
                                    matches |= 1 << (offset+m)
                computed[key] = matches
                evaluations.append(dict(mode=mode, helper=value.hex(), tables=tables))
            failed = alive & ~computed[key]
            after = alive & computed[key]
            remaining += len(helpers)*after.bit_count()
            for hid in helpers:
                survivors[hid] = after
            if failed:
                kills.write(json.dumps(dict(row=row_index, helpers=helpers,
                    helper=value.hex(), consumers=hex(failed)))+'\n')
                hid=helpers[0];cid=(failed&-failed).bit_length()-1
                h=hp[hid];c=cp[cid]
                ho,hr,hv=native.evaluate(hs[h['template']]['graph'],inputs,h['cw'],h['mask'],
                                        trace=True,inner_cw=h['inner_cw'])
                co,cr,cv=native.evaluate(cs[c['template']]['graph'],inputs,c['cw'],c['mask'],value,trace=True,inner_cw=c['inner_cw'])
                traces.write(json.dumps(dict(row=row_index,box=boxes,helper_profile=hid,consumer_profile=cid,
                    helper_output=ho.hex(),helper_raw=b64(hr),helper_values=b64(hv),
                    consumer_output=co.hex(),consumer_raw=b64(cr),consumer_values=b64(cv)))+'\n')
                boxes += 1
                kills_here += len(helpers)*failed.bit_count()
        write(output/'rows'/f'{row_index:03d}.json', dict(witness_id=row['id'],
            helper_tables=helper_tables, consumer_evaluations=evaluations))
        killed_total += kills_here
        used += 1
        print(json.dumps(dict(row=row_index, id=row['id'], killed=kills_here,
            remaining=remaining, groups=len(groups), seconds=time.monotonic()-start)), flush=True)
        if killed_total+remaining != total:
            raise ValueError('coverage arithmetic differs')
        if not remaining:
            break
    kills.close()
    traces.close()
    write(output/'survivors.json', [dict(helper=i, consumers=hex(x))
                                   for i, x in enumerate(survivors) if x])
    if any(digest(ROOT/p) != expected_digest for p, expected_digest in sources.items()):  # ubs:ignore - Public file integrity hashes, never secrets or authentication.
        raise ValueError('research sources changed during search')
    write(output/'summary.json', dict(complete=True, helper_profiles=len(hp),
        consumer_profiles=len(cp), declared_assignments=total, excluded_assignments=killed_total,
        surviving_assignments=total-killed_total, first_kill_boxes=boxes,
        rows_used=used, scoped_rows=len(rows), excluded_rows=len(excluded),
        trace_sha256=digest(trace), library_sha256=digest(library), source_sha256=sources,
        seconds=time.monotonic()-start,
        limits='Finite instruction grammar, native CPU microcode and original recorded-input assumptions; no Excel kernel identified, no fitted intermediates, no product gains. Counts include dead or redundant physical store labels.'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    b = sub.add_parser('build'); b.add_argument('library')
    r = sub.add_parser('run')
    r.add_argument('trace'); r.add_argument('output'); r.add_argument('library')
    a = p.parse_args()
    if a.command == 'build':
        print(json.dumps(build_library(a.library)))
    else:
        run(a.trace, a.output, a.library)


if __name__ == '__main__':
    main()
