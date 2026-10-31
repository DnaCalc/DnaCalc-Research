"""R8 structure discovery; measurements, never a calculation implementation.

Implementer, tier 1 annuity algebra/IEEE 754 and original Excel observations.
Consumes existing xtask witness traces, preserving typed identities and sources.
Groups using inferred inverse inputs are labelled separately from direct inputs.
"""
import collections as C
import hashlib
import itertools
import json
import math
from pathlib import Path
import struct
import sys

FUNCTIONS = ['PV', 'FV', 'NPER', 'PMT', 'IPMT', 'PPMT', 'RATE']


def lines(path):
    with open(path) as f:
        for line in f:
            yield json.loads(line)


def num(x):
    if isinstance(x, str) and x.startswith('0x') and len(x) == 18:
        return struct.unpack('>d', bytes.fromhex(x[2:]))[0]
    return None


def bits(x):
    return '0x' + struct.pack('>d', x).hex()


def argtype(x):
    if isinstance(x, bool):
        return 'boolean'
    if num(x) is not None:
        return 'number'
    if x is None:
        return 'omitted'
    return type(x).__name__ + ':' + str(x)[:40]


def coordinate(row):
    """NPER/RATE's output is a projection, never a direct input pair."""
    f, a = row['function'], row['args']
    if f in ['PV', 'FV', 'PMT'] and len(a) >= 2:
        return a[0], a[1], 'direct'
    if f in ['IPMT', 'PPMT'] and len(a) >= 3:
        return a[0], a[2], 'direct'
    if f == 'NPER' and len(a) >= 1:
        return a[0], row['expected'], 'inverse_output_projection'
    if f == 'RATE' and len(a) >= 1:
        return row['expected'], a[0], 'inverse_output_projection'
    return None


def region(r, n):
    if r is None or n is None or not math.isfinite(r) or not math.isfinite(n):
        return {'rate_region': 'unavailable', 'nper_region': 'unavailable', 'tau_region': 'unavailable'}
    rr = ('zero' if r == 0 else 'base_nonpositive' if r <= -1 else
          'rounded_base_one' if 1 + r == 1 else 'abs_rate_below_1e-12' if abs(r) < 1e-12 else
          'abs_rate_below_2^-20' if abs(r) < 2**-20 else 'ordinary_positive' if r > 0 else 'ordinary_negative')
    nr = 'one' if n == 1 else 'zero' if n == 0 else 'positive_integer' if n > 0 and n.is_integer() else 'negative_integer' if n.is_integer() else 'fractional'
    tau = n * math.log1p(r) if r > -1 else math.nan
    tr = 'unavailable' if not math.isfinite(tau) else 'zero' if tau == 0 else ('discount_' if tau > 0 else 'growth_') + ('abs_tau_lt_1' if abs(tau) < 1 else 'abs_tau_1_to_700' if abs(tau) < 700 else 'abs_tau_ge_700')
    return dict(rate_region=rr, nper_region=nr, tau_region=tr)


def gap_clusters(round6, round7, out):
    summaries = {}
    for name, root, frozen, expected in [
        ('r6_frozen_expm1', round6/'native-final', round6/'frozen/scoped.jsonl', 38432),
        ('r7_helper_family', round7/'wide-native', round7/'frozen/scoped.jsonl', 11932),
    ]:
        gaps = set()
        for path in sorted(root.glob('unreachable-roots-*.jsonl')):
            gaps.update(r['trace_index'] for r in lines(path))
        assert len(gaps) == expected, (name, len(gaps))
        counts = C.defaultdict(lambda: C.Counter())
        sample = []
        for r in lines(frozen):
            raw = r['raw']
            hit = r['trace_index'] in gaps
            # The two families have different scopes. Every rate below uses its
            # own family's denominator, never the other family's population.
            if name.startswith('r7') and num(raw['args'][3]) != 0:
                continue
            c = coordinate(raw)
            tags = dict(build=raw['reading'].get('excel_build'), source=raw['source'],
                        types=','.join(map(argtype, raw['args'])),
                        timing='zero' if num(raw['args'][4]) == 0 else 'nonzero')
            tags.update(region(num(c[0]), num(c[1])))
            for axis, val in tags.items():
                counts[(axis, str(val))]['rows'] += 1
                counts[(axis, str(val))]['gaps'] += hit
            if hit and len(sample) < 12:
                sample.append(dict(trace_index=r['trace_index'], witness=raw))
        table = [dict(axis=a, value=v, **c, gap_rate=c['gaps']/c['rows']) for (a, v), c in sorted(counts.items())]
        summaries[name] = dict(gaps=len(gaps), clusters=table, examples=sample,
            limitation='Conditioned on original numeric finite-discount scope; typed or omitted arguments outside that scope cannot explain these gaps. Build/source/region confounding is retained.')
    (out/'gap-clusters.json').write_text(json.dumps(summaries, indent=2)+'\n')
    return {name: data['gaps'] for name, data in summaries.items()}


def discover(trace_dir, out):
    out.mkdir(parents=True, exist_ok=True)
    # Compact per-function outcome counts within exact rate/nper/build cells.
    groups = C.defaultdict(lambda: C.defaultdict(lambda: C.Counter()))
    fullenv = C.defaultdict(set)
    totals, regions, errors = {}, C.defaultdict(C.Counter), C.defaultdict(C.Counter)
    cases = []
    boundaries = C.defaultdict(lambda: C.defaultdict(C.Counter))
    for fn in FUNCTIONS:
        path = trace_dir / (fn + '.jsonl')
        if not path.exists():
            continue
        total = C.Counter()
        for idx, row in enumerate(lines(path)):
            assert row['function'] == fn
            total['rows'] += 1
            status = 'matched' if row['expected'] == row['actual'] else 'refused' if row['actual'].startswith('fault:') else 'mismatched'
            total[status] += 1
            types = ','.join(map(argtype, row['args']))
            errors[(fn, 'types', types)][status] += 1
            errors[(fn, 'build', str(row['reading'].get('excel_build')))][status] += 1
            c = coordinate(row)
            if not c:
                continue
            r, n, projection = num(c[0]), num(c[1]), c[2]
            tags = region(r, n)
            for a, v in tags.items():
                regions[(fn, a, v)][status] += 1
            errors[(fn, 'expected_kind', 'number' if num(row['expected']) is not None else row['expected'])][status] += 1
            if r is None or n is None or not math.isfinite(r) or not math.isfinite(n):
                total['unavailable_coordinates'] += 1
                continue
            build = str(row['reading'].get('excel_build'))
            gkey = (c[0], c[1], build, projection)
            groups[gkey][fn][status] += 1
            fullenv[gkey].add(json.dumps(row['reading'], sort_keys=True))
            # Record exact observed kind transitions. Adjacent observed rates do
            # not establish an exact branch boundary unless their bits are adjacent.
            if projection == 'direct':
                # A rate transition is meaningful only with every other input
                # and the complete recorded environment fixed.
                key = (fn, c[1], build, projection, json.dumps(row['args'][1:], sort_keys=True),
                       json.dumps(row['reading'], sort_keys=True))
                boundaries[key][c[0]][row['expected'] if num(row['expected']) is None else 'number'] += 1
            if projection == 'direct' and fn in ['PV', 'FV'] and len(row['args']) >= 4:
                a = [num(x) for x in row['args']]
                if all(x is not None for x in a) and a[2] == 0 and a[3] != 0 and num(row['expected']) is not None:
                    cases.append(dict(trace_index=idx, witness=row, isolation='growth_only' if fn=='FV' else 'inverse_growth_only'))
        total['fail'] = total['rows'] - total['matched']
        total['pass_rate_decoded'] = total['matched']/total['rows'] if total['rows'] else 0
        total['trace_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        totals[fn] = dict(total)
    pair_counts = C.defaultdict(C.Counter)
    examples = C.defaultdict(list)
    with (out/'joint-groups.jsonl').open('x') as f:
        for (rb, nb, build, projection), funcs in sorted(groups.items()):
            if len(funcs) < 2:
                continue
            data = dict(rate=rb, nper=nb, build=build, projection=projection,
                        functions={k:dict(v) for k,v in sorted(funcs.items())},
                        distinct_full_environments=len(fullenv[(rb,nb,build,projection)]))
            f.write(json.dumps(data,separators=(',',':'))+'\n')
            for a,b in itertools.combinations(sorted(funcs),2):
                sa = 'has_error' if funcs[a]['mismatched']+funcs[a]['refused'] else 'all_matched'
                sb = 'has_error' if funcs[b]['mismatched']+funcs[b]['refused'] else 'all_matched'
                pair_counts[(a,b,projection)][sa+' / '+sb] += 1
                if len(examples[(a,b,sa,sb,projection)]) < 3:
                    examples[(a,b,sa,sb,projection)].append(data)
    transitions=[]
    for (fn,nb,build,projection,other,environment), rates in boundaries.items():
        vals=sorted(rates,key=lambda x:num(x))
        for a,b in zip(vals,vals[1:]):
            ka,kb=sorted(rates[a]),sorted(rates[b])
            if ka!=kb:
                transitions.append(dict(function=fn,nper=nb,build=build,projection=projection,
                    fixed_other_arguments=json.loads(other),environment=json.loads(environment),
                    lower=a,upper=b,lower_kinds=ka,upper_kinds=kb,
                    adjacent_binary64=math.nextafter(num(a),math.inf)==num(b)))
    (out/'growth-only.jsonl').write_text(''.join(json.dumps(c,separators=(',',':'))+'\n'for c in cases))
    summary=dict(totals=totals,
        pair_cooccurrence=[dict(functions=[a,b],projection=p,group_counts=dict(c))for(a,b,p),c in sorted(pair_counts.items())],
        pair_examples=[dict(functions=[a,b],states=[sa,sb],projection=p,groups=x)for(a,b,sa,sb,p),x in sorted(examples.items())],
        regions=[dict(function=f,axis=a,value=v,counts=dict(c))for(f,a,v),c in sorted(regions.items())],
        typed_build_and_kind=[dict(function=f,axis=a,value=v,counts=dict(c))for(f,a,v),c in sorted(errors.items())],
        growth_isolation_rows=len(cases),observed_kind_transitions=transitions,
        limitations=['Counts cover decoded traces; unread failures are in original judge logs.',
            'Co-occurrence uses exact input bits and build. Full environment differences are counted; build-only groups cannot prove the same machine.',
            'NPER and RATE coordinate projections use recorded outputs, are separate from direct input groups, and do not expose RATE solver evaluation rates.',
            'Observed adjacent sample transitions are brackets, not identified internal branch points.'])
    (out/'structure-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(dict(totals=totals,growth_isolation_rows=len(cases),joint_groups=sum(1 for x in groups.values() if len(x)>1),kind_transitions=len(transitions))))


if __name__ == '__main__':
    out=Path(sys.argv[-1]);out.mkdir(exist_ok=True,parents=True)
    if sys.argv[1]=='gaps':
        print(json.dumps(gap_clusters(Path(sys.argv[2]),Path(sys.argv[3]),out)))
    else:
        discover(Path(sys.argv[1]),out)
