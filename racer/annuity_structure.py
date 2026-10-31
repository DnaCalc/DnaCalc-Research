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
from fractions import Fraction as Q

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


def relationships(trace_dir,out):
    """Loan identities and externally probed RATE/FV companions, not a solver trace."""
    from intervals import B64
    from intersect import rnd
    out.mkdir(parents=True,exist_ok=True)
    data={fn:list(lines(trace_dir/(fn+'.jsonl')))for fn in FUNCTIONS}
    env=lambda r:json.dumps(r['reading'],sort_keys=True)
    full=lambda a:json.dumps(a,sort_keys=True,separators=(',',':'))
    pmt={(full(r['args']),env(r)):r for r in data['PMT']}
    ipmt={(full(r['args']),env(r)):r for r in data['IPMT']}
    links=[];sums=[]
    for fn in ['IPMT','PPMT']:
        for row in data[fn]:
            a=row['args'];key=full([a[0],a[2],*a[3:]])if len(a)==6 else None
            mate=pmt.get((key,env(row)))if key else None
            links.append(dict(function=fn,id=row['id'],PMT_companion=mate['id']if mate else None,
                own_error=row['expected']!=row['actual'],companion_error=mate['expected']!=mate['actual']if mate else None))
            if fn=='PPMT':
                interest=ipmt.get((full(a),env(row)))
                if interest and mate and all(num(x['expected'])is not None for x in [row,interest,mate]):
                    q=Q(num(row['expected']))+Q(num(interest['expected']));rounded=B64.rn(q)
                    sums.append(dict(PPMT=row['id'],IPMT=interest['id'],PMT=mate['id'],
                        rounded_recorded_sum=bits(float(rounded))if rounded is not None else'error:Num',
                        PMT_expected=mate['expected'],same=rounded==Q(num(mate['expected']))))
    first=[];first_models=C.Counter()
    for row in data['IPMT']:
        a=[num(v)for v in row['args']]
        if len(a)!=6 or any(v is None for v in a):continue
        rate,period,n,pv,fv,timing=a
        if period!=1 or rate<=-1 or n<1:continue
        target=row['expected'];models=[]
        for m in range(5):
            q=Q(0)if rate==0 or math.trunc(timing)!=0 else rnd(-Q(pv)*Q(rate),m)
            published=B64.rn(q)if q is not None else None
            got=bits(float(published))if published is not None else'error:Num'
            if got=='0x8000000000000000':got='0x0000000000000000'
            first_models[(m,got==target)]+=1;models.append(dict(model=m,prediction=got,match=got==target))
        first.append(dict(witness=row,product_models=models))
    fv_cash=C.defaultdict(list)
    for row in data['FV']:
        if len(row['args'])==5:
            fv_cash[(full(row['args'][1:]),env(row))].append(row)
    rate_rows=[];companion_totals=C.Counter();residuals=[]
    for row in data['RATE']:
        a=row['args'];future=a[3]if len(a)>3 else None
        if len(a)<5:
            companion_totals['RATE_shape_outside_direct_companion_key']+=1;continue
        companions=fv_cash.get((full([a[0],a[1],a[2],a[4]]),env(row)),[])
        companion_totals['RATE_rows']+=1;companion_totals['RATE_rows_with_FV_companions']+=bool(companions)
        companion_totals['FV_companion_links']+=len(companions)
        counted=C.Counter();records=[]
        for fvrow in companions:
            rate=fvrow['args'][0]
            at='recorded_RATE_result'if rate==row['expected']else'current_RATE_result'if rate==row['actual']else'explicit_input_guess'if len(a)>5 and rate==a[5]else'other_external_probe'
            counted[at]+=1;companion_totals[at]+=1
            rec=dict(FV_id=fvrow['id'],rate=rate,at=at,FV_error=fvrow['expected']!=fvrow['actual'],source=fvrow['source'])
            if at in ['recorded_RATE_result','current_RATE_result']and num(future)is not None and num(fvrow['expected'])is not None:
                residual=Q(num(fvrow['expected']))-Q(num(future))
                rec['exact_recorded_FV_residual']=[str(residual.numerator),str(residual.denominator)]
                rec['residual_float_diagnostic']=float(residual)
                residuals.append(dict(RATE_id=row['id'],**rec))
            records.append(rec)
        rate_rows.append(dict(RATE_id=row['id'],RATE_error=row['expected']!=row['actual'],FV_companion_counts=dict(counted),companions=records))
    result=dict(loan_links=links,exact_recorded_payment_partitions=sums,
        loan_link_totals={fn:dict(rows=sum(x['function']==fn for x in links),with_PMT_companion=sum(x['function']==fn and x['PMT_companion']is not None for x in links))for fn in ['IPMT','PPMT']},
        first_period_IPMT_observations=first,
        first_period_product_model_counts=[dict(model=m,matched=hit,rows=n)for(m,hit),n in sorted(first_models.items())],
        RATE_FV_companion_totals=dict(companion_totals),RATE_FV_companions=rate_rows,RATE_residual_observations=residuals,
        limitations=['Full recorded environment and every cashflow/type bit must match for a companion link.',
            'Externally chosen FV guesses/perturbations are observations of FV, not observed internal Excel RATE evaluations. No claim about which points Excel visited.',
            'Residuals subtract the future-value target from a recorded binary64 FV value exactly; they do not expose the internal RATE residual or its stopping tolerance.',
            'Missing companions are coverage limits, not matched calculations; first-period identities alone do not identify later IPMT/PPMT helper programmes.'])
    (out/'family-relationships.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k]for k in ['loan_link_totals','first_period_product_model_counts','RATE_FV_companion_totals']}))


def homogeneity(trace_dir,out):
    """Exact power-of-two/sign scaling diagnostics on original recorded bits."""
    result={}
    for fn in FUNCTIONS:
        table={};counts=C.Counter();failures=[];wanted=set()
        invariant=fn in ['NPER','RATE']
        cash_indices=[1,2,3]if invariant else[3,4]if fn in ['IPMT','PPMT']else[2,3]
        for idx,row in enumerate(lines(trace_dir/(fn+'.jsonl'))):
            counts['recorded_rows']+=1;a=[num(v)for v in row['args']];y=num(row['expected'])
            if any(v is None or not math.isfinite(v)for v in a)or y is None or not math.isfinite(y)or y==0:
                counts['typed_or_nonfinite_or_zero_output_control']+=1;continue
            if len(a)!=(6 if fn in ['IPMT','PPMT','RATE']else 5):counts['shape_control']+=1;continue
            if fn=='RATE':rate=a[5];n=a[0]
            elif fn=='NPER':rate=a[0];n=y
            else:rate=a[0];n=a[2]if fn in ['IPMT','PPMT']else a[1]
            if not(2**-40<=rate<=0.25 and 1<=n<=2**20 and 0<n*math.log1p(rate)<1):
                counts['outside_bounded_rate_period_region']+=1;continue
            if fn in ['IPMT','PPMT']and not(1<=a[1]<=n):counts['period_domain_control']+=1;continue
            amounts=[a[j]for j in cash_indices]
            if any(x and not(2**-128<=abs(x)<=2**128)for x in amounts)or not any(amounts):
                counts['outside_bounded_amount_region']+=1;continue
            rawbits=[int(row['args'][j],16)for j in cash_indices]
            anchor=next(b for b in rawbits if b&0x7fffffffffffffff)
            ae=((anchor>>52)&2047)-1023;sign=anchor>>63
            def normalize(b):
                if b&0x7fffffffffffffff==0:return ('zero',b>>63)
                return (b>>63^sign,b&0x000fffffffffffff,((b>>52)&2047)-1023-ae)
            def output(text):
                value=num(text)
                if value is None or not math.isfinite(value)or value==0:return None
                b=int(text,16)
                return b if invariant else normalize(b)
            key=(tuple(row['args'][j]for j in range(len(a))if j not in cash_indices),tuple(normalize(b)for b in rawbits),json.dumps(row['reading'],sort_keys=True),row['source'])
            ey=output(row['expected']);ay=output(row['actual']);counts['eligible_rows']+=1
            if key not in table:
                table[key]=dict(first=idx,scale=(ae,sign),expected=ey,actual=ay,count=1,contradiction=False,actual_contradiction=False)
                continue
            rec=table[key];rec['count']+=1
            if rec['count']==2:counts['compared_groups']+=1
            counts['additional_scaled_points']+=1
            counts['distinct_scale_or_sign_pairs']+=(ae,sign)!=rec['scale']
            if ey!=rec['expected']:
                counts['contradicting_points']+=1
                if not rec['contradiction']:counts['contradicting_groups']+=1;rec['contradiction']=True
                if len(failures)<20:
                    wanted.add(rec['first']);failures.append(dict(first_trace_index=rec['first'],new_trace_index=idx,new_witness=row,first_normalized_output=rec['expected'],new_normalized_output=ey,source=row['source']))
            if ay is not None and rec['actual']is not None and ay!=rec['actual']:
                counts['current_kernel_contradicting_points']+=1
                if not rec['actual_contradiction']:counts['current_kernel_contradicting_groups']+=1;rec['actual_contradiction']=True
        witnesses={idx:row for idx,row in enumerate(lines(trace_dir/(fn+'.jsonl')))if idx in wanted}
        for f in failures:f['first_witness']=witnesses[f['first_trace_index']]
        counts['unique_normalized_groups']=len(table)
        result[fn]=dict(counts=dict(counts),examples=failures,
            property='Output is invariant under common cashflow power-of-two/sign scaling'if invariant else'Output scales by the same exact power-of-two/sign factor as every cashflow',
            limitations=['Original typed arguments/full environment/source and all noncashflow arguments are fixed. Comparisons use exact significand/exponent/sign identities, not approximate ratios.',
                'Only bounded rate/period/amount, finite nonzero-output cases are compared. No absence claim outside these controls.',
                'A recorded violation suggests amount-dependent arithmetic/entry or solver stopping; it does not by itself prove a broken capture or identify an operation.',
                'RATE can break mathematical scaling through absolute stopping criteria; these records are not internal solver traces.'])
    (out/'scale-homogeneity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({fn:x['counts']for fn,x in result.items()}))


if __name__ == '__main__':
    out=Path(sys.argv[-1]);out.mkdir(exist_ok=True,parents=True)
    if sys.argv[1]=='gaps':
        print(json.dumps(gap_clusters(Path(sys.argv[2]),Path(sys.argv[3]),out)))
    elif sys.argv[1]=='relationships':relationships(Path(sys.argv[2]),out)
    elif sys.argv[1]=='homogeneity':homogeneity(Path(sys.argv[2]),out)
    else:
        discover(Path(sys.argv[1]),out)
