"""Complete value closures for arbitrary-site three-input two-pass variance.

Every addition tree/permutation and every five-model arithmetic annotation is
included. Corrected two-pass bounds deliberately drop the dependency between
sum(delta) and sum(delta^2); an exclusion of that larger set is still valid.
No arithmetic value is inferred from an expected answer.
"""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from bessel import operation,MODELS,pair,digest,require
from intervals import B64,E64,E53,Interval,number,pow2


def encode(values):return [pair(x)for x in sorted(values)]


def sums(leaves):
    """Subset closure includes both ordered children at every binary node."""
    table={1<<i:set(values)for i,values in enumerate(leaves)}
    for mask in range(1,8):
        if mask in table:continue
        result=set();left=(mask-1)&mask
        while left:
            right=mask^left
            if right:
                for a,b,model in itertools.product(table[left],table[right],MODELS):
                    result.add(operation('add',model,a,b))
            left=(left-1)&mask
        table[mask]=result
    return table


def serialized(table):return {str(mask):encode(table[mask])for mask in range(1,8)}


def close(row):
    def flatten(arg):return list(itertools.chain.from_iterable(flatten(x)for x in arg))if isinstance(arg,list)else[arg]
    args=flatten(row['args'])
    require(len(args)==3 and all(isinstance(x,str)and x.startswith('0x')for x in args), 'outside numeric three-input scope')
    xs=[Q(number(x))for x in args]
    require(all(x==0 or pow2(-100)<=abs(x)<=pow2(100)for x in xs), 'outside input magnitude bound')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'), 'unknown variance function')
    population=row['function']in('VAR.P','VARP','VARPA')
    divisor=Q(3 if population else 2)
    total=sums([{x}for x in xs])
    means={operation('div',model,s,Q(3))for s,model in itertools.product(total[7],MODELS)}
    plain=set();corrected=set();details=[]
    for mean in sorted(means):
        ds=[{operation('add',model,x,-mean)for model in MODELS}for x in xs]
        squares=[{operation('mul',model,d,d)for d,model in itertools.product(values,MODELS)}for values in ds]
        q=sums(squares);t=sums(ds)
        correction={operation('div',md,operation('mul',mm,d,d),Q(3))
                    for d,mm,md in itertools.product(t[7],MODELS,MODELS)}
        outputs={operation('div',model,value,divisor)for value,model in itertools.product(q[7],MODELS)}
        adjusted={operation('div',md,operation('add',ms,value,-c),divisor)
                  for value,c,ms,md in itertools.product(q[7],correction,MODELS,MODELS)}
        details.append(dict(mean=pair(mean),deltas=list(map(encode,ds)),squares=list(map(encode,squares)),
                            square_sum=serialized(q),delta_sum=serialized(t),corrections=encode(correction),
                            plain_roots=encode(outputs),corrected_roots=encode(adjusted)))
        plain.update(outputs);corrected.update(adjusted)
    target=Q(number(row['expected']));span=B64.preimage(Interval(target,target))
    reachable=lambda roots:any(span.contains(value)for value in roots)
    return dict(raw=row,mean_sum=serialized(total),means=encode(means),per_mean=details,
                plain_roots=encode(plain),corrected_roots=encode(corrected),target_span=span.json(),
                plain_reachable=reachable(plain),corrected_reachable=reachable(corrected))


def run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    counts=Counter()
    with(dest/'closures.jsonl').open('x')as out:
        for row in rows:
            proof=close(row);out.write(json.dumps(proof)+'\n');counts['observations']+=1
            counts['plain_unreachable']+=not proof['plain_reachable'];counts['corrected_unreachable']+=not proof['corrected_reachable']
            counts['mean_values']+=len(proof['means']);counts['plain_root_values']+=len(proof['plain_roots'])
            counts['corrected_root_values']+=len(proof['corrected_roots'])
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,counts=counts,seconds=time.monotonic()-start,
                source_sha256=source,input_sha256={str(input_path):digest(input_path)},
                claim='All binary sum trees/permutations and arbitrary five-model assignments at every arithmetic site for three-input plain two-pass variance. '
                      'Corrected two-pass uses a conservative dependency-free closure containing every consistent delta/square program.',
                limits='Three-input discovery observations and listed algebra only. Reachable roots need not share a reusable program; conservative corrected roots may not have consistent leaf choices. '
                       'Zero signs are immaterial to squares and division by positive integers. No branch, fit or complete Excel kernel follows.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def online_close(row):
    def flatten(v):return list(itertools.chain.from_iterable(map(flatten,v)))if isinstance(v,list)else[v]
    bits=flatten(row['args']);require(len(bits)==3 and all(isinstance(v,str)and v.startswith('0x')for v in bits),'outside three-input online scope')
    xs=list(map(lambda v:Q(number(v)),bits));require(all(v==0 or pow2(-100)<=abs(v)<=pow2(100)for v in xs),'online magnitude scope')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'),'wrong online variance function')
    divisor=Q(3 if row['function']in('VAR.P','VARP','VARPA')else 2)
    def values(kind,a,b):return{operation(kind,m,a,b)for m in MODELS}
    details=[];roots={'welford':set(),'youngs_cramer':set()}
    for perm in itertools.permutations(range(3)):
        record=dict(indices=list(perm))
        for family in roots:
            states={(xs[perm[0]],Q(0))};stages=[]
            for k,index in enumerate(perm[1:],2):
                x=xs[index];following=set()
                for accumulator,moment in states:
                    if family=='welford':
                        for delta in values('add',x,-accumulator):
                            means={v for q in values('div',delta,Q(k))for v in values('add',accumulator,q)}
                            for mean in means:
                                increments={v for e in values('add',x,-mean)for v in values('mul',delta,e)}
                                for t in increments:
                                    following.update((mean,q)for q in values('add',moment,t))
                    else:
                        for total in values('add',accumulator,x):
                            deviations={v for scaled in values('mul',Q(k),x)for v in values('add',scaled,-total)}
                            increments={v for delta in deviations for square in values('mul',delta,delta)for v in values('div',square,Q(k*(k-1)))}
                            for t in increments:following.update((total,q)for q in values('add',moment,t))
                states=following;stages.append([[pair(a),pair(b)]for a,b in sorted(states)])
            outputs={v for _,q in states for v in values('div',q,divisor)}
            record[family]=dict(stages=stages,roots=encode(outputs));roots[family].update(outputs)
        details.append(record)
    target=Q(number(row['expected']));span=B64.preimage(Interval(target,target))
    return dict(raw=row,permutations=details,roots={k:encode(v)for k,v in roots.items()},target_span=span.json(),reachable={k:any(span.contains(v)for v in values)for k,values in roots.items()})


def online_run(input_path,destination):
    start=time.monotonic();dest=Path(destination);dest.mkdir();source=digest(Path(__file__))
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()];require(rows,'empty online observation scope')
    counts=Counter()
    with(dest/'closures.jsonl').open('x')as out:
        for row in rows:
            proof=online_close(row);out.write(json.dumps(proof)+'\n');counts['observations']+=1
            for family in('welford','youngs_cramer'):
                counts[family+'_unreachable']+=not proof['reachable'][family];counts[family+'_roots']+=len(proof['roots'][family])
                counts[family+'_state_pairs']+=sum(len(stage)for p in proof['permutations']for stage in p[family]['stages'])
    require(digest(Path(__file__))==source,'source changed during online closure')  # ubs:ignore — public integrity hashes, not secrets
    summary=dict(complete=True,counts=counts,seconds=time.monotonic()-start,source_sha256=source,input_sha256={str(input_path):digest(input_path)},
                 limits='Three-input Welford and Youngs-Cramer unrolled updates, all6 permutations and independent five-model annotations at every site. Boundary states preserve mean/sum-second-moment correlations. Reachability is not a reusable producer, fitted value or product gain.')
    with(dest/'summary.json').open('x')as out:json.dump(summary,out,indent=2);out.write('\n')
    print(json.dumps(summary))


def raw_moment_close(row, fused=False):
    """All-site closures of six fixed algebraic variance identities."""
    def flatten(v):
        return list(itertools.chain.from_iterable(map(flatten, v))) if isinstance(v, list) else [v]
    bits = flatten(row['args'])
    require(len(bits) == 3 and all(isinstance(v, str) and v.startswith('0x') for v in bits), 'raw-moment three-input scope')
    xs = [Q(number(v)) for v in bits]
    require(all(v == 0 or pow2(-100) <= abs(v) <= pow2(100) for v in xs), 'raw-moment magnitude scope')
    require(row['function'] in ('VAR', 'VAR.S', 'VAR.P', 'VARP', 'VARA', 'VARPA'), 'raw-moment variance function')
    d = Q(3 if row['function'] in ('VAR.P', 'VARP', 'VARPA') else 2)
    def combine(kind, left, right):
        return {operation(kind, model, a, b) for a, b, model in itertools.product(left, right, MODELS)}
    def subtract(left, right):
        return combine('add', left, {-v for v in right})
    def square(values):
        return {operation('mul', model, v, v) for v, model in itertools.product(values, MODELS)}
    def fma(a, b, c):
        result = set()
        for model in MODELS:
            if model == 'strict':
                a0, b0, c0 = map(B64.rn, (a, b, c))
                result.add(B64.rn(a0*b0+c0))
            else:
                value = (E53 if model.endswith('pc53') else E64).rn(a*b+c)
                result.add(B64.rn(value) if model.startswith('stored') else value)
        return result
    total = sums([{v} for v in xs])[7]
    squares = sums([square({v}) for v in xs])[7]
    if fused:
        for order in itertools.permutations(xs):
            states = {Q(0)}
            for value in order:
                states = set().union(*(fma(value, value, state) for state in states))
            squares.update(states)
    total_square = square(total)
    nd = combine('mul', {Q(3)}, {d})
    roots = []
    # Independent sums and their correlations are conservatively projected.
    roots.append(combine('div', subtract(squares, combine('div', total_square, {Q(3)})), {d}))
    roots.append(combine('div', subtract(combine('mul', {Q(3)}, squares), total_square), nd))
    roots.append(subtract(combine('div', squares, {d}), combine('div', total_square, nd)))
    divided_squares = combine('div', squares, {d})
    mean_products = set()
    for value in total:
        mean_products.update(combine('mul', combine('div', {value}, {Q(3)}), combine('div', {value}, {d})))
    roots.append(subtract(divided_squares, mean_products))
    means = combine('div', total, {Q(3)})
    raw_population = subtract(combine('div', squares, {Q(3)}), square(means))
    roots.append(combine('mul', raw_population, combine('div', {Q(3)}, {d})))
    differences = [square(subtract({xs[i]}, {xs[j]})) for i, j in ((0, 1), (0, 2), (1, 2))]
    roots.append(combine('div', sums(differences)[7], nd))
    if fused:
        numerators = [set() for _ in range(4)]
        for value in total:
            for mean in combine('div', {value}, {Q(3)}):
                for q in squares:
                    numerators[0].update(fma(-mean, value, q))
                for q in combine('div', squares, {Q(3)}):
                    numerators[3].update(fma(-mean, mean, q))
            for q in combine('mul', {Q(3)}, squares):
                numerators[2].update(fma(-value, value, q))
        for q, s2 in itertools.product(squares, total_square):
            numerators[1].update(fma(Q(3), q, -s2))
        roots.extend((combine('div', numerators[0], {d}),
                      combine('div', numerators[1], nd),
                      combine('div', numerators[2], nd),
                      combine('mul', numerators[3], combine('div', {Q(3)}, {d}))))
    target = Q(number(row['expected']))
    span = B64.preimage(Interval(target, target))
    import hashlib
    forms = []
    for values in roots:
        payload = json.dumps(encode(values), separators=(',', ':')).encode()
        forms.append(dict(root_count=len(values), roots_sha256=hashlib.sha256(payload).hexdigest(),
                          reachable=any(span.contains(v) for v in values)))
    return dict(raw=row, target_span=span.json(), forms=forms,
                total_values=len(total), square_sum_values=len(squares))


def raw_moment_run(input_path, destination, fused=False):
    start = time.monotonic()
    dest = Path(destination)
    dest.mkdir()
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__), Path(__file__).with_name('bessel.py'), Path(__file__).with_name('intervals.py'))}
    counts, by_function = Counter(), {}
    rows = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    require(rows, 'empty raw-moment scope')
    with (dest/'closures.jsonl').open('x') as out:
        for row in rows:
            result = raw_moment_close(row, fused)
            out.write(json.dumps(result)+'\n')
            out.flush()
            function = by_function.setdefault(row['function'], Counter())
            counts['observations'] += 1
            function['observations'] += 1
            for form, values in enumerate(result['forms']):
                counts['form_'+str(form)+'_unreachable'] += not values['reachable']
                function['form_'+str(form)+'_unreachable'] += not values['reachable']
                counts['root_values'] += values['root_count']
    require(frozen == {p: digest(p) for p in frozen}, 'raw-moment source/input changed')
    summary = dict(complete=True, fused=fused, counts=counts, by_function=by_function,
        seconds=time.monotonic()-start, input_sha256=frozen,
        forms=['(Q-S*S/n)/d', '(n*Q-S*S)/(n*d)', 'Q/d-S*S/(n*d)', 'Q/d-(S/n)*(S/d)', '(Q/n-(S/n)^2)*(n/d)', 'sum_pairs((xi-xj)^2)/(n*d)']
            + (['FMA(-S/n,S,Q)/d', 'FMA(n,Q,-S*S)/(n*d)', 'FMA(-S,S,n*Q)/(n*d)', 'FMA(-S/n,S/n,Q/n)*(n/d)'] if fused else []),
        limits='Three numeric inputs; '+('ten forms including FMA' if fused else 'six identities')
            +'; all sum trees/permutations and independent five-model arithmetic choices at every site. Shared-sum correlations are discarded where documented, enlarging the sets. Root digests are checked by complete independent regeneration. Reachability supplies no reusable program, fitted value, product gain or tolerance claim.')
    (dest/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


def raw_moment_scope(input_path, destination):
    """Select canonical bank rows for the existing three-input closure."""
    dest = Path(destination)
    dest.mkdir()
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__),
        Path(__file__).with_name('intervals.py'))}
    counts, selected = Counter(), []
    def flat(value):
        return [leaf for child in value for leaf in flat(child)] if isinstance(value, list) else [value]
    with Path(input_path).open() as source, (dest/'cases.jsonl').open('x') as out, (dest/'excluded.jsonl').open('x') as excluded:
        for index, line in enumerate(source):
            raw = json.loads(line)
            reason = None
            values = flat(raw['args'])
            if raw['function'] not in ('VAR', 'VAR.S', 'VAR.P', 'VARP', 'VARA', 'VARPA'):
                reason = 'other_function'
            elif len(values) != 3:
                reason = 'other_input_count'
            elif not all(isinstance(v, str) and v.startswith('0x') and len(v) == 18 for v in values):
                reason = 'non_numeric_input'
            elif not isinstance(raw['expected'], str) or not raw['expected'].startswith('0x'):
                reason = 'non_numeric_expected'
            else:
                try:
                    xs = [Q(number(v)) for v in values]
                    Q(number(raw['expected']))
                    if not all(x == 0 or pow2(-100) <= abs(x) <= pow2(100) for x in xs):
                        reason = 'outside_magnitude_scope'
                except (ValueError, OverflowError):
                    reason = 'non_finite_value'
            if reason:
                counts[reason] += 1
                excluded.write(json.dumps(dict(trace_index=index, reason=reason))+'\n')
            else:
                counts['selected:'+raw['function']] += 1
                selected.append(index)
                out.write(json.dumps(raw)+'\n')
    require(selected and frozen == {p: digest(p) for p in frozen}, 'empty or changed raw-moment selection')
    summary = dict(complete=True, counts=counts, selected_indices=selected, input_sha256=frozen,
        output_sha256={p.name:digest(p) for p in dest.iterdir()},
        limits='Only original decoded observations are selected; all omitted indices and reasons are retained. '
               'Original arguments, shapes, metadata and answers are unchanged. This scope is not a bank score or new capture.')
    (dest/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(dict(complete=True, counts=counts)))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='raw-moment-scope':raw_moment_scope(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='raw-moments-fma':raw_moment_run(*sys.argv[2:], fused=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='raw-moments':raw_moment_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='online':online_run(*sys.argv[2:])
    else:run(*sys.argv[1:])
