"""Independent full-set proof using explicit permutations/trees, not subset DP."""
from collections import Counter
from fractions import Fraction as R
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys
import time

from exact import BINARY64 as B,EXTENDED64 as E,EXTENDED53 as P53,Span,binary64,operation,two
from check import records,require,digest

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


@lru_cache(maxsize=262144)
def op(kind,model,a,b):return operation(a,b,kind,model)


def encode(values):return [[str(x.numerator),str(x.denominator)]for x in sorted(values)]


def all_sums(leaves):
    table={}
    @lru_cache(maxsize=None)
    def ordered(indices):
        if len(indices) == 1:
            return frozenset(leaves[indices[0]])
        values = set()
        for split in range(1, len(indices)):
            for a, b, model in itertools.product(ordered(indices[:split]), ordered(indices[split:]), MODELS):
                values.add(op('add', model, a, b))
        return frozenset(values)
    for mask in range(1,1 << len(leaves)):
        indices=[i for i in range(len(leaves))if mask&(1<<i)];values=set()
        for perm in itertools.permutations(indices):
            selected=[leaves[i]for i in perm]
            if len(perm)==1:values.update(selected[0])
            elif len(perm)==2:
                for a,b,m in itertools.product(*selected,MODELS):values.add(op('add',m,a,b))
            else:
                if len(perm) == 3:
                    for a,b,c,first,last in itertools.product(*selected,MODELS,MODELS):
                        values.add(op('add',last,op('add',first,a,b),c))
                        values.add(op('add',last,a,op('add',first,b,c)))
                else:
                    values.update(ordered(perm))
        table[str(mask)]=values
    return table


def encoded(table):return {mask:encode(values)for mask,values in table.items()}


def verify(row,proof):
    require(proof['raw']==row,'changed input, source or answer')
    def flatten(v):return list(itertools.chain.from_iterable(map(flatten,v)))if isinstance(v,list)else[v]
    bits=flatten(row['args']);require(len(bits)==3,'outside three-input scope')
    xs=list(map(binary64,bits));require(all(x==0 or two(-100)<=abs(x)<=two(100)for x in xs),'magnitude scope')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'),'wrong function')
    divisor=R(3)if row['function']in('VAR.P','VARP','VARPA')else R(2)
    total=all_sums([{x}for x in xs]);require(encoded(total)==proof['mean_sum'],'incomplete mean-sum closure')
    means={op('div',m,s,R(3))for s,m in itertools.product(total['7'],MODELS)}
    require(encode(means)==proof['means']and len(means)==len(proof['per_mean']),'incomplete shared-mean closure')
    plain=set();corrected=set();nodes=0
    for mean,detail in zip(sorted(means),proof['per_mean']):
        require(encode({mean})[0]==detail['mean'],'wrong mean')
        delta=[{op('add',m,x,-mean)for m in MODELS}for x in xs]
        square=[{op('mul',m,d,d)for d,m in itertools.product(values,MODELS)}for values in delta]
        q,t=all_sums(square),all_sums(delta)
        corrections={op('div',b,op('mul',a,d,d),R(3))for d,a,b in itertools.product(t['7'],MODELS,MODELS)}
        unadjusted={op('div',m,v,divisor)for v,m in itertools.product(q['7'],MODELS)}
        adjusted={op('div',b,op('add',a,v,-c),divisor)for v,c,a,b in itertools.product(q['7'],corrections,MODELS,MODELS)}
        require(list(map(encode,delta))==detail['deltas']and list(map(encode,square))==detail['squares'],'incomplete correlated leaf closure')
        require(encoded(q)==detail['square_sum']and encoded(t)==detail['delta_sum'],'incomplete sum-tree closure')
        require(encode(corrections)==detail['corrections']and encode(unadjusted)==detail['plain_roots']and
                encode(adjusted)==detail['corrected_roots'],'incomplete reconstruction closure')
        plain.update(unadjusted);corrected.update(adjusted)
        nodes+=sum(len(v)for v in q.values())+sum(len(v)for v in t.values())+len(corrections)+len(unadjusted)+len(adjusted)
    require(encode(plain)==proof['plain_roots']and encode(corrected)==proof['corrected_roots'],'incomplete root union')
    target=binary64(row['expected']);span=B.pullback(Span(target,target))
    require(span.json()==proof['target_span'],'incorrect exact target rounding interval')
    # Comparing exact rational roots to the open/closed rounding interval also
    # verifies every final publication, including tie endpoint parity.
    reachable=lambda values:any(span.contains(v)for v in values)
    require(proof['plain_reachable']==reachable(plain)and proof['corrected_reachable']==reachable(corrected),'invalid reachability claim')
    return dict(plain_unreachable=int(not reachable(plain)),corrected_unreachable=int(not reachable(corrected)),
                mean_values=len(means),plain_root_values=len(plain),corrected_root_values=len(corrected)),nodes


def check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    sources={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));proofs=list(records(directory/'closures.jsonl'))
    require(len(rows)==len(proofs),'missing original obligation')
    counts=Counter();states=0
    for row,proof in zip(rows,proofs):
        delta,nodes=verify(row,proof);counts.update(delta);counts['observations']+=1;states+=nodes
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']is True and summary['counts']==counts,'incorrect closure summary')
    require(all(digest(p)==h for p,h in sources.items()),'checker changed during proof')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,counts=counts,full_width_value_states=states,seconds=time.monotonic()-start,
                 arithmetic_cache=op.cache_info()._asdict(),source_sha256=sources,
                 input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 limits='Every specified three-input arithmetic tree/site assignment is covered. Corrected closures are dependency-free supersets. '
                        'No target is used to produce an intermediate, and reachability is not a reusable program or a product gain.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def online_verify(row,proof):
    require(proof['raw']==row,'changed online source observation')
    def flatten(v):return [y for x in v for y in flatten(x)]if isinstance(v,list)else[v]
    bits=flatten(row['args']);require(len(bits)==3 and all(isinstance(v,str)and v.startswith('0x')for v in bits),'online numeric scope')
    xs=list(map(binary64,bits));require(all(v==0 or two(-100)<=abs(v)<=two(100)for v in xs),'online magnitude bound')
    require(row['function']in('VAR','VAR.S','VAR.P','VARP','VARA','VARPA'),'online function scope')
    divisor=R(3)if row['function']in('VAR.P','VARP','VARPA')else R(2)
    def rounded(kind,a,b):return{op(kind,m,a,b)for m in MODELS}
    def advance(family,state,x,k):
        a,q=state;outputs=set()
        if family=='welford':
            differences=rounded('add',x,-a)
            for delta in differences:
                for quotient in rounded('div',delta,R(k)):
                    for mean in rounded('add',a,quotient):
                        for residual in rounded('add',x,-mean):
                            for product in rounded('mul',delta,residual):
                                outputs.update((mean,value)for value in rounded('add',q,product))
        else:
            for total in rounded('add',a,x):
                for multiple in rounded('mul',R(k),x):
                    for residual in rounded('add',multiple,-total):
                        for square in rounded('mul',residual,residual):
                            for contribution in rounded('div',square,R(k*(k-1))):
                                outputs.update((total,value)for value in rounded('add',q,contribution))
        return outputs
    permutations=list(itertools.permutations(range(3)))
    require([p['indices']for p in proof['permutations']]==[list(p)for p in permutations],'incomplete online permutations')
    unions={k:set()for k in('welford','youngs_cramer')};counts=Counter()
    for indices,record in zip(permutations,proof['permutations']):
        for family in unions:
            states={(xs[indices[0]],R(0))};stages=[]
            for k in(2,3):
                following=set()
                for state in states:following.update(advance(family,state,xs[indices[k-1]],k))
                states=following;stages.append([[encode({a})[0],encode({b})[0]]for a,b in sorted(states)])
                counts[family+'_state_pairs']+=len(states)
            require(record[family]['stages']==stages,'incomplete correlated online boundary states')
            values={v for _,q in states for v in rounded('div',q,divisor)}
            require(record[family]['roots']==encode(values),'incomplete online per-order roots')
            unions[family].update(values)
    span=B.pullback(Span(binary64(row['expected']),binary64(row['expected'])))
    require(span.json()==proof['target_span'],'changed online target preimage')
    require(set(proof['roots'])==set(unions)and set(proof['reachable'])==set(unions),'changed online family coverage')
    for family,values in unions.items():
        require(proof['roots'][family]==encode(values),'incomplete online root union')
        reachable=any(span.contains(v)for v in values);require(proof['reachable'][family]==reachable,'incorrect online exclusion')
        counts[family+'_roots']+=len(values);counts[family+'_unreachable']+=not reachable
    counts['observations']+=1
    return counts


def online_check(input_path,directory,output):
    start=time.monotonic();directory=Path(directory);source=Path(__file__)
    hashes={str(p):digest(p)for p in(source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));proofs=list(records(directory/'closures.jsonl'));require(rows and len(rows)==len(proofs),'missing online obligations')
    counts=Counter()
    for row,proof in zip(rows,proofs):counts.update(online_verify(row,proof))
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True and summary['counts']==counts,'incorrect online summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during online closure')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,counts=counts,seconds=time.monotonic()-start,source_sha256=hashes,
                 input_sha256={str(p):digest(p)for p in(Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 limits='Independent exact boundary-state reproduction for the specified two unrolled online algorithms. All per-site choices and6 index permutations covered; no target-derived state, complete kernel identification or product gain.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def raw_moment_result(row, fused=False, root=False):
    """Independent explicit-tree enumeration of the six uncentered identities."""
    require(not root or fused and row.get('result_index') is None, 'STDEV scope requires fused forms and a scalar call')
    import hashlib
    def flatten(v):
        return list(itertools.chain.from_iterable(map(flatten, v))) if isinstance(v, list) else [v]
    bits = flatten(row['args'])
    require(len(bits) == (4 if root else 3) and all(isinstance(v, str) and v.startswith('0x') for v in bits), 'raw arity scope')
    xs = list(map(binary64, bits))
    require(all(x == 0 or two(-100) <= abs(x) <= two(100) for x in xs), 'raw magnitude scope')
    names = ('STDEV','STDEV.S','STDEV.P','STDEVP','STDEVA','STDEVPA') if root else ('VAR','VAR.S','VAR.P','VARP','VARA','VARPA')
    require(row['function'] in names, 'raw function scope')
    size = len(xs)
    denominator = R(size if row['function'] in ('VAR.P','VARP','VARPA','STDEV.P','STDEVP','STDEVPA') else size-1)
    full = str((1 << size)-1)
    def cross(kind, left, right):
        result = set()
        for model in MODELS:
            for a in left:
                for b in right:
                    result.add(op(kind, model, a, b))
        return result
    def minus(left, right):
        return cross('add', left, {-b for b in right})
    def squared(values):
        return {op('mul', model, a, a) for model in MODELS for a in values}
    def fused_results(a, b, c):
        result = set()
        for model, fmt in [('strict', B), ('stored_pc64', E), ('resident_pc64', E),
                           ('stored_pc53', P53), ('resident_pc53', P53)]:
            if model == 'strict':
                exact = B.round(a)*B.round(b)+B.round(c)
            else:
                exact = a*b+c
            value = fmt.round(exact)
            result.add(B.round(value) if model.startswith('stored') else value)
        return result
    s = all_sums([{x} for x in xs])[full]
    q = all_sums([squared({x}) for x in xs])[full]
    if fused:
        def square_tree(tree):
            if isinstance(tree, int):
                return squared({xs[tree]})
            left, right = tree
            a, b = square_tree(left), square_tree(right)
            values = cross('add', a, b)
            for leaf, sibling in ((left, b), (right, a)):
                if isinstance(leaf, int):
                    for accumulator in sibling:
                        values |= fused_results(xs[leaf], xs[leaf], accumulator)
            return values
        q = set()
        def trees(ordering):
            if len(ordering) == 1:
                yield ordering[0]
            else:
                for split in range(1, len(ordering)):
                    for a, b in itertools.product(trees(ordering[:split]), trees(ordering[split:])):
                        yield (a, b)
        for ordering in itertools.permutations(range(size)):
            for tree in trees(ordering):
                q |= square_tree(tree)
    s2 = squared(s)
    n, d = {R(size)}, {denominator}
    nd = cross('mul', n, d)
    results = [cross('div', minus(q, cross('div', s2, n)), d),
               cross('div', minus(cross('mul', n, q), s2), nd),
               minus(cross('div', q, d), cross('div', s2, nd))]
    products = set()
    for total in s:
        products |= cross('mul', cross('div', {total}, n), cross('div', {total}, d))
    results.append(minus(cross('div', q, d), products))
    centered = minus(cross('div', q, n), squared(cross('div', s, n)))
    results.append(cross('mul', centered, cross('div', n, d)))
    if not root:
        pair_terms = [squared(minus({a}, {b})) for a, b in itertools.combinations(xs, 2)]
        results.append(cross('div', all_sums(pair_terms)['7'], nd))
    if fused:
        extra = [set() for _ in range(4)]
        for total in s:
            mean_values = cross('div', {total}, n)
            for mean, squares in itertools.product(mean_values, q):
                extra[0] |= fused_results(-mean, total, squares)
            for mean, average_square in itertools.product(mean_values, cross('div', q, n)):
                extra[3] |= fused_results(-mean, mean, average_square)
            for scaled in cross('mul', n, q):
                extra[2] |= fused_results(-total, total, scaled)
        for squares, square_sum in itertools.product(q, s2):
            extra[1] |= fused_results(R(size), squares, -square_sum)
        results += [cross('div', extra[0], d), cross('div', extra[1], nd),
                    cross('div', extra[2], nd), cross('mul', extra[3], cross('div', n, d))]
    target = binary64(row['expected'])
    span = B.pullback(Span(target, target))
    forms = []
    for values in results:
        hasher = hashlib.sha256(json.dumps(encode(values), separators=(',', ':')).encode())
        form = dict(root_count=len(values), roots_sha256=hasher.hexdigest(),
                    reachable=any(span.contains(v) for v in values))
        if root:
            if __package__:
                from .geomean import root_publication_domains
            else:
                from geomean import root_publication_domains
            domains = root_publication_domains(values, 2, target, True)
            form.update(root_roundings=domains, reachable=any(d['reachable'] for d in domains.values()))
        forms.append(form)
    return dict(raw=row, target_span=span.json(), forms=forms,
                total_values=len(s), square_sum_values=len(q))


def raw_moment_check(input_path, directory, output, fused=False, root=False):
    started = time.monotonic()
    directory = Path(directory)
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__),
        Path(__file__).with_name('exact.py'), Path(__file__).with_name('check.py'), *sorted(directory.iterdir())) if p.is_file()}
    if root:
        frozen[str(Path(__file__).with_name('geomean.py'))] = digest(Path(__file__).with_name('geomean.py'))
    rows = list(records(input_path))
    proofs = list(records(directory/'closures.jsonl'))
    require(rows and len(rows) == len(proofs), 'missing raw-moment obligations')
    counts, by_function = Counter(), {}
    for row, proof in zip(rows, proofs):
        independently = raw_moment_result(row, fused, root)
        require(independently == proof, 'raw-moment source, root set, digest or reachability differs')
        function = by_function.setdefault(row['function'], Counter())
        counts['observations'] += 1
        function['observations'] += 1
        for form, result in enumerate(independently['forms']):
            counts['form_'+str(form)+'_unreachable'] += not result['reachable']
            function['form_'+str(form)+'_unreachable'] += not result['reachable']
            counts['root_values'] += result['root_count']
    summary = json.loads((directory/'summary.json').read_text())
    require(summary['complete'] is True and summary['fused'] == fused and summary['counts'] == counts and summary['by_function'] == by_function,
            'incorrect raw-moment summary')
    require(summary.get('root_output', False) == root, 'changed variance/sqrt grammar')
    require(frozen == {p: digest(p) for p in frozen}, 'raw-moment evidence changed during check')
    receipt = dict(complete=True, fused=fused, counts=counts, by_function=by_function,
        seconds=time.monotonic()-started, input_sha256=frozen,
        limits='Complete independent explicit-tree/permutation regeneration of every root-set digest, using separate exact lattice arithmetic. '
               + ('Ten declared forms, including fused multiply-add.' if fused else 'Six declared uncentered identities.')
               + (' Q covers ordinary square-sum trees with optional fusion of any leaf square and its sibling sum, including mixed trees.' if fused else '')
               + ' Conservative correlation projections; reachability never supplies a reusable kernel or product match.')
    if root:
        receipt.update(root_output=True, arity=4,
            limits='Independent explicit permutation/tree regeneration of nine four-input uncentered/mixed-FMA variance forms, then correctly rounded real sqrt directly to binary64 or via extended64/extended53 before publication. '
                   'Exact inverse root/publication endpoints are checked; Q includes mixed leaf-square/sibling-sum fusion. '
                   'Pairwise six-term sums, other fusion placements/arities and native approximations are not covered. Counts sum local pre-sqrt variance states, not programs or product matches.')
    with open(output, 'x') as out:
        json.dump(receipt, out, indent=2)
        out.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='raw-moments-stdev4':raw_moment_check(*sys.argv[2:], fused=True, root=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='raw-moments-fma':raw_moment_check(*sys.argv[2:], fused=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='raw-moments':raw_moment_check(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='online':online_check(*sys.argv[2:])
    else:check(*sys.argv[1:])
