"""Exact necessary recurrence constraints with unpublished common Ext80 seeds.

Original mathematics, MIT. DLMF 10.6.1 and 10.29.1. No seed, coefficient,
threshold or answer is fitted or selected. Bounds quantify over all Ext80 seeds
whose binary64 publications agree with independently recorded orders 0 and 1.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import time

try:
    from .intervals import B64, E64, E53, Interval, number, pow2
except ImportError:
    from intervals import B64, E64, E53, Interval, number, pow2

MODELS = ('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pair(value):
    return [str(value.numerator), str(value.denominator)]


def context(row):
    return json.dumps([row['function'],row['args'][0],row['reading'],row['source']],
                      sort_keys=True,separators=(',',':'))


def groups(rows):
    grouped = defaultdict(dict)
    for i,row in enumerate(rows):
        args = row['args']
        if len(args) != 2 or row['result_index'] is not None:
            continue
        if not all(isinstance(a,str) and a.startswith('0x') for a in args):
            continue
        x, n = map(number,args)
        if x <= 0 or n not in (0,1,2) or not row['expected'].startswith('0x'):
            continue
        y = Q(number(row['expected']))
        if abs(y) < pow2(-1022):
            continue
        key = context(row)
        if int(n) in grouped[key]:
            require(rows[grouped[key][int(n)]]['expected'] == row['expected'], 'conflicting same-context observation')
        else:
            grouped[key][int(n)] = i
    return [dict(context=k, indices=[v[n]for n in range(3)])
            for k,v in sorted(grouped.items()) if set(v)=={0,1,2}]


def seed_bounds(target):
    span = B64.preimage(Interval(target,target))
    lo,hi = E64.ceil(span.lo),E64.floor(span.hi)
    if lo==span.lo and not span.lc:lo=E64.next(lo)
    if hi==span.hi and not span.hc:hi=E64.prev(hi)
    require(lo is not None and hi is not None and lo<=hi, 'empty finite seed bounds')
    return lo,hi


@lru_cache(maxsize=131072)
def operation(kind, model, a, b):
    if model=='strict':
        a,b=B64.rn(a),B64.rn(b)
    exact={'mul':lambda:a*b,'div':lambda:a/b,'add':lambda:a+b}[kind]()
    fmt=B64 if model=='strict'else E53 if model.endswith('pc53')else E64
    value=fmt.rn(exact)
    require(value is not None, 'nonfinite recurrence outside scope')
    if model.startswith('stored'):value=B64.rn(value)
    require(value is not None, 'nonfinite stored recurrence outside scope')
    return value


def recurrence(x, zero, one, function, shape, models):
    m0,m1,m2=models
    if shape==0:
        term=operation('div',m1,operation('mul',m0,Q(2),one),x)
    elif shape==1:
        term=operation('mul',m1,Q(2),operation('div',m0,one,x))
    elif shape==2:
        term=operation('mul',m1,operation('div',m0,Q(2),x),one)
    else:
        term=operation('div',m1,one,operation('div',m0,x,Q(2)))
    if function in ('BESSELJ','BESSELY'):
        return operation('add',m2,term,-zero)
    if function=='BESSELI':
        return operation('add',m2,zero,-term)
    require(function=='BESSELK', 'unknown Bessel consumer')
    return operation('add',m2,zero,term)


def profiles():
    for identity,(shape,models) in enumerate(itertools.product(range(4),itertools.product(MODELS,repeat=3))):
        yield dict(id=identity,shape=shape,models=list(models))


def bounds(rows,group,profile):
    a,b,target=[rows[i]for i in group['indices']]
    function=a['function'];x=Q(number(a['args'][0]))
    zero=seed_bounds(Q(number(a['expected'])))
    one=seed_bounds(Q(number(b['expected'])))
    if function in ('BESSELJ','BESSELY'):
        corners=((zero[1],one[0]),(zero[0],one[1]))
    elif function=='BESSELI':
        corners=((zero[0],one[1]),(zero[1],one[0]))
    else:
        corners=((zero[0],one[0]),(zero[1],one[1]))
    roots=[recurrence(x,z,o,function,profile['shape'],profile['models'])for z,o in corners]
    published=list(map(B64.rn,roots))
    require(roots[0]<=roots[1]and published[0]<=published[1], 'nonmonotone corner orientation')
    expected=Q(number(target['expected']))
    return dict(seed0=list(map(pair,zero)),seed1=list(map(pair,one)),
                corner_seeds=[[pair(z),pair(o)]for z,o in corners],
                root_bounds=list(map(pair,roots)),published_bounds=list(map(pair,published)),
                target=pair(expected),excluded=not(published[0]<=expected<=published[1]))


def run(input_path,destination):
    started=time.monotonic();dest=Path(destination);dest.mkdir()
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    gs=groups(rows);ps=list(profiles());counts={};source=digest(Path(__file__))
    for name,obj in [('groups.json',gs),('profiles.json',ps)]:
        with(dest/name).open('x')as out:json.dump(obj,out,indent=2);out.write('\n')
    for function in ('BESSELJ','BESSELY','BESSELI','BESSELK'):
        scoped=[(i,g)for i,g in enumerate(gs)if rows[g['indices'][0]]['function']==function]
        require(scoped, 'empty recurrence cohort')
        survivors=[];executions=0
        with(dest/f'{function}-first-kills.jsonl').open('x')as kills:
            for p in ps:
                for order,(index,g)in enumerate(scoped):
                    result=bounds(rows,g,p);executions+=1
                    if result['excluded']:
                        kills.write(json.dumps(dict(candidate_id=p['id'],order=order,group_index=index,proof=result))+'\n')
                        break
                else:survivors.append(p['id'])
        with(dest/f'{function}-necessary-survivors.json').open('x')as out:json.dump(survivors,out);out.write('\n')
        counts[function]=dict(profiles=len(ps),groups=len(scoped),executions=executions,
                              rejected=len(ps)-len(survivors),necessary_survivors=len(survivors))
    require(digest(Path(__file__))==source, 'source changed during run')  # ubs:ignore — public integrity hashes, not secrets
    result=dict(complete=True,counts=counts,seconds=time.monotonic()-started,
                input_sha256={str(input_path):digest(input_path)},source_sha256=source,
                claim='Necessary constraints on order-two recurrence consumers with arbitrary common Ext80 order-zero/one seeds. '
                      'Monotone corner bounds cover all seed pairs; a surviving bound need not contain any actual seed pair or producer.',
                limits='Same original source/environment and positive input; nonzero normal publications at orders 0,1,2 only. '
                       'No seed generator, higher-order branch, coefficient, cutoff or complete Excel program is identified.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def domain_run(input_path, destination):
    """Test fixed published guards as sole explanations, without inventing a seed body."""
    from constants import atan, function_bounds
    dest=Path(destination);dest.mkdir()
    a,b=atan(Q(1,5),64),atan(Q(1,239),32)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    maximum=Q(2**1024-2**971)
    logmax=function_bounds(maximum,maximum,'ln')
    def load(span,lattice):
        low,high=map(lattice.rn,span);require(low==high,'ambiguous universal constant load');return low
    constants={name:{kind:pair(load(span,lattice)) for kind,lattice in [('binary64',B64),('extended64',E64)]}
        for name,span in [('pi',pi),('logmax',logmax),('ribesl_exparg',(Q(709),Q(709))),('cody_i0_xmax',(Q('713.986'),Q('713.986')))]}
    rows=[];scope=Counter()
    for index,raw in enumerate(map(json.loads,Path(input_path).read_text().splitlines())):
        if raw['function']not in ('BESSELI','BESSELJ','BESSELY'):
            scope['other_function']+=1;continue
        try:
            require(raw['result_index']is None and len(raw['args'])==2,'shape')
            x,n=map(lambda v:Q(number(v)),raw['args']);require(n.denominator==1 and 0<=n<=100,'order')
            require(raw['expected'].startswith('0x')or raw['expected']=='error:Num','result')
            if raw['function']=='BESSELI':require(abs(x)>=500,'large_argument_scope')
            else:require(abs(x)>=2**24 and n<=4 and (raw['function']!='BESSELY'or x>0),'large_argument_scope')
        except (ValueError,TypeError,OverflowError)as error:
            scope[str(error)]+=1;continue
        rows.append(dict(trace_index=index,raw=raw));scope['scoped']+=1
    require(rows and any(r['raw']['function']=='BESSELI'for r in rows)and any(r['raw']['function']=='BESSELJ'for r in rows),'empty domain group')
    profiles=[]
    for family,forms in [('I',range(3)),('JY',range(5))]:
        models=('strict',)if family=='I'else MODELS
        for form,model,load_kind,inclusive in itertools.product(forms,models,('binary64','extended64'),(False,True)):
            profiles.append(dict(id=len(profiles),family=family,form=form,model=model,load=load_kind,inclusive=inclusive))
    (dest/'constants.json').write_text(json.dumps(constants,indent=2)+'\n')
    (dest/'scoped.json').write_text(json.dumps(rows)+'\n');(dest/'profiles.json').write_text(json.dumps(profiles)+'\n')
    def rejects(p,raw):
        x,n=map(lambda v:Q(number(v)),raw['args']);x=abs(x);n=int(n)
        if p['family']=='I':
            bound=Q(*map(int,constants[('ribesl_exparg','logmax','cody_i0_xmax')[p['form']]][p['load']]))
            return x>=bound if p['inclusive']else x>bound
        pi_value=Q(*map(int,constants['pi'][p['load']]))
        def phase(order,sine):
            multiplier=2*order+(-1 if sine else 1)
            numerator=operation('mul',p['model'],Q(multiplier),pi_value)
            correction=operation('div',p['model'],numerator,Q(4))
            return abs(operation('add',p['model'],x,-correction))
        if p['form']==0:arguments=[x]
        elif p['form']in(1,3):arguments=[phase(n,p['form']==3)]
        else:arguments=[phase(order,p['form']==4)for order in (0,1)]
        return any(v>=2**27 if p['inclusive']else v>2**27 for v in arguments)
    results=[]
    for p in profiles:
        observed=early=sole=0;first_early=first_sole=None
        for row in rows:
            raw=row['raw']
            if (raw['function']=='BESSELI')!=(p['family']=='I'):continue
            rejection=rejects(p,raw);error=raw['expected']=='error:Num';observed+=1
            if rejection and not error:
                early+=1
                if first_early is None:first_early=row['trace_index']
            if rejection!=error:
                sole+=1
                if first_sole is None:first_sole=row['trace_index']
        require(observed,'empty guard profile scope')
        results.append(dict(candidate_id=p['id'],observations=observed,
            rejects_recorded_numeric=early,sole_explanation_mismatches=sole,
            first_early_contradiction=first_early,first_sole_contradiction=first_sole))
    (dest/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    result=dict(complete=True,programs=len(profiles),observations=len(rows),scope=scope,
        sole_explanation_survivors=sum(r['sole_explanation_mismatches']==0 for r in results),
        early_guard_compatible=sum(r['rejects_recorded_numeric']==0 for r in results),
        source_sha256=digest(Path(__file__)),input_sha256={str(input_path):digest(input_path)},
        sources=['https://www.netlib.org/specfun/ribesl','https://www.netlib.org/specfun/i0',
                 'A&S 9.2 asymptotic phases; Machin identity; IEEE754',
                 'DnaCalc real_policy: measured circular trig limit 2^27'],
        limits='A guard rejecting a recorded numeric result is impossible as an early guard in this scope. A guard admitting a recorded NUM fails as its sole explanation, but may still precede a body that fails. No Bessel seed producer, fitted cutoff, new Excel answer or numerical repair is inferred.')
    (dest/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def normalized_i(input_path, producers_path, destination):
    """Miller's recurrence normalized to the existing public I0 producer.

    DLMF 10.29.1 and 10.74(iv). This is an original fixed program. Starting
    n+ceil(abs(x))+53 gives at least 53 tail steps above x; the ratio-map
    derivative there is bounded by 1/4. Rescaling uses binary64's exponent
    range, with no measured coefficient, threshold, depth or answer selector.
    """
    import math
    import struct
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    raw_rows = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    producers = [json.loads(line) for line in Path(producers_path).read_text().splitlines()]
    as_float = lambda value: struct.unpack('>d', bytes.fromhex(value[2:]))[0]
    i0 = {r['args'][0]: as_float(r['calc_func_prediction']) for r in producers
          if r['calc_func_prediction'].startswith('0x')}
    counts = Counter()
    decisions = []
    def flush(x):
        return math.copysign(0., x) if 0 < abs(x) < 2.**-1022 else x
    for raw in raw_rows:
        try:
            x, order = map(as_float, raw['args'])
        except (ValueError, TypeError):
            counts['outside_typed_scope'] += 1
            continue
        if raw['function'] != 'BESSELI' or raw['args'][0] not in i0 or not (2 <= order < 10000001) or not x:
            counts['outside_normalization_scope'] += 1
            continue
        n, x = int(order), abs(x)
        before, current, answer = 0., 1., 0.
        failed = False
        for k in range(n+math.ceil(x)+53, 0, -1):
            if k == n: answer = current
            product = flush(flush((2.*k)/x)*current)
            following = flush(before+product)
            if not math.isfinite(following):
                failed = True
                break
            before, current = current, following
            if current > 2.**512:
                before = flush(before*2.**-512)
                current = flush(current*2.**-512)
                answer = flush(answer*2.**-512)
        ratio = flush(answer/current) if not failed else 0.
        if not ratio:
            counts['unrepresented_ratio_use_series'] += 1
            continue
        value = flush(i0[raw['args'][0]]*ratio)
        if as_float(raw['args'][0]) < 0 and n % 2: value = -value
        actual = '0x'+struct.pack('>d', value).hex()
        matched, previous = actual == raw['expected'], raw['actual'] == raw['expected']
        counts['scoped'] += 1
        counts['matched'] += matched
        counts['gains'] += matched and not previous
        counts['losses'] += previous and not matched
        corpus = str(raw.get('id', '')).startswith('corpus:')
        counts['corpus_scoped' if corpus else 'bank_scoped'] += 1
        counts['corpus_gains' if corpus else 'bank_gains'] += matched and not previous
        decisions.append(dict(raw=raw, calculated=actual, matched=matched, previous_matched=previous))
    (destination/'records.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in decisions))
    result = dict(complete=True, program='backward Miller, I0*(I_n/I0), strict binary64 with FTZ',
        counts=counts, input_sha256={str(p): digest(p) for p in (input_path, producers_path, Path(__file__))},
        limits='Fixed mathematical program with the current cleared I0 producer. No fitted constant '
               'or Excel seed. Unrepresentable normalization ratios retain the existing positive-series '
               'path. Prototype comparisons are discovery, not new Excel readings or a complete Bessel implementation.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


def forward_groups(rows):
    """Positive normal publications, with J restricted to its forward regime."""
    grouped = defaultdict(dict)
    for index, row in enumerate(rows):
        a = row['args']
        if (row['function'] not in ('BESSELJ', 'BESSELY') or row['result_index'] is not None
                or len(a) != 2 or not all(isinstance(v, str) and v.startswith('0x') for v in a)
                or not row['expected'].startswith('0x')):
            continue
        x, n = map(number, a)
        if not (0 < x and n == int(n) and 0 <= n <= 32):
            continue
        if abs(Q(number(row['expected']))) < pow2(-1022):
            continue
        if row['function'] == 'BESSELJ' and n >= 2 and n >= x:
            continue
        key = context(row)
        if int(n) in grouped[key]:
            require(rows[grouped[key][int(n)]]['expected'] == row['expected'], 'conflicting forward publication')
        else:
            grouped[key][int(n)] = index
    return [dict(context=key, indices=[values[n] for n in sorted(values)])
            for key, values in sorted(grouped.items())
            if 0 in values and 1 in values and any(n >= 2 for n in values)]


def forward_profiles():
    # The first four forms have three sites. The reusable 2/x form has four.
    result = []
    for form in range(5):
        for models in itertools.product(MODELS, repeat=4 if form == 4 else 3):
            result.append(dict(id=len(result), form=form, models=list(models)))
    return result


def forward_step(x, n, previous, current, profile):
    m, form = profile['models'], profile['form']
    def op(kind, a, b, site):
        return operation(kind, m[site], a, b)
    coefficient = Q(2*n)
    if form == 0:
        term = op('div', op('mul', coefficient, current, 0), x, 1)
    elif form == 1:
        term = op('mul', coefficient, op('div', current, x, 0), 1)
    elif form == 2:
        term = op('mul', op('div', coefficient, x, 0), current, 1)
    elif form == 3:
        term = op('div', current, op('div', x, coefficient, 0), 1)
    else:
        reciprocal = op('div', Q(2), x, 0)
        term = op('mul', op('mul', Q(n), reciprocal, 1), current, 2)
    return op('add', term, -previous, len(m)-1)


def forward_bounds(rows, group, profile):
    observations = {int(number(rows[i]['args'][1])): i for i in group['indices']}
    x = Q(number(rows[observations[0]]['args'][0]))
    before = seed_bounds(Q(number(rows[observations[0]]['expected'])))
    current = seed_bounds(Q(number(rows[observations[1]]['expected'])))
    steps = []
    for n in range(1, max(observations)):
        # Positive coefficient: increasing in current, decreasing in previous.
        roots = (forward_step(x, n, before[1], current[0], profile),
                 forward_step(x, n, before[0], current[1], profile))
        require(roots[0] <= roots[1], 'forward monotonicity')
        following = roots
        step = dict(order=n+1, previous=list(map(pair, before)), current=list(map(pair, current)),
                    roots=list(map(pair, roots)), observation=observations.get(n+1))
        if n+1 in observations:
            target = Q(number(rows[observations[n+1]]['expected']))
            allowed = seed_bounds(target)
            following = max(roots[0], allowed[0]), min(roots[1], allowed[1])
            step.update(target=pair(target), allowed=list(map(pair, allowed)))
            if following[0] > following[1]:
                step['excluded'] = True
                steps.append(step)
                return dict(excluded=True, steps=steps)
        step.update(excluded=False, retained=list(map(pair, following)))
        steps.append(step)
        # This projection discards correlations and can only enlarge the set.
        before, current = current, following
    return dict(excluded=False, steps=steps)


def forward_run(input_path, destination):
    started = time.monotonic()
    dest = Path(destination)
    dest.mkdir()
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__), Path(__file__).with_name('intervals.py'))}
    rows = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    gs, ps = forward_groups(rows), forward_profiles()
    require(len(ps) == 1125 and gs, 'empty or incomplete forward grammar')
    for name, data in [('groups.json', gs), ('profiles.json', ps)]:
        (dest/name).write_text(json.dumps(data, indent=2)+'\n')
    counts = {}
    for function in ('BESSELJ', 'BESSELY'):
        scoped = [(i, g) for i, g in enumerate(gs) if rows[g['indices'][0]]['function'] == function]
        require(scoped, 'empty forward function scope')
        alive, executions = [], 0
        with (dest/(function+'-first-kills.jsonl')).open('x') as out:
            for p in ps:
                for order, (index, g) in enumerate(scoped):
                    result = forward_bounds(rows, g, p)
                    executions += 1
                    if result['excluded']:
                        out.write(json.dumps(dict(candidate_id=p['id'], order=order, group_index=index, proof=result))+'\n')
                        break
                else:
                    alive.append(p['id'])
        (dest/(function+'-necessary-survivors.json')).write_text(json.dumps(alive)+'\n')
        counts[function] = dict(profiles=len(ps), groups=len(scoped), executions=executions,
                                rejected=len(ps)-len(alive), necessary_survivors=len(alive))
    require(frozen == {p: digest(p) for p in frozen}, 'forward source/input changed')
    summary = dict(complete=True, counts=counts, seconds=time.monotonic()-started, input_sha256=frozen,
        claim='Fixed forward recurrence policies through all recorded orders up to32, with exact conservative Ext80 publication-preimage intersections at each recorded step.',
        limits='Y uses positive x; J additionally requires x>n for target orders>=2. Shared original source/environment and seed helper publications are assumptions. Correlations discarded after each step; a survivor need not have a consistent seed pair or producer. No seed, cutoff, coefficient or answer is fitted.')
    (dest/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


def miller_rows(rows):
    selected = []
    for index, row in enumerate(rows):
        a = row['args']
        if (row['function'] != 'BESSELJ' or row['result_index'] is not None or len(a) != 2
                or not all(isinstance(v, str) and v.startswith('0x') for v in a)
                or not row['expected'].startswith('0x')):
            continue
        x, order = map(number, a)
        n = int(order)
        if Q(1, 8) <= x <= n and 2 <= n <= 100 and abs(Q(number(row['expected']))) >= pow2(-1022):
            selected.append(index)
    return selected


def miller_profiles(expanded=False, sites=False):
    identity = 0
    for margin, last, models in itertools.product(range(65), (False, True), itertools.product(MODELS, repeat=6 if sites else 4)):
        for normalization in range(3 if expanded else 1):
            result = dict(id=identity, margin=margin, include_boundary=last, models=list(models))
            if expanded:
                result['normalization'] = normalization
            yield result
            identity += 1


@lru_cache(maxsize=8192)
def miller_prefix(x, n, margin, include_boundary, models):
    start = n+margin
    start += start % 2
    ratio_model, coefficient_model, product_model, subtract_model, sum_model = models
    ratio = operation('div', ratio_model, Q(2), x)
    following, current = Q(0), Q(1)
    selected = Q(1) if start == n else Q(0)
    even_sum = Q(1) if include_boundary else Q(0)
    scalings = 0
    for k in range(start, 0, -1):
        coefficient = operation('mul', coefficient_model, Q(k), ratio)
        previous = operation('add', subtract_model,
            operation('mul', product_model, coefficient, current), -following)
        following, current = current, previous
        if k-1 == n:
            selected = current
        if k-1 >= 2 and (k-1) % 2 == 0:
            even_sum = operation('add', sum_model, even_sum, current)
        if max(abs(following), abs(current), abs(even_sum)) > pow2(512):
            # Exact powers of two keep the trial sequence in its represented
            # range. No target or measured error enters this rescaling.
            following *= pow2(-512)
            current *= pow2(-512)
            selected *= pow2(-512)
            even_sum *= pow2(-512)
            def represented(value, model):
                fmt = B64 if model == 'strict' or model.startswith('stored') else E53 if model.endswith('pc53') else E64
                require(fmt.rn(value) == value, 'Miller rescaling loses a represented state')
            for value in (following, current, selected):
                represented(value, subtract_model)
            represented(even_sum, sum_model)
            scalings += 1
    return start, current, even_sum, selected, scalings


def miller_value(row, profile):
    """Trial recurrence normalized by J0+2*sum(J_even)=1 (DLMF10.12.4)."""
    x = Q(number(row['args'][0]))
    n = int(number(row['args'][1]))
    models = profile['models']
    if len(models) == 4:
        ratio_model, recurrence_model, sum_model, finish_model = models
        prefix = (ratio_model, recurrence_model, recurrence_model, recurrence_model, sum_model)
    else:
        require(len(models) == 6, 'Miller site policy count')
        prefix, sum_model, finish_model = tuple(models[:5]), models[4], models[5]
    start, current, even_sum, selected, scalings = miller_prefix(
        x, n, profile['margin'], profile['include_boundary'], prefix)
    form = profile.get('normalization', 0)
    if form == 0:
        denominator = operation('add', sum_model, current, operation('mul', sum_model, Q(2), even_sum))
    elif form == 1:
        total = operation('add', sum_model, current, even_sum)
        denominator = operation('add', sum_model, operation('mul', sum_model, Q(2), total), -current)
    else:
        denominator = operation('add', sum_model, operation('add', sum_model, current, even_sum), even_sum)
    require(denominator != 0, 'zero Miller normalizer')
    root = operation('div', finish_model, selected, denominator)
    return dict(start=start, scalings=scalings, root=pair(root), published=pair(B64.rn(root)))


def miller_run(input_path, destination, expanded=False, sites=False):
    started = time.monotonic()
    dest = Path(destination)
    dest.mkdir()
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__), Path(__file__).with_name('intervals.py'))}
    rows = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    scope = miller_rows(rows)
    require(scope, 'empty Miller scope')
    (dest/'scope.json').write_text(json.dumps(scope)+'\n')
    alive, executions, kill_counts = [], 0, Counter()
    import gzip
    with gzip.open(dest/'first-kills.jsonl.gz', 'xt', encoding='utf-8') as out:
        for p in miller_profiles(expanded, sites):
            for order, index in enumerate(scope):
                result = miller_value(rows[index], p)
                executions += 1
                if Q(*map(int, result['published'])) != Q(number(rows[index]['expected'])):
                    out.write(json.dumps(dict(candidate_id=p['id'], order=order, trace_index=index, proof=result))+'\n')
                    kill_counts[order] += 1
                    break
            else:
                alive.append(p['id'])
    (dest/'survivors.json').write_text(json.dumps(alive)+'\n')
    require(frozen == {p: digest(p) for p in frozen}, 'Miller source/input changed')
    programs = 65*2*5**(6 if sites else 4)*(3 if expanded else 1)
    summary = dict(complete=True, expanded=expanded, sites=sites, programs=programs, observations=len(scope), executions=executions,
        survivors=len(alive), first_kills=kill_counts, seconds=time.monotonic()-started, input_sha256=frozen,
        grammar='Even starting order ceil_even(n+margin), margin0..64; two truncated even-sum endpoints; five arithmetic models. Basic mode has ratio/recurrence/sum/final policies; sites mode splits recurrence coefficient/product/subtraction policies. Expanded mode has three fixed normalizer associations (t0+2*S, 2*(t0+S)-t0, (t0+S)+S). Cached prefixes depend only on operands and physical policies, never expected answers.',
        limits='Positive x in[1/8,n], truncated integer order2..100, normal numeric publications. DLMF10.6.1/10.12.4 identities; finite starting-order hypotheses only. No starting-depth rule, tolerance, seed coefficient or row-dependent selector is recovered or adopted. All observations are discovery.')
    (dest/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary))


def real_reference(row):
    """Real J/Y enclosures from DLMF10.2.2,10.8.1,10.17(iii), not Excel predictions."""
    from math import factorial, isqrt
    from constants import atan
    from intervals import log2floor
    x,order = map(lambda a:Q(number(a)),row['args'])
    n = int(order)
    require(n == order and 0 <= n <= 64 and row['function'] in ('BESSELJ','BESSELY'), 'real reference scope')
    negative = x < 0
    if row['function'] == 'BESSELJ':x = abs(x)
    require(x >= 0 and (x > 0 or row['function'] == 'BESSELJ'), 'real reference domain')
    def add(a,b):return a[0]+b[0],a[1]+b[1]
    def neg(a):return -a[1],-a[0]
    def mul(a,b):
        values = [u*v for u in a for v in b]
        return min(values),max(values)
    def inv(a):
        require(a[0]*a[1] > 0,'reference inverse crosses zero')
        return 1/a[1],1/a[0]
    a,b = atan(Q(1,5),64),atan(Q(1,239),32)
    pi = 16*a[0]-4*b[1],16*a[1]-4*b[0]
    def logarithm(q):
        exponent = log2floor(q)
        reduced = q/pow2(exponent)
        def atanh_log(t):
            z = (t-1)/(t+1)
            total = sum((2*z**(2*k+1)/Q(2*k+1) for k in range(192)),Q(0))
            tail = 2*abs(z)**385/Q(385)/(1-z*z)
            return total-tail,total+tail
        return add(atanh_log(reduced),mul((Q(exponent),Q(exponent)),atanh_log(Q(2))))
    if x <= 32:
        require(x == 0 or x >= Q(1,1024), 'tiny argument outside reference scope')
        t = (x/2)**n/factorial(n)
        total,weighted = Q(0),Q(0)
        h,higher = Q(0),sum((Q(1,j)for j in range(1,n+1)),Q(0))
        for k in range(512):
            total += t
            weighted += (h+higher)*t
            next_term = -t*(x*x/4)/((k+1)*(n+k+1))
            next_weight = h+higher+Q(1,k+1)+Q(1,n+k+1)
            if k > x and abs(next_term)*max(Q(1),2*next_weight) < pow2(-192):
                j = total-abs(next_term),total+abs(next_term)
                w = weighted-2*abs(next_term)*next_weight,weighted+2*abs(next_term)*next_weight
                break
            t = next_term
            h += Q(1,k+1)
            higher += Q(1,n+k+1)
        else:raise ValueError('series enclosure did not converge')
        if row['function'] == 'BESSELJ':value=j
        else:
            # DLMF5.2.3 prints these twenty digits as the prefix of gamma.
            gamma = Q('0.57721566490153286060'),Q('0.57721566490153286061')
            normalizer = mul((Q(2),Q(2)),inv(pi))
            value = mul(normalizer,mul(add(logarithm(x/2),gamma),j))
            finite = (x/2)**(-n)*sum((Q(factorial(n-i-1),factorial(i))*(x*x/4)**i for i in range(n)),Q(0))
            value = add(value,neg(mul(inv(pi),add((finite,finite),w))))
        method = 'defining_series'
        terms = k+1
    else:
        terms = max(16,(n+1)//2)
        coefficients = [Q(1)]
        for k in range(1,2*terms+2):
            coefficients.append(coefficients[-1]*Q(4*n*n-(2*k-1)**2,8*k))
        def expansion(parity):
            total = sum(((-1)**k*coefficients[2*k+parity]/x**(2*k+parity)for k in range(terms)),Q(0))
            tail = (-1)**terms*coefficients[2*terms+parity]/x**(2*terms+parity)
            return min(total,total+tail),max(total,total+tail)
        p,q = expansion(0),expansion(1)
        phase = add((x,x),neg(mul((Q(2*n+1,4),Q(2*n+1,4)),pi)))
        turns = int(((phase[0]+phase[1])/(4*((pi[0]+pi[1])/2))+Q(1,2))//1)
        reduced = add(phase,neg(mul((Q(2*turns),Q(2*turns)),pi)))
        midpoint = sum(reduced)/2
        # A dyadic Taylor center avoids huge odd-denominator powers. Its exact
        # displacement is included in the radius, so this only widens a proof.
        center = (midpoint/pow2(-192)//1)*pow2(-192)
        radius = (reduced[1]-reduced[0])/2+abs(center-midpoint)
        require(abs(center) < 4,'trigonometric reduction outside Taylor bound')
        def trig(parity):
            total = sum(((-1)**k*center**(2*k+parity)/factorial(2*k+parity)for k in range(64)),Q(0))
            error = abs(center)**(128+parity)/factorial(128+parity)+radius
            return total-error,total+error
        sine,cosine = trig(1),trig(0)
        square = mul((Q(2)/x,Q(2)/x),inv(pi))
        unit = 1<<256
        floor = lambda q: Q(isqrt(q.numerator*unit*unit//q.denominator),unit)
        amplitude = floor(square[0]),floor(square[1])+Q(1,unit)
        body = add(mul(cosine,p),neg(mul(sine,q))) if row['function'] == 'BESSELJ' else add(mul(sine,p),mul(cosine,q))
        value = mul(amplitude,body)
        method = 'hankel_remainder_bound'
    if negative and n%2:value = neg(value)
    lower,upper = E64.floor(value[0]),E64.ceil(value[1])
    require(lower is not None and upper is not None,'reference output outside Ext80')
    return dict(method=method,terms=terms,enclosure=list(map(pair,(lower,upper))))


def regime_run(input_path,destination):
    dest = Path(destination);dest.mkdir()
    source = digest(Path(__file__))
    rows = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    counts,regions,excluded = Counter(),defaultdict(Counter),[]
    cache = {}
    with (dest/'reference-bounds.jsonl').open('x') as out:
        for index,row in enumerate(rows):
            if row['function'] not in ('BESSELJ','BESSELY'):
                counts['other_function'] += 1;continue
            if row['result_index'] is not None or len(row['args']) != 2 or not all(isinstance(v,str)and v.startswith('0x')for v in row['args']):
                counts['surface'] += 1;continue
            if not row['expected'].startswith('0x'):
                counts['domain_obligation'] += 1;continue
            try:
                key = row['function'],tuple(row['args'])
                if key not in cache:cache[key] = real_reference(row)
                proof = cache[key]
            except ValueError as error:
                counts['outside_reference_scope'] += 1;excluded.append(dict(trace_index=index,reason=str(error)));continue
            lo,hi = [Q(*map(int,v))for v in proof['enclosure']]
            published = B64.rn(lo),B64.rn(hi)
            target = Q(number(row['expected']))
            impossible = not published[0] <= target <= published[1]
            x,n = map(number,row['args'])
            region = ('seed' if n <= 1 else 'above_order' if abs(x) > n else 'at_or_below_order')+'_'+('small' if abs(x) <= 32 else 'large')
            regions[row['function']+':'+region]['observations'] += 1
            regions[row['function']+':'+region]['outside_real_rounding_bound'] += impossible
            counts['scoped'] += 1;counts['outside_real_rounding_bound'] += impossible
            out.write(json.dumps(dict(trace_index=index,id=row['id'],function=row['function'],args=row['args'],
                expected=row['expected'],region=region,proof=proof,outside_real_rounding_bound=impossible))+'\n')
    # ubs:ignore — public source-integrity hashes, never secret authentication material.
    require(source == digest(Path(__file__)),'reference source changed during run')
    result = dict(complete=True,counts=counts,regions=regions,excluded=excluded,
        source_sha256=source,input_sha256=digest(input_path),bounds_sha256=digest(dest/'reference-bounds.jsonl'),
        sources=['https://dlmf.nist.gov/10.2.E2','https://dlmf.nist.gov/10.8.E1','https://dlmf.nist.gov/10.17.iii','https://dlmf.nist.gov/5.2.E3'],
        claim='Rational enclosures exclude correctly rounded real J/Y kernels at the listed recorded inputs. '
              'The series and Hankel labels name reference mathematics, not identified Excel implementation branches.',
        limits='Integer orders0..64; argument zero for J or magnitude>=1/1024. '
               'These constraints do not exclude inaccurately truncated/rounded series, approximate seed recurrences or asymptotic approximations. '
               'No complete Excel program, cutoff, seed or coefficient is selected; no product gain follows.')
    with (dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='regimes':regime_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='miller-sites':miller_run(*sys.argv[2:], expanded=True, sites=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='miller-forms':miller_run(*sys.argv[2:], expanded=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='miller':miller_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='forward':forward_run(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='normalized-i':normalized_i(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='domains':domain_run(*sys.argv[2:]);raise SystemExit(0)
    parser=argparse.ArgumentParser();parser.add_argument('input');parser.add_argument('destination')
    args=parser.parse_args();run(args.input,args.destination)
