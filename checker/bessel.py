"""Independent exact proof of necessary Bessel order-two recurrence bounds."""
from collections import defaultdict,Counter
from fractions import Fraction as R
import itertools
import json
from pathlib import Path
import sys
import time
from functools import lru_cache

from exact import BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P53, Span, binary64, operation, two
from check import records, require, digest

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def decode_pair(pair):
    require(len(pair)==2 and all(isinstance(x,str)for x in pair),'rational certificate encoding')
    return R(*map(int,pair))


def groups(rows):
    result=defaultdict(dict)
    for index,r in enumerate(rows):
        args=r['args']
        if len(args)!=2 or r['result_index']is not None or not all(isinstance(x,str)and x.startswith('0x')for x in args):continue
        x,n=map(binary64,args)
        if x<=0 or n not in (0,1,2)or not r['expected'].startswith('0x'):continue
        if abs(binary64(r['expected']))<two(-1022):continue
        key=json.dumps([r['function'],args[0],r['reading'],r['source']],sort_keys=True,separators=(',',':'))
        if int(n)in result[key]:
            require(rows[result[key][int(n)]]['expected']==r['expected'],'conflicting source/environment')
        else:result[key][int(n)]=index
    return [dict(context=k,indices=[v[n]for n in (0,1,2)])for k,v in sorted(result.items())if set(v)=={0,1,2}]


def seed_bounds(target):
    span=B.pullback(Span(target,target))
    lo,hi=E.coordinates(span)
    return E.value(lo),E.value(hi)


def calculate(profile,function,x,seed0,seed1):
    m=profile['models'];form=profile['shape']
    if form==0:
        t=operation(R(2),seed1,'mul',m[0]);t=operation(t,x,'div',m[1])
    elif form==1:
        t=operation(seed1,x,'div',m[0]);t=operation(R(2),t,'mul',m[1])
    elif form==2:
        t=operation(R(2),x,'div',m[0]);t=operation(t,seed1,'mul',m[1])
    else:
        d=operation(x,R(2),'div',m[0]);t=operation(seed1,d,'div',m[1])
    if function in ('BESSELJ','BESSELY'):return operation(t,-seed0,'add',m[2])
    return operation(seed0,t if function=='BESSELK'else-t,'add',m[2])


def proof(rows,g,p,certificate=None):
    zero,one,target=[rows[i]for i in g['indices']]
    function=zero['function'];x=binary64(zero['args'][0]);y=binary64(target['expected'])
    z=seed_bounds(binary64(zero['expected']));o=seed_bounds(binary64(one['expected']))
    if function in ('BESSELJ','BESSELY'):corners=((z[1],o[0]),(z[0],o[1]))
    elif function=='BESSELI':corners=((z[0],o[1]),(z[1],o[0]))
    else:corners=((z[0],o[0]),(z[1],o[1]))
    roots=[calculate(p,function,x,*pair)for pair in corners]
    published=[B.round(root)for root in roots]
    require(roots[0]<=roots[1]and published[0]<=published[1],'monotone orientation failed')
    excluded=not(published[0]<=y<=published[1])
    if certificate is not None:
        for field,expected in [('seed0',z),('seed1',o),('root_bounds',roots),('published_bounds',published)]:
            require(list(map(decode_pair,certificate[field]))==list(expected),'incorrect '+field)
        require([[decode_pair(a),decode_pair(b)]for a,b in certificate['corner_seeds']]==[list(c)for c in corners],'incorrect extremal seed pairing')
        require(decode_pair(certificate['target'])==y and certificate['excluded']==excluded,'changed target or exclusion')
    return excluded


def check(input_path,directory,output):
    started=time.monotonic();directory=Path(directory)
    source=Path(__file__);hashes={str(p):digest(p)for p in (source,source.with_name('exact.py'),source.with_name('check.py'))}
    rows=list(records(input_path));gs=groups(rows)
    require(gs==json.loads((directory/'groups.json').read_text()),'changed group selection')
    ps=[dict(id=i,shape=s,models=list(m))for i,(s,m)in enumerate(itertools.product(range(4),itertools.product(MODELS,repeat=3)))]
    require(ps==json.loads((directory/'profiles.json').read_text())and len(ps)==500,'changed grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial search')
    result={};counts=Counter()
    for function in ('BESSELJ','BESSELY','BESSELI','BESSELK'):
        cohort=[(i,g)for i,g in enumerate(gs)if rows[g['indices'][0]]['function']==function]
        killed=set();executions=0
        for row in records(directory/f'{function}-first-kills.jsonl'):
            identity=row['candidate_id'];order=row['order']
            require(identity not in killed and 0<=identity<len(ps),'duplicate or unknown profile')
            require(0<=order<len(cohort)and row['group_index']==cohort[order][0],'changed first-kill group')
            for previous in range(order):
                require(not proof(rows,cohort[previous][1],ps[identity]),'not first rejection');executions+=1
            require(proof(rows,cohort[order][1],ps[identity],row['proof']),'invalid rejection');executions+=1
            killed.add(identity);counts['independently_checked_seed_box_exclusions']+=1
        alive=json.loads((directory/f'{function}-necessary-survivors.json').read_text())
        require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(500)),'incomplete profile accounting')
        for identity in alive:
            for _,g in cohort:
                require(not proof(rows,g,ps[identity]),'excluded necessary survivor');executions+=1
        actual=dict(profiles=500,groups=len(cohort),executions=executions,rejected=len(killed),necessary_survivors=len(alive))
        require(summary['counts'][function]==actual,'incorrect summary')
        result[function]=actual
        counts['independently_checked_seed_boxes']+=executions
        counts['full_width_corner_results']+=2*executions
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during run')  # ubs:ignore — public integrity hashes, not secrets
    receipt=dict(complete=True,counts=counts,functions=result,seconds=time.monotonic()-started,
                 source_sha256=hashes,input_sha256={str(p):digest(p)for p in (Path(input_path),*sorted(directory.iterdir()))if p.is_file()},
                 claim='Exact monotone corner proof over all Ext80 seed pairs publishing to recorded orders zero and one. Necessary constraints only; no seed pair is fitted or selected.',
                 limits='Positive x, normal nonzero publications, fixed specified consumers and same original source/environment. A necessary survivor is not an executable complete program.')
    with open(output,'x')as out:json.dump(receipt,out,indent=2);out.write('\n')
    print(json.dumps(receipt))


def check_domains(input_path,directory,output):
    from reference import log_bounds,LIB
    directory=Path(directory)
    frozen={str(p):digest(p)for p in (Path(input_path),*directory.iterdir(),Path(__file__),
            Path(__file__).with_name('reference.py'),Path(__file__).with_name('exact.py'))if p.is_file()}
    def atan_bounds(denominator,terms):
        total=sum((R((-1)**i,(2*i+1)*denominator**(2*i+1))for i in range(terms)),R(0))
        next_term=R((-1)**terms,(2*terms+1)*denominator**(2*terms+1))
        return min(total,total+next_term),max(total,total+next_term)
    a,b=atan_bounds(5,80),atan_bounds(239,40)
    pi=(16*a[0]-4*b[1],16*a[1]-4*b[0])
    logmax=log_bounds(0x7fefffffffffffff)
    constants={}
    for name,span in [('pi',pi),('logmax',logmax),('ribesl_exparg',(R(709),R(709))),('cody_i0_xmax',(R('713.986'),R('713.986')))]:
        constants[name]={}
        for kind,lattice in [('binary64',B),('extended64',E)]:
            lo,hi=map(lattice.round,span);require(lo==hi,'ambiguous independent constant load')
            constants[name][kind]=[str(lo.numerator),str(lo.denominator)]
    require(json.loads((directory/'constants.json').read_text())==constants,'incorrect domain constants')
    rows=[];scope=Counter()
    for index,raw in enumerate(records(input_path)):
        if raw['function']not in ('BESSELI','BESSELJ','BESSELY'):scope['other_function']+=1;continue
        try:
            require(raw['result_index']is None and len(raw['args'])==2,'shape')
            x,n=map(binary64,raw['args']);require(n.denominator==1 and 0<=n<=100,'order')
            require(raw['expected'].startswith('0x')or raw['expected']=='error:Num','result')
            require(abs(x)>=500 if raw['function']=='BESSELI'else abs(x)>=two(24)and n<=4 and(raw['function']!='BESSELY'or x>0),'large_argument_scope')
        except(ValueError,TypeError,OverflowError)as error:scope[str(error)]+=1;continue
        rows.append(dict(trace_index=index,raw=raw));scope['scoped']+=1
    require(rows and any(r['raw']['function']=='BESSELI'for r in rows)and any(r['raw']['function']=='BESSELJ'for r in rows),'empty domain group')
    require(json.loads((directory/'scoped.json').read_text())==rows,'changed domain scope')
    ps=[]
    for family,forms in [('I',range(3)),('JY',range(5))]:
        for form,model,load,inclusive in itertools.product(forms,('strict',)if family=='I'else MODELS,('binary64','extended64'),(False,True)):
            ps.append(dict(id=len(ps),family=family,form=form,model=model,load=load,inclusive=inclusive))
    require(json.loads((directory/'profiles.json').read_text())==ps,'changed domain grammar')
    def rejects(p,raw):
        x,n=map(binary64,raw['args']);x=abs(x);n=int(n)
        if p['family']=='I':
            value=decode_pair(constants[('ribesl_exparg','logmax','cody_i0_xmax')[p['form']]][p['load']])
            return x>=value if p['inclusive']else x>value
        pi_value=decode_pair(constants['pi'][p['load']])
        orders=[n]if p['form']in(1,3)else[0,1]
        angles=[]
        for order in orders:
            multiple=2*order+(-1 if p['form']in(3,4)else 1)
            product=operation(R(multiple),pi_value,'mul',p['model'])
            offset=operation(product,R(4),'div',p['model'])
            angles.append(abs(operation(x,-offset,'add',p['model'])))
        if p['form']==0:angles=[x]
        return any(v>=two(27)if p['inclusive']else v>two(27)for v in angles)
    expected=[]
    for p in ps:
        hits=[];early=[];observed=0
        for row in rows:
            raw=row['raw']
            if(raw['function']=='BESSELI')!=(p['family']=='I'):continue
            refusal=rejects(p,raw);numeric=raw['expected']!='error:Num';observed+=1
            if refusal and numeric:early.append(row['trace_index'])
            if refusal==numeric:hits.append(row['trace_index'])
        require(observed,'empty domain profile')
        expected.append(dict(candidate_id=p['id'],observations=observed,rejects_recorded_numeric=len(early),
            sole_explanation_mismatches=len(hits),first_early_contradiction=early[0]if early else None,
            first_sole_contradiction=hits[0]if hits else None))
    require(json.loads((directory/'results.json').read_text())==expected,'incorrect domain contradiction')
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['complete']and summary['programs']==len(ps)and summary['observations']==len(rows)and summary['scope']==scope and
            summary['sole_explanation_survivors']==sum(r['sole_explanation_mismatches']==0 for r in expected)and
            summary['early_guard_compatible']==sum(r['rejects_recorded_numeric']==0 for r in expected),'wrong domain summary')
    require(all(digest(p)==h for p,h in frozen.items()),'domain evidence/source changed')
    result=dict(complete=True,programs=len(ps),observations=len(rows),sole_explanation_survivors=summary['sole_explanation_survivors'],
        early_guard_compatible=summary['early_guard_compatible'],input_sha256=frozen,mpfr_version=LIB.mpfr_get_version().decode(),
        limits='Independent rational arithmetic, Machin brackets and directed MPFR logmax loads. Early-guard impossibilities and incomplete sole explanations are counted separately. No numeric Bessel body, fitted cutoff or domain repair is identified.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def forward_cohorts(rows):
    by_context = defaultdict(dict)
    for index, row in enumerate(rows):
        a = row['args']
        if (row['function'] not in ('BESSELJ', 'BESSELY') or row['result_index'] is not None
                or len(a) != 2 or not all(isinstance(v, str) and v.startswith('0x') for v in a)
                or not row['expected'].startswith('0x')):
            continue
        x, n = map(binary64, a)
        if x <= 0 or n.denominator != 1 or not 0 <= n <= 32:
            continue
        if abs(binary64(row['expected'])) < two(-1022):
            continue
        if row['function'] == 'BESSELJ' and n >= 2 and x <= n:
            continue
        key = json.dumps([row['function'], a[0], row['reading'], row['source']], sort_keys=True, separators=(',', ':'))
        n = int(n)
        if n in by_context[key]:
            require(rows[by_context[key][n]]['expected'] == row['expected'], 'conflicting forward answer')
        else:
            by_context[key][n] = index
    return [dict(context=k, indices=[v[n] for n in sorted(v)])
            for k, v in sorted(by_context.items()) if 0 in v and 1 in v and max(v) >= 2]


def forward_proof(rows, group, profile):
    observations = {int(binary64(rows[i]['args'][1])): i for i in group['indices']}
    x = binary64(rows[observations[0]]['args'][0])
    previous = seed_bounds(binary64(rows[observations[0]]['expected']))
    current = seed_bounds(binary64(rows[observations[1]]['expected']))
    model = profile['models']
    def apply(kind, left, right, site):
        values = [operation(a, b, kind, model[site]) for a, b in itertools.product(left, right)]
        return min(values), max(values)
    pair = lambda q: [str(q.numerator), str(q.denominator)]
    encoded = lambda bounds: list(map(pair, bounds))
    steps = []
    for order in range(2, max(observations)+1):
        factor = (R(2*(order-1)),)*2
        argument = (x, x)
        form = profile['form']
        if form == 0:
            term = apply('div', apply('mul', factor, current, 0), argument, 1)
        elif form == 1:
            term = apply('mul', factor, apply('div', current, argument, 0), 1)
        elif form == 2:
            term = apply('mul', apply('div', factor, argument, 0), current, 1)
        elif form == 3:
            term = apply('div', current, apply('div', argument, factor, 0), 1)
        else:
            reciprocal = apply('div', (R(2),)*2, argument, 0)
            coefficient = apply('mul', (R(order-1),)*2, reciprocal, 1)
            term = apply('mul', coefficient, current, 2)
        roots = apply('add', term, (-previous[1], -previous[0]), len(model)-1)
        step = dict(order=order, previous=encoded(previous), current=encoded(current),
                    roots=encoded(roots), observation=observations.get(order))
        retained = roots
        if order in observations:
            target = binary64(rows[observations[order]]['expected'])
            allowed = seed_bounds(target)
            step.update(target=pair(target), allowed=encoded(allowed))
            retained = max(roots[0], allowed[0]), min(roots[1], allowed[1])
            if retained[0] > retained[1]:
                step['excluded'] = True
                return dict(excluded=True, steps=steps+[step])
        step.update(excluded=False, retained=encoded(retained))
        steps.append(step)
        previous, current = current, retained
    return dict(excluded=False, steps=steps)


def check_forward(input_path, directory, output):
    started = time.monotonic()
    directory = Path(directory)
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__),
        Path(__file__).with_name('exact.py'), Path(__file__).with_name('check.py'), *sorted(directory.iterdir())) if p.is_file()}
    rows = list(records(input_path))
    gs = forward_cohorts(rows)
    require(gs and gs == json.loads((directory/'groups.json').read_text()), 'empty or changed higher-order scope')
    profiles = []
    for form in range(5):
        for choices in itertools.product(MODELS, repeat=4 if form == 4 else 3):
            profiles.append(dict(id=len(profiles), form=form, models=list(choices)))
    require(len(profiles) == 1125 and profiles == json.loads((directory/'profiles.json').read_text()), 'changed higher-order grammar')
    summary = json.loads((directory/'summary.json').read_text())
    require(summary['complete'] is True, 'partial forward run')
    counts, functions = Counter(), {}
    for function in ('BESSELJ', 'BESSELY'):
        scoped = [(i, g) for i, g in enumerate(gs) if rows[g['indices'][0]]['function'] == function]
        require(scoped, 'empty higher-order function cohort')
        killed = set()
        executions = 0
        for row in records(directory/(function+'-first-kills.jsonl')):
            identity, order = row['candidate_id'], row['order']
            require(identity not in killed and 0 <= identity < len(profiles), 'duplicate/unknown forward profile')
            require(0 <= order < len(scoped) and row['group_index'] == scoped[order][0], 'wrong forward first kill')
            for earlier in range(order):
                require(not forward_proof(rows, scoped[earlier][1], profiles[identity])['excluded'], 'earlier contradiction omitted')
                executions += 1
            independently = forward_proof(rows, scoped[order][1], profiles[identity])
            require(independently == row['proof'] and independently['excluded'], 'incorrect forward proof')
            executions += 1
            killed.add(identity)
        alive = json.loads((directory/(function+'-necessary-survivors.json')).read_text())
        require(len(alive) == len(set(alive)) and killed.isdisjoint(alive)
                and killed | set(alive) == set(range(len(profiles))), 'incomplete forward coverage')
        for identity in alive:
            for _, group in scoped:
                require(not forward_proof(rows, group, profiles[identity])['excluded'], 'invalid forward survivor')
                executions += 1
        totals = dict(profiles=len(profiles), groups=len(scoped), executions=executions,
                      rejected=len(killed), necessary_survivors=len(alive))
        require(summary['counts'][function] == totals, 'incorrect forward totals')
        functions[function] = totals
        counts['independent_cohort_replays'] += executions
    require(frozen == {p: digest(p) for p in frozen}, 'forward evidence changed during check')
    receipt = dict(complete=True, functions=functions, counts=counts, seconds=time.monotonic()-started,
        input_sha256=frozen, limits='Independent exact lattice arithmetic and interval-corner replay. Publication constraints are intersected at every recorded order. Necessary conditions only; no shared seed producer or complete Bessel program.')
    with open(output, 'x') as out:
        json.dump(receipt, out, indent=2)
        out.write('\n')
    print(json.dumps(receipt))


def miller_scope(rows):
    selected = []
    for index, row in enumerate(rows):
        args = row['args']
        if (row['function'] != 'BESSELJ' or row['result_index'] is not None or len(args) != 2
                or not all(isinstance(v, str) and v.startswith('0x') for v in args)
                or not row['expected'].startswith('0x')):
            continue
        x, order = map(binary64, args)
        n = int(order)
        if R(1, 8) <= x <= n and 2 <= n <= 100 and abs(binary64(row['expected'])) >= two(-1022):
            selected.append(index)
    return selected


@lru_cache(maxsize=131072)
def miller_op(kind, model, a, b):
    return operation(a, b, kind, model)


@lru_cache(maxsize=8192)
def miller_state(x, n, margin, boundary, policies):
    top = 2*((n+margin+1)//2)
    ratio_policy, coefficient_policy, multiplication_policy, subtraction_policy, reduction_policy = policies
    ratio = miller_op('div', ratio_policy, R(2), x)
    state = [R(1), R(0)]  # t_k, t_(k+1)
    target_term = R(1) if top == n else R(0)
    total = R(int(boundary))
    shifts = 0
    for order in reversed(range(top)):
        factor = miller_op('mul', coefficient_policy, R(order+1), ratio)
        product = miller_op('mul', multiplication_policy, factor, state[0])
        computed = miller_op('add', subtraction_policy, product, -state[1])
        state = [computed, state[0]]
        if order == n:
            target_term = state[0]
        if order != 0 and order % 2 == 0:
            total = miller_op('add', reduction_policy, total, state[0])
        if any(abs(value) > two(512) for value in (*state, total)):
            state = [value*two(-512) for value in state]
            total *= two(-512)
            target_term *= two(-512)
            def represented(value, policy):
                fmt = B if policy == 'strict' or policy.startswith('stored') else E if policy.endswith('pc64') else P53
                require(fmt.round(value) == value, 'Miller rescaling loses a represented state')
            for value in (*state, target_term):
                represented(value, subtraction_policy)
            represented(total, reduction_policy)
            shifts += 1
    return top, state[0], total, target_term, shifts


def miller_result(row, profile):
    x = binary64(row['args'][0])
    n = int(binary64(row['args'][1]))
    policies = profile['models']
    if len(policies) == 4:
        ratio, recurrence, reduction_policy, finish_policy = policies
        prefix = (ratio, recurrence, recurrence, recurrence, reduction_policy)
    else:
        require(len(policies) == 6, 'wrong independent Miller site count')
        prefix, reduction_policy, finish_policy = tuple(policies[:5]), policies[4], policies[5]
    top, first, total, target_term, shifts = miller_state(
        x, n, profile['margin'], profile['include_boundary'], prefix)
    form = profile.get('normalization', 0)
    if form == 0:
        doubled = miller_op('mul', reduction_policy, R(2), total)
        normalizer = miller_op('add', reduction_policy, first, doubled)
    elif form == 1:
        added = miller_op('add', reduction_policy, first, total)
        doubled = miller_op('mul', reduction_policy, R(2), added)
        normalizer = miller_op('add', reduction_policy, doubled, -first)
    else:
        added = miller_op('add', reduction_policy, first, total)
        normalizer = miller_op('add', reduction_policy, added, total)
    require(normalizer != 0, 'zero independent Miller normalizer')
    root = miller_op('div', finish_policy, target_term, normalizer)
    pair = lambda v: [str(v.numerator), str(v.denominator)]
    return dict(start=top, scalings=shifts, root=pair(root), published=pair(B.round(root)))


def miller_profile(identity, expanded=False, sites=False):
    width, associations = (6 if sites else 4), (3 if expanded else 1)
    association, rest = identity % associations, identity // associations
    model_index, outer = rest % (5**width), rest // (5**width)
    policies = [MODELS[(model_index//5**(width-1-i)) % 5] for i in range(width)]
    result = dict(id=identity, margin=outer//2, include_boundary=bool(outer%2), models=policies)
    if expanded:
        result['normalization'] = association
    return result


def check_miller(input_path, directory, output, expanded=False, sites=False):
    import gzip
    started = time.monotonic()
    directory = Path(directory)
    frozen = {str(p): digest(p) for p in (Path(input_path), Path(__file__),
        Path(__file__).with_name('exact.py'), Path(__file__).with_name('check.py'), *sorted(directory.iterdir())) if p.is_file()}
    rows = list(records(input_path))
    selected = miller_scope(rows)
    require(selected and selected == json.loads((directory/'scope.json').read_text()), 'changed Miller obligations')
    programs = 65*2*5**(6 if sites else 4)*(3 if expanded else 1)
    executions, killed, first_kills = 0, bytearray(programs), Counter()
    killed_count = 0
    with gzip.open(directory/'first-kills.jsonl.gz', 'rt', encoding='utf-8') as source:
        for line in source:
            row = json.loads(line)
            identity, order = row['candidate_id'], row['order']
            require(0 <= identity < programs and not killed[identity], 'duplicate or unknown Miller program')
            require(0 <= order < len(selected) and row['trace_index'] == selected[order], 'changed Miller first kill')
            for earlier in range(order+1):
                independently = miller_result(rows[selected[earlier]], miller_profile(identity, expanded, sites))
                expected = binary64(rows[selected[earlier]]['expected'])
                matched = decode_pair(independently['published']) == expected
                executions += 1
                if earlier == order:
                    require(independently == row['proof'] and not matched, 'incorrect Miller contradiction')
                else:
                    require(matched, 'earlier Miller contradiction omitted')
            killed[identity] = 1
            killed_count += 1
            first_kills[str(order)] += 1
    alive = json.loads((directory/'survivors.json').read_text())
    require(len(alive) == len(set(alive)) and all(0 <= i < programs and not killed[i] for i in alive)
            and killed_count+len(alive) == programs, 'incomplete Miller accounting')
    for identity in alive:
        for index in selected:
            require(decode_pair(miller_result(rows[index], miller_profile(identity, expanded, sites))['published']) == binary64(rows[index]['expected']), 'invalid Miller survivor')
            executions += 1
    summary = json.loads((directory/'summary.json').read_text())
    require(summary['complete'] is True and summary['expanded'] == expanded and summary.get('sites', False) == sites and summary['programs'] == programs
            and summary['observations'] == len(selected) and summary['executions'] == executions
            and summary['survivors'] == len(alive) and summary['first_kills'] == first_kills, 'incorrect Miller totals')
    require(frozen == {p: digest(p) for p in frozen}, 'Miller evidence changed during check')
    receipt = dict(complete=True, expanded=expanded, sites=sites, programs=programs, observations=len(selected), executions=executions,
        survivors=len(alive), seconds=time.monotonic()-started, input_sha256=frozen,
        limits='Independent rational recurrence, normalization and first-kill coverage for the specified finite starting-depth family. No fitted depth rule, elementary leaf, complete Excel algorithm or product improvement.')
    with open(output, 'x') as out:
        json.dump(receipt, out, indent=2)
        out.write('\n')
    print(json.dumps(receipt))


def asymptotic_inventory(directory):
    directory=Path(directory)
    summary=json.loads((directory/'summary.json').read_text())
    grammar=json.loads((directory/'grammar.json').read_text())
    survivors=json.loads((directory/'BESSELJ-survivors.json').read_text())
    axes=[['phase',10],['scale',10],['body',5],['recurrence',5],['trig',3],['terms',3],['association',3]]
    require(grammar['programs']==67500 and grammar['axes_in_loop_order']==axes,'asymptotic axes differ')
    counts=summary['counts']
    require(len(counts)==2 and {r['function']for r in counts}=={'BESSELJ','BESSELY'},'asymptotic function inventory differs')
    j=next(r for r in counts if r['function']=='BESSELJ')
    require(j['programs']==67500 and len(set(survivors))==len(survivors)==j['survivors']>0
            and all(isinstance(i,int)and 0<=i<67500 for i in survivors),'bad asymptotic survivor inventory')
    return summary,grammar,survivors


def check_asymptotic_survivors(input_path,tables_path,directory,executable,output):
    """Exact J survivor replay with separately calculated registered trig leaves."""
    import csv,struct,subprocess
    from functools import lru_cache
    directory=Path(directory);executable=Path(executable).resolve();output=Path(output)
    rows=list(records(input_path));tables=json.loads(Path(tables_path).read_text())
    summary,grammar,survivors=asymptotic_inventory(directory)
    frozen={str(p):digest(p)for p in (Path(input_path),Path(tables_path),directory/'summary.json',directory/'grammar.json',
        directory/'BESSELJ-survivors.json',executable,Path(__file__),Path(__file__).with_name('exact.py'))}
    def decode(identity):
        result={}
        for key,size in reversed((('phase',10),('scale',10),('body',5),('recurrence',5),('trig',3),('terms',3),('association',3))):
            result[key]=identity%size;identity//=size
        require(identity==0,'out-of-range asymptotic program');return result
    programs=[decode(i)for i in survivors]
    require(all(p['trig']==0 for p in programs),'this independent replay requires registered trig leaves')
    # Independently reconstruct each universal-constant quantizer.
    def atan(d,n):return sum((R((-1)**k,(2*k+1)*d**(2*k+1))for k in range(n)),R(0))
    pi=16*atan(5,128)-4*atan(239,64)
    pi_error=R(16,257*5**257)+R(4,129*239**129)
    def rounded_decimal(value,digits,mode):
        scale=10**digits;q=value*scale;k=q.numerator//q.denominator
        if mode=='ceil':k+=q!=k
        if mode=='nearest_even':k+=q-k>R(1,2)or q-k==R(1,2)and k%2==1
        return B.round(R(k,scale))
    phases=[[B.round(R(2*n+1,4)*pi)for n in (0,1)]]
    scales=[B.round(2/pi)]
    require(all(B.round(R(2*n+1,4)*(pi-pi_error))==B.round(R(2*n+1,4)*(pi+pi_error))==phases[0][n]for n in (0,1))
            and B.round(2/(pi-pi_error))==B.round(2/(pi+pi_error))==scales[0],'unresolved universal binary64 load')
    for digits in (8,9,10):
        for mode in ('floor','nearest_even','ceil'):
            phases.append([rounded_decimal(R(2*n+1,4)*pi,digits,mode)for n in (0,1)])
            scales.append(rounded_decimal(2/pi,digits,mode))
            require(all(rounded_decimal(R(2*n+1,4)*(pi-pi_error),digits,mode)==
                        rounded_decimal(R(2*n+1,4)*(pi+pi_error),digits,mode)==phases[-1][n]for n in (0,1))
                    and rounded_decimal(2/(pi-pi_error),digits,mode)==rounded_decimal(2/(pi+pi_error),digits,mode)==scales[-1],
                    'unresolved universal decimal quantizer')
    require(len(tables['phases'])==len(tables['scales'])==10,'universal table inventory differs')
    require(all([binary64(v['binary64'])for v in p['values']]==phases[i]for i,p in enumerate(tables['phases']))
            and all(binary64(p['value']['binary64'])==scales[i]for i,p in enumerate(tables['scales'])),'universal constants changed')
    scoped=[]
    for index,raw in enumerate(rows):
        if raw['function']!='BESSELJ'or raw['result_index']is not None or len(raw['args'])!=2 or not raw['expected'].startswith('0x'):continue
        try:x,n=map(binary64,raw['args'])
        except(ValueError,TypeError):continue
        if abs(x)>=two(24)and n in range(5):scoped.append((index,x,int(n),raw))
    require(len(scoped)==next(r for r in summary['counts']if r['function']=='BESSELJ')['observations'],'J replay scope differs')
    killed=set();native_executions=len(survivors)*len(scoped);orders={at:i for i,(at,*_)in enumerate(scoped)}
    path=directory/'BESSELJ-first-kills.tsv';frozen[str(path)]=digest(path)
    with path.open()as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        require(reader.fieldnames==['program','trace_index','got'],'asymptotic kill header differs')
        for kill in reader:
            identity,at=int(kill['program']),int(kill['trace_index'])
            require(0<=identity<67500 and identity not in killed and identity not in survivors and at in orders,
                    'duplicate or invalid asymptotic kill identity')
            require(kill['got']!=rows[at]['expected'],'invalid asymptotic contradiction')
            native_executions+=orders[at]+1;killed.add(identity)
    require(killed|set(survivors)==set(range(67500)) and
            native_executions==next(r for r in summary['counts']if r['function']=='BESSELJ')['executions'],
            'asymptotic search coverage or executions differ')
    def bits(q):return '0x'+struct.pack('>d',float(q)).hex()
    # Only seeds actually consumed by the declared recurrence need leaves.
    arguments=set()
    for p in programs:
        model=MODELS[p['body']]
        for _,x,n,_ in scoped:
            for seed in ((n,)if n<=1 else(0,1)):
                angle=operation(abs(x),-phases[p['phase']][seed],'add',model)
                arguments.add(bits(B.round(angle)))
    leaves={}
    for function in ('SIN','COS'):
        request=output.with_name(output.stem+'-'+function+'-requests.json')
        prediction=output.with_name(output.stem+'-'+function+'-predictions.jsonl')
        questions=[dict(id=f'r19-{function.lower()}-leaf-{i:03d}',args=[arg],input_transport='cell_value2_bulk')for i,arg in enumerate(sorted(arguments))]
        with request.open('x')as out:json.dump(dict(function=function,witnesses=questions),out);out.write('\n')
        subprocess.run([str(executable),'invoke-inputs',str(request),str(prediction)],check=True,timeout=300)
        actual=list(records(prediction));require(len(actual)==len(questions),'missing registered trig leaf')
        for q,r in zip(questions,actual):
            require(q['id']==r['id']and q['args']==r['args']and r['calc_func_prediction'].startswith('0x'),'trig leaf identity or domain differs')
            leaves[function,q['args'][0]]=binary64(r['calc_func_prediction'])
        frozen[str(request)]=digest(request);frozen[str(prediction)]=digest(prediction)
    @lru_cache(maxsize=4096)
    def sqrt_round(q,model):
        if model=='strict':q=B.round(q)
        fmt=B if model=='strict'else P53 if model.endswith('pc53')else E
        estimate=R(__import__('math').sqrt(float(q)));at=fmt.floor_coordinate(estimate)
        low,high=at-4096,at+4096
        require(fmt.value(low)**2<=q<=fmt.value(high)**2,'sqrt search does not enclose root')
        while low+1<high:
            middle=(low+high)//2
            if fmt.value(middle)**2<=q:low=middle
            else:high=middle
        a,b=fmt.value(low),fmt.value(high);midpoint=(a+b)/2
        result=a if q<midpoint*midpoint or q==midpoint*midpoint and low%2==0 else b
        return B.round(result)if model.startswith('stored')else result
    @lru_cache(maxsize=65536)
    def evaluate(x,n,phase_id,scale_id,body,recurrence,terms,association):
        model=MODELS[body];rec=MODELS[recurrence]
        def op(kind,a,b):return operation(a,b,kind,model)
        def seed(n):
            z=op('div',R(1),abs(x));z2=op('mul',z,z)
            angle=B.round(op('add',abs(x),-phases[phase_id][n]))
            sine,cosine=leaves['SIN',bits(angle)],leaves['COS',bits(angle)]
            a1,a2,a3=(R(-1,8),R(9,128),R(-75,1024))if n==0 else(R(3,8),R(-15,128),R(105,1024))
            p=op('add',R(1),-op('mul',a2,z2))if terms>=2 else R(1)
            q=R(0)if terms==0 else op('mul',op('add',a1,-op('mul',a3,z2))if terms>=2 else a1,z)
            expression=op('add',op('mul',cosine,p),-op('mul',sine,q))
            amplitude=sqrt_round(op('div',scales[scale_id],abs(x)),model)
            return op('mul',amplitude,expression)
        if n<=1:result=seed(n)
        else:
            previous,current=seed(0),seed(1);tox=operation(R(2),abs(x),'div',rec)
            for k in range(1,n):
                if association==0:step=operation(operation(R(k),tox,'mul',rec),current,'mul',rec)
                elif association==1:step=operation(operation(R(2*k),abs(x),'div',rec),current,'mul',rec)
                else:step=operation(operation(R(2*k),current,'mul',rec),abs(x),'div',rec)
                previous,current=current,operation(step,-previous,'add',rec)
            result=current
        return bits(B.round(-result if x<0 and n%2 else result))
    executions=0
    for p in programs:
        for _,x,n,row in scoped:
            got=evaluate(x,n,p['phase'],p['scale'],p['body'],p['recurrence'],p['terms'],p['association']);executions+=1
            require(got==row['expected'],'independent asymptotic survivor mismatch')
    require(frozen=={p:digest(p)for p in frozen},'asymptotic input/source changed')
    result=dict(complete=True,programs=len(programs),observations=len(scoped),executions=executions,registered_trig_arguments=len(arguments),input_sha256=frozen,
        limits='Independent exact arithmetic and rational sqrt-rounding replay of every declared J numeric survivor. '
               'SIN/COS leaves are current registered application calculations, not new Excel observations or identified private Bessel helpers. '
               'Domain obligations are excluded; no complete Excel kernel or product gain follows.')
    with output.open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def real_bound(raw):
    """Independent directed MPFR references, no producer arithmetic imports."""
    from reference import bessel_bounds
    return bessel_bounds(raw['function'],raw['args'][0],int(binary64(raw['args'][1])))


def check_regimes(input_path,directory,output):
    from reference import LIB
    start=time.monotonic();directory=Path(directory)
    rows=list(records(input_path));proofs=records(directory/'reference-bounds.jsonl')
    summary=json.loads((directory/'summary.json').read_text())
    frozen={str(p):digest(p)for p in (Path(input_path),directory/'reference-bounds.jsonl',directory/'summary.json',Path(__file__),
                                    Path(__file__).with_name('reference.py'),Path(__file__).with_name('exact.py'))}
    require(summary['input_sha256']==digest(input_path)and summary['bounds_sha256']==digest(directory/'reference-bounds.jsonl'),'reference hashes differ')
    seen=set();counts=Counter();cache={};regions=defaultdict(Counter)
    for proof in proofs:
        index=proof['trace_index'];require(index not in seen and 0<=index<len(rows),'duplicate reference identity');seen.add(index)
        raw=rows[index]
        require(all(raw[k]==proof[k]for k in ('id','function','args','expected')),'reference identity changed')
        key=raw['function'],tuple(raw['args'])
        if key not in cache:cache[key]=real_bound(raw)
        low,high=cache[key];declared=list(map(decode_pair,proof['proof']['enclosure']))
        require(declared[0]<=low<=high<=declared[1],'independent real bound escapes declared enclosure')
        target=binary64(raw['expected'])
        impossible=not B.round(declared[0])<=target<=B.round(declared[1])
        require(proof['outside_real_rounding_bound']==impossible,'incorrect exact-reference exclusion')
        x,n=map(binary64,raw['args'])
        region=('seed'if n<=1 else'above_order'if abs(x)>n else'at_or_below_order')+('_small'if abs(x)<=32 else'_large')
        require(proof['region']==region,'reference region differs')
        regions[raw['function']+':'+region]['observations']+=1
        regions[raw['function']+':'+region]['outside_real_rounding_bound']+=impossible
        counts['observations']+=1;counts['excluded']+=impossible
    require(counts['observations']==summary['counts']['scoped']and counts['excluded']==summary['counts']['outside_real_rounding_bound'],'reference counts differ')
    require(dict(regions)==summary['regions'],'reference region counts differ')
    eligible=set()
    for index,raw in enumerate(rows):
        if raw['function']not in ('BESSELJ','BESSELY')or raw['result_index']is not None or len(raw['args'])!=2 or not raw['expected'].startswith('0x'):continue
        try:x,n=map(binary64,raw['args'])
        except (ValueError,TypeError):continue
        if (n.denominator==1 and 0<=n<=64 and (raw['function']=='BESSELJ'or x>0)
                and (abs(x)>=R(1,1024)or x==0 and raw['function']=='BESSELJ')):eligible.add(index)
    require(seen==eligible,'reference scope omits or adds eligible observations')
    require(frozen=={p:digest(p)for p in frozen},'reference inputs changed')
    result=dict(complete=True,counts=counts,regions=regions,seconds=time.monotonic()-start,input_sha256=frozen,mpfr_version=LIB.mpfr_get_version().decode(),
        method='Independent directed256-bit MPFR real Jn/Yn bounds through the official shared-library API; no library source read.',
        limits='Certifies the supplied real-function exclusions; does not identify Excel series, seed approximation, recurrence direction or asymptotic program. No product answer is emitted.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='asymptotic-survivors':check_asymptotic_survivors(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='regimes':check_regimes(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='miller-sites':check_miller(*sys.argv[2:], expanded=True, sites=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='miller-forms':check_miller(*sys.argv[2:], expanded=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='miller':check_miller(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='forward':check_forward(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='domains':check_domains(*sys.argv[2:])
    else:check(*sys.argv[1:])
