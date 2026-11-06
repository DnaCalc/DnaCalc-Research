"""Isolated PMT timing/numerator stages constrained by public type-zero calls.

Original MIT tooling; no logarithm or exponential producer is assumed. Shared
binary64 em values are exact necessary constraints, never fitted helpers.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path
import random
import subprocess

from fingerprint import CachedService, load_service, bits, value, sha, write
from exact import BINARY64 as B, Span, binary64, inverse, operation, round_model

SOURCE_BINDINGS={str(p.resolve()):sha(p) for p in (
    Path(__file__),Path(__file__).with_name('x87_programs.py'),
    Path(__file__).with_name('x87_native.c'),Path(__file__).parents[1]/'checker/exact.py')}


def run(checkout, output, cache, offline=False, negative=False):
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    # Fixed ordinary rates, periods and scaled/adjacent balances; no bank
    # failures or current application predictions select these questions.
    rates = (.0025, .0075, .03, .1, .125, .25, .5, 2.**-10)
    periods = (1., 12., 360., 1024.)
    if negative:
        rates=tuple(-r for r in rates)
        periods=(1.,12.,60.,120.)
    final_model='strict' if negative else 'stored_pc64'
    balances = [.125, 1., 1.1, 1000., 125000., 2.**20]
    balances += [value(f'{int(bits(1.),16)+k:016x}') for k in (-2,-1,1,2)]
    rows = []
    for r in rates:
        for n in periods:
            for pv in balances:
                for timing in (0., 1.):
                    rows.append([r,n,pv,0.,timing])
            for pv,fv in ((0.,1.), (0.,1.1), (1.,1.), (1.,.7), (1.,-1.), (1000.,-700.)):
                rows.append([r,n,pv,fv,0.])
    questions = [[client.number(x) for x in row] for row in rows]
    answers = service.ask('PMT', questions)['rows']
    observed = {tuple(row): answer for row, answer in zip(rows, answers)}
    totals, groups = Counter(), []
    factor_groups = {}
    for rate in rates:
        for n in periods:
            r = binary64(bits(rate))
            domain = None
            contradicted = False
            constraints = []
            for pv in balances:
                answer = observed[(rate,n,pv,0.,0.)]
                if answer['type'] != 'number':
                    raise ValueError('positive finite loan must have numeric payment')
                y = binary64(answer['bits'])
                q = inverse(Span(y,y), r, 'mul', True, final_model)
                em = inverse(q, binary64(bits(pv)), 'div', False, 'stored_pc64')
                constraints.append(em.json() if em else None)
                if not contradicted:
                    domain = em if domain is None else domain.meet(em)
                    if domain is None:
                        contradicted = True
            coordinates = B.coordinates(domain) if domain else None
            em_values = [] if coordinates is None else [B.value(k) for k in
                range(coordinates[0], min(coordinates[1]+1, coordinates[0]+8))]
            total_values = 0 if coordinates is None else coordinates[1]-coordinates[0]+1
            if total_values > 8:
                raise ValueError('unresolved broad helper interval; extend balances before comparing stages')
            totals['groups'] += 1
            totals['unique_em_groups'] += total_values == 1
            totals['empty_em_groups'] += total_values == 0
            group = dict(rate=bits(rate), periods=bits(n), type_zero_constraints=constraints,
                         em_domain=domain.json() if domain else None,
                         em_values=[bits(float(e)) for e in em_values], timing=[], numerator=[])
            tf = round_model(1+r, 'strict')
            for pv in balances:
                p = binary64(bits(pv))
                answer = observed[(rate,n,pv,0.,1.)]
                zero = observed[(rate,n,pv,0.,0.)]
                candidates = {name:[] for name in ('timing_before_rate','outer_timing','due_coefficient')}
                candidates['outer_timing'] = [operation(binary64(zero['bits']),tf,'div','stored_pc64')]
                for em in em_values:
                    q = operation(p,em,'div','stored_pc64')
                    candidates['timing_before_rate'].append(operation(operation(q,tf,'div','stored_pc64'),r,'mul',final_model))
                    candidates['due_coefficient'].append(operation(q,operation(r,tf,'div','stored_pc64'),'mul',final_model))
                target = binary64(answer['bits'])
                if len(em_values) == 1:
                    q = operation(p,em_values[0],'div','stored_pc64')
                    factor_domain = inverse(Span(target,target),q,'mul',True,final_model)
                    fg = factor_groups.setdefault(bits(rate),dict(domain=None,contradicted=False,constraints=0))
                    fg['constraints'] += 1
                    if not fg['contradicted']:
                        fg['domain'] = factor_domain if fg['domain'] is None else fg['domain'].meet(factor_domain)
                        if fg['domain'] is None:
                            fg['contradicted'] = True
                            fg['first_contradiction'] = dict(periods=bits(n),pv=bits(pv))
                possible = {name:target in predictions for name,predictions in candidates.items()}
                group['timing'].append(dict(pv=bits(pv), zero=zero, due=answer,
                    predictions={k:[bits(float(v)) for v in vv] for k,vv in candidates.items()}, possible=possible))
                totals['timing_rows'] += 1
                for name, v in possible.items():
                    totals[name+'_possible'] += v
            for pv,fv in ((0.,1.), (0.,1.1), (1.,1.), (1.,.7), (1.,-1.), (1000.,-700.)):
                answer = observed[(rate,n,pv,fv,0.)]
                target = binary64(answer['bits'])
                predictions = {}
                for product_model in ('strict','stored_pc64'):
                    for sum_model in ('strict','stored_pc64'):
                        label = product_model+'/'+sum_model
                        values = []
                        for em in em_values:
                            inv_factor = round_model(1+em,'strict')
                            num = operation(binary64(bits(pv)), operation(binary64(bits(fv)),inv_factor,
                                'mul',product_model), 'add', sum_model)
                            values.append(operation(operation(num,em,'div','stored_pc64'),r,'mul',final_model))
                        predictions[label] = values
                for add_model in ('strict','stored_pc64'):
                    values = []
                    for em in em_values:
                        p, f = binary64(bits(pv)), binary64(bits(fv))
                        balance = operation(p,f,'add',add_model)
                        q = operation(balance,em,'div','stored_pc64')
                        corrected = operation(q,f,'add',add_model)
                        values.append(operation(corrected,r,'mul',final_model))
                    predictions['balance_quotient/'+add_model] = values
                possible = {name:target in values for name,values in predictions.items()}
                group['numerator'].append(dict(pv=bits(pv),fv=bits(fv),answer=answer,possible=possible,
                    predictions={k:[bits(float(v)) for v in vv] for k,vv in predictions.items()}))
                totals['numerator_rows'] += 1
                for name,v in possible.items():
                    totals['numerator_'+name+'_possible'] += v
            groups.append(group)
    for rate, group in factor_groups.items():
        domain = group.pop('domain')
        coordinates = B.coordinates(domain) if domain else None
        r = binary64(rate)
        tf = round_model(1+r,'stored_pc64')
        ratio = operation(r,tf,'div','stored_pc64')
        candidates = dict(due_ratio=ratio,
            reciprocal=operation(1,operation(operation(1,r,'div','stored_pc64'),1,'add','stored_pc64'),'div','stored_pc64'),
            subtract_product=operation(r,-operation(r,ratio,'mul','stored_pc64'),'add','stored_pc64'),
            complement_product=operation(r,operation(1,-ratio,'add','stored_pc64'),'mul','stored_pc64'))
        group.update(domain=domain.json() if domain else None,
            binary64_count=0 if coordinates is None else coordinates[1]-coordinates[0]+1,
            candidates={k:bits(float(v)) for k,v in candidates.items()},
            possible={k:domain.contains(v) if domain else False for k,v in candidates.items()})
    result = dict(format='dnacalc-annuity-paired-stages/1', complete=True, counts=totals,
        branch='negative_rate_growth' if negative else 'positive_rate_discount', final_model=final_model,
        groups=groups, factor_groups=factor_groups, receipts=service.receipts, source_bindings=SOURCE_BINDINGS,
        limits='Positive-rate finite PMT only. One shared binary64 em per rate/period is constrained by ten independent type-zero balances. Timing and nonzero-future-value stages are tested separately against all remaining em values. Nonempty constraints do not identify the logarithm or expm1; no parameter or selector is fitted. Opposing balances retain their separately measured cancellation route.')
    write(output,result)
    print(json.dumps(totals,indent=2),flush=True)
    return result


def discriminants(checkout, prior, output, cache, offline=False):
    """Select isolated sum ties and reciprocal differences before asking Excel."""
    data=json.loads(Path(prior).read_text())
    if not data.get('complete'):
        raise ValueError('incomplete prior helper constraints')
    final_model=data.get('final_model','stored_pc64')
    negative=data.get('branch')=='negative_rate_growth'
    questions=[]
    def predict(pv,fv,rate,em,first,last):
        s=operation(pv,fv,'add',first)
        q=operation(s,em,'div','stored_pc64')
        v=operation(q,fv,'add',last)
        return operation(v,rate,'mul',final_model)
    models=('strict','stored_pc64')
    for g in data['groups']:
        if len(g['em_values']) != 1:
            continue
        r,n,em=map(binary64,(g['rate'],g['periods'],g['em_values'][0]))
        for stage in ('first_sum','last_sum'):
            for exponent in (-10,0,10):
                scale=2.**exponent
                for sign in (-1.,1.):
                    if stage=='first_sum':
                        pv=scale
                        fv=sign*scale*(2.**-53+2.**-66)
                    else:
                        q=operation(binary64(bits(scale)),em,'div','stored_pc64')
                        fv=sign*(math.ulp(float(q))*.5+math.ulp(float(q))*2.**-14)
                        pv=scale-fv
                    pp,ff=map(binary64,(bits(pv),bits(fv)))
                    predictions={a+'/'+b:bits(float(predict(pp,ff,r,em,a,b))) for a in models for b in models}
                    # Only the named sum distinguishes the two predictions.
                    if stage=='first_sum':
                        selected=(predictions['strict/strict']==predictions['strict/stored_pc64']
                            and predictions['stored_pc64/strict']==predictions['stored_pc64/stored_pc64']
                            and predictions['strict/strict']!=predictions['stored_pc64/strict'])
                    else:
                        selected=(predictions['strict/strict']==predictions['stored_pc64/strict']
                            and predictions['strict/stored_pc64']==predictions['stored_pc64/stored_pc64']
                            and predictions['strict/strict']!=predictions['strict/stored_pc64'])
                    if selected:
                        questions.append(dict(stage=stage,args=[g['rate'],g['periods'],bits(pv),bits(fv),bits(0.)],
                            em=g['em_values'][0],predictions=predictions))
    if not all(any(q['stage']==s for q in questions) for s in ('first_sum','last_sum')):
        raise ValueError('no isolated discriminant for one of the sums')
    # No Excel outputs participate in selecting new rate discriminants.
    rng=random.Random(29)
    rate_questions=[]
    wanted={'strict_reciprocal':8,'resident_ratio':8,'resident_reciprocal':8}
    selected=Counter()
    for trial in range(100000):
        rh=f'{((1023-rng.randint(2,20))<<52)|rng.getrandbits(52):016x}'
        if negative:
            rh=f'{int(rh,16)|(1<<63):016x}'
        r=binary64(rh)
        def recip(model):
            a=operation(1,r,'div',model)
            b=operation(1,a,'add',model)
            return operation(1,b,'div',model)
        compiled=recip('stored_pc64')
        variants=dict(stored_reciprocal=compiled,strict_reciprocal=recip('strict'),
            resident_ratio=round_model(r/round_model(1+r,'resident_pc64'),'stored_pc64'),
            resident_reciprocal=round_model(1/round_model(1+round_model(1/r,'resident_pc64'),'resident_pc64'),'stored_pc64'),
            strict_inner=round_model(1/round_model(1+round_model(1/r,'strict'),'stored_pc64'),'stored_pc64'),
            strict_return=round_model(1/round_model(1+round_model(1/r,'stored_pc64'),'stored_pc64'),'strict'),
            resident_return=round_model(1/round_model(1+round_model(1/r,'stored_pc64'),'stored_pc64'),'resident_pc64'))
        classes=[name for name in wanted if variants[name]!=compiled and selected[name]<wanted[name]]
        if not classes:
            continue
        for name in classes:
            selected[name]+=1
        rate_questions.append(dict(rate=rh,classes=classes,
            factors={k:bits(float(v)) for k,v in variants.items()},
            exact_factors={k:[str(v.numerator),str(v.denominator)] for k,v in variants.items()}))
        if all(selected[name]>=wanted[name] for name in wanted):
            break
    if not all(selected[name]>=wanted[name] for name in wanted):
        raise ValueError('reciprocal discriminator search incomplete')
    rows=[q['args'] for q in questions]
    balances=[.125,1.,1.1,1000.,125000.,2.**20]
    balances += [value(f'{int(bits(1.),16)+k:016x}') for k in (-2,-1,1,2)]
    for rq in rate_questions:
        for pv in balances:
            for timing in (0.,1.):
                rows.append([rq['rate'],bits(12.),bits(pv),bits(0.),bits(timing)])
    # Persist predictions and selection before submitting any fresh call.
    manifest=dict(format='dnacalc-annuity-stage-discriminants-questions/1',complete=True,
        prior_sha256=sha(prior),sums=questions,rates=rate_questions,rows=rows,
        selection_counts=selected,source_sha256=sha(__file__))
    write(str(output)+'.questions.json',manifest)
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('PMT',[[client.number(x) for x in row] for row in rows])['rows']
    totals=Counter()
    for i,q in enumerate(questions):
        answer=answers[i]
        q['answer']=answer
        totals[q['stage']+'_rows']+=1
        for name,prediction in q['predictions'].items():
            totals[q['stage']+'_'+name]+=answer==dict(type='number',bits=prediction)
    pos=len(questions)
    for rq in rate_questions:
        r=binary64(rq['rate']);em_domain=None;factor_domain=None;targets=[]
        factor_failed=False
        for pv in balances:
            zero,due=answers[pos:pos+2];pos+=2
            if zero['type']!=due['type'] or zero['type']!='number':
                raise ValueError('positive finite rate stage calls must be numeric')
            p=binary64(bits(pv)); y=binary64(zero['bits'])
            domain=inverse(inverse(Span(y,y),r,'mul',True,final_model),p,'div',False,'stored_pc64')
            em_domain=domain if em_domain is None else em_domain.meet(domain)
            if em_domain is None:
                raise ValueError('fresh type-zero calls contradict the shared denominator stage')
            targets.append((pv,zero,due))
        coords=B.coordinates(em_domain) if em_domain else None
        if coords is None or coords[0]!=coords[1]:
            raise ValueError('new rate did not have one constrained denominator')
        em=B.value(coords[0])
        rq.update(em=bits(float(em)),em_domain=em_domain.json(),observations=[])
        for pv,zero,due in targets:
            q=operation(binary64(bits(pv)),em,'div','stored_pc64')
            target=binary64(due['bits'])
            domain=inverse(Span(target,target),q,'mul',True,final_model)
            if not factor_failed:
                factor_domain=domain if factor_domain is None else factor_domain.meet(domain)
                factor_failed=factor_domain is None
            rq['observations'].append(dict(pv=bits(pv),zero=zero,due=due))
        rq['factor_domain']=factor_domain.json() if factor_domain else None
        from fractions import Fraction
        rq['possible']={k:factor_domain.contains(Fraction(int(v[0]),int(v[1]))) if factor_domain else False
                       for k,v in rq['exact_factors'].items()}
        totals['new_rates']+=1
        for name,possible in rq['possible'].items():
            totals[name+'_possible']+=possible
    if pos!=len(answers):
        raise ValueError('unaccounted answer')
    result=dict(format='dnacalc-annuity-stage-discriminants/1',complete=True,counts=totals,
        branch=data.get('branch','positive_rate_discount'), final_model=final_model,
        sums=questions,rates=rate_questions,receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        source_bindings=SOURCE_BINDINGS,
        limits='Questions selected only by exact arithmetic predictions and earlier unique denominator constraints. Each balance tie isolates one addition; rate questions distinguish fixed reciprocal programs before observing Excel. Denominators for fresh rates are independently constrained from public type-zero calls. No fitted helper or coefficient is installed.')
    write(output,result)
    print(json.dumps(totals,indent=2),flush=True)
    return result


def compose(checkout,prior,output,cache,offline=False):
    """Fresh coupled calls test composition of the independently pinned stages."""
    data=json.loads(Path(prior).read_text())
    require_complete=data.get('complete') is True
    if not require_complete:
        raise ValueError('incomplete paired stage study')
    rows=[]
    for group in data['groups']:
        if len(group['em_values'])!=1:
            raise ValueError('composed stage needs one independently constrained denominator')
        r=binary64(group['rate']);em=binary64(group['em_values'][0])
        coefficient=operation(1,operation(1,operation(1,r,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64')
        for original in group['numerator']:
            pv,fv=map(binary64,(original['pv'],original['fv']))
            q=operation(operation(pv,fv,'add','stored_pc64'),em,'div','stored_pc64')
            corrected=operation(q,fv,'add','stored_pc64')
            prediction=operation(corrected,coefficient,'mul',data.get('final_model','stored_pc64'))
            rows.append(dict(args=[group['rate'],group['periods'],original['pv'],original['fv'],'3ff0000000000000'],
                em=group['em_values'][0],prediction=bits(float(prediction)),opposing=pv!=0 and fv==-pv))
    manifest=dict(format='dnacalc-annuity-composed-questions/1',complete=True,rows=rows,
        prior_sha256=sha(prior),branch=data.get('branch'),final_model=data.get('final_model','stored_pc64'),source_sha256=sha(__file__))
    write(str(output)+'.questions.json',manifest)
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('PMT',[[client.number(x) for x in row['args']] for row in rows])['rows']
    counts=Counter()
    for row,answer in zip(rows,answers):
        row['answer']=answer
        row['matched']=answer==dict(type='number',bits=row['prediction'])
        counts['rows']+=1;counts['matched']+=row['matched']
        counts[('opposing' if row['opposing'] else 'ordinary')+'_rows']+=1
        counts[('opposing' if row['opposing'] else 'ordinary')+'_matched']+=row['matched']
    result=dict(format='dnacalc-annuity-composed-stages/1',complete=True,counts=counts,records=rows,
        prior_sha256=sha(prior),branch=data.get('branch'),final_model=data.get('final_model','stored_pc64'),
        receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        source_bindings=SOURCE_BINDINGS,
        limits='New nonzero-future-value beginning-payment calls, predicted before observing Excel. The denominator comes only from earlier independent type-zero calls. Opposing balances are identified mathematically and counted separately because the product already has a separately measured cancellation route. No answer-selected coefficient or selector is introduced.')
    write(output,result)
    print(json.dumps(counts),flush=True)
    return result


def edges(checkout,output,cache,offline=False):
    """Adjacent rate-domain and balance-overflow controls for the same fold."""
    maximum=value('7fefffffffffffff')
    rows=[]
    for k in range(1,17):
        r=value(f'{0xbff0000000000000-k:016x}')
        for n in (1.,2.,12.):
            for fv in (0.,1.):
                for timing in (0.,1.):
                    rows.append([r,n,1.,fv,timing])
    for r in (.1,.5,1e308):
        for n in (1.,12.,512.):
            for pv,fv in ((.75*maximum,.75*maximum),(.75*maximum,-.75*maximum),
                          (maximum,0.),(maximum,-maximum),(0.,maximum)):
                for timing in (0.,1.):
                    rows.append([r,n,pv,fv,timing])
    questions=[[dict(type='number',bits=bits(x)) for x in row] for row in rows]
    write(str(output)+'.questions.json',dict(format='dnacalc-annuity-exponent-edge-questions/1',complete=True,
        rows=questions,source_sha256=sha(__file__),selection='Sixteen adjacent doubles above rate -1; three fixed positive rates, three periods, five fixed large-balance patterns, both timing values.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answer=service.ask('PMT',questions)['rows']
    result=dict(format='dnacalc-annuity-exponent-edges/1',complete=True,
        records=[dict(args=[a['bits'] for a in q],answer=a) for q,a in zip(questions,answer)],
        receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        source_bindings=SOURCE_BINDINGS,
        limits='Original numeric/error readings at domain and overflow edges; no predicted value is treated as an answer. These controls measure the changed finite fold and reciprocal coefficient outside the ordinary-rate stage probes.')
    write(output,result);print(Counter(a['type'] for a in answer),flush=True)
    return result


def products(checkout,prior,output,cache,offline=False):
    """Mine final-product ties with a denominator pinned by earlier public calls."""
    data=json.loads(Path(prior).read_text())
    if data.get('complete') is not True:
        raise ValueError('incomplete denominator study')
    rng=random.Random(2901);rows=[];excluded=[]
    for group in data['groups']:
        if len(group['em_values'])!=1:
            raise ValueError('unresolved denominator')
        rate=binary64(group['rate']);em=binary64(group['em_values'][0])
        magnitude=abs(rate)
        if magnitude.numerator.bit_count()==1 and magnitude.denominator.bit_count()==1:
            excluded.append(dict(rate=group['rate'],periods=group['periods'],reason='Power-of-two multiplication is exact on these normal operands.'))
            continue
        selected=0
        for trial in range(100000):
            pv=f'{((1023+rng.randint(-10,10))<<52)|rng.getrandbits(52):016x}'
            q=operation(binary64(pv),em,'div','stored_pc64')
            single=operation(q,rate,'mul','strict')
            double=operation(q,rate,'mul','stored_pc64')
            if single==double:
                continue
            rows.append(dict(args=[group['rate'],group['periods'],pv,'0000000000000000','0000000000000000'],
                em=group['em_values'][0],strict=bits(float(single)),stored_pc64=bits(float(double))))
            selected+=1
            if selected==2:
                break
        if selected!=2:
            raise ValueError('no complete final-product discriminant set')
    write(str(output)+'.questions.json',dict(format='dnacalc-annuity-product-questions/1',complete=True,
        rows=rows,excluded=excluded,prior_sha256=sha(prior),source_bindings=SOURCE_BINDINGS,
        selection='First two exact strict/PC64 product discriminants per non-power-of-two rate/period group; seeded binary64 balance draws, before Excel answers.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('PMT',[[client.number(x) for x in row['args']] for row in rows])['rows']
    counts=Counter()
    for row,answer in zip(rows,answers):
        row['answer']=answer
        row['classification']=('strict' if answer==dict(type='number',bits=row['strict'])
            else 'stored_pc64' if answer==dict(type='number',bits=row['stored_pc64']) else 'neither')
        counts[row['classification']]+=1
    result=dict(format='dnacalc-annuity-product-discriminants/1',complete=True,counts=counts,
        records=rows,excluded=excluded,prior_sha256=sha(prior),receipts=service.receipts,
        questions_sha256=sha(str(output)+'.questions.json'),source_bindings=SOURCE_BINDINGS,
        limits='One final multiplication, type0/fv0, with independently pinned binary64 em and a stored-PC64 quotient. Neither results refute the conditional earlier-stage model; no result changes a denominator or selector. Power-of-two factors cannot distinguish the two final roundings here.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def timing_edges(checkout,output,cache,offline=False):
    """Fixed exponent controls for a reciprocal timing coefficient."""
    rates=[2.**e for e in (-1022,-1000,-64,-53,-52,-10,0,52,53,500,1023)]
    rates += [1e308,value('7fefffffffffffff')]
    periods=(.5,1.,2.,12.,360.,-1.,-12.)
    balances=((.125,0.),(1.,0.),(1.1,0.),(1.,.7),(1.,-1.),(0.,1.),(1000.,-700.))
    rows=[[r,n,p,f,t] for r in rates for n in periods for p,f in balances for t in (0.,1.)]
    questions=[[dict(type='number',bits=bits(v)) for v in row] for row in rows]
    write(str(output)+'.questions.json',dict(format='dnacalc-annuity-timing-exponent-questions/1',
        complete=True,rows=questions,source_bindings=SOURCE_BINDINGS,
        selection='Fixed normal-rate exponent edges, maximum finite rates, ordinary balances and future values, both timings, and negative-period unchanged-route controls.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('PMT',questions)['rows']
    result=dict(format='dnacalc-annuity-timing-exponent-edges/1',complete=True,
        records=[dict(args=[a['bits'] for a in q],answer=a) for q,a in zip(questions,answers)],
        receipts=service.receipts,source_bindings=SOURCE_BINDINGS,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Original numeric/error readings at preselected exponent and mathematical branch controls. '
        'No edge is an inferred approximation cutoff and no answer is a model prediction. '
        'Negative-period controls do not establish a new timing model for that route.')
    write(output,result);print(json.dumps(Counter(a['type'] for a in answers)),flush=True);return result


def cumulative_pairs(checkout,output,cache,offline=False):
    """Separate the cumulative singleton route from its summation stage."""
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    groups=[];points=[];spans=[]
    for rate in (.0025,.03,.125,2.**-10):
        for n in (12,60,360):
            for pv in (1.,1.1,1000.):
                for timing in (0.,1.):
                    group=dict(rate=bits(rate),periods=bits(n),pv=bits(pv),timing=bits(timing))
                    groups.append(group)
                    intervals=sorted({(start,min(start+length-1,n)) for start in (1,2,5,n//2,n-1)
                                      for length in (2,3)})
                    ks=sorted({k for start,end in intervals for k in range(start,end+1)}|{n})
                    points.extend(dict(group,period=bits(k)) for k in ks)
                    spans.extend(dict(group,start=start,end=end) for start,end in intervals)
    def pmt_args(g):return [g['rate'],g['periods'],g['pv'],bits(0.),g['timing']]
    def part_args(p):return [p['rate'],p['period'],p['periods'],p['pv'],bits(0.),p['timing']]
    def cumulative_args(r,start,end):return [r['rate'],r['periods'],r['pv'],bits(start),bits(end),r['timing']]
    def ask(name,rows):return service.ask(name,[[client.number(v) for v in row] for row in rows])['rows']
    key=lambda g:tuple(g[k] for k in ('rate','periods','pv','timing'))
    write(str(output)+'.questions.json',dict(format='dnacalc-cumulative-public-route-questions/1',complete=True,
        groups=groups,points=points,spans=spans,source_bindings=SOURCE_BINDINGS,
        selection='Fixed rates, periods, balances and timings. Short intervals expose singleton publications and ascending sums separately.'))
    payments=dict(zip(map(key,groups),ask('PMT',[pmt_args(g) for g in groups])))
    parts={name:ask(name,[part_args(p) for p in points]) for name in ('PPMT','IPMT')}
    point_predictions=[]
    for i,p in enumerate(points):
        payment,principal,interest=payments[key(p)],parts['PPMT'][i],parts['IPMT'][i]
        if any(a['type']!='number' for a in (payment,principal,interest)):
            raise ValueError('ordinary financial controls must be numeric')
        m,pr,it=map(lambda a:binary64(a['bits']),(payment,principal,interest))
        point_predictions.append(dict(p,payment=payment,public_ppmt=principal,public_ipmt=interest,
            principal_predictions=dict(public_ppmt=principal['bits'],
                payment_minus_interest_strict=bits(float(operation(m,-it,'add','strict'))),
                payment_minus_interest_stored_pc64=bits(float(operation(m,-it,'add','stored_pc64'))))))
    write(str(output)+'.singleton-predictions.json',dict(format='dnacalc-cumulative-singleton-predictions/1',
        complete=True,records=point_predictions,receipts=service.receipts,source_bindings=SOURCE_BINDINGS))
    principals=ask('CUMPRINC',[cumulative_args(p,value(p['period']),value(p['period'])) for p in points])
    interests=ask('CUMIPMT',[cumulative_args(p,value(p['period']),value(p['period'])) for p in points])
    lookup={};counts=Counter()
    for i,p in enumerate(point_predictions):
        p['cumprinc']=principals[i];p['cumipmt']=interests[i]
        if principals[i]['type']!='number' or interests[i]['type']!='number':
            raise ValueError('ordinary cumulative singleton must be numeric')
        lookup[key(p)+(int(value(p['period'])),)]=p
        counts['singleton_rows']+=1
        for model,want in p['principal_predictions'].items():
            counts['singleton_principal/'+model]+=principals[i]==dict(type='number',bits=want)
        counts['singleton_interest/public_ipmt']+=interests[i]==p['public_ipmt']
        for model in ('strict','stored_pc64'):
            predicted=bits(float(operation(binary64(p['payment']['bits']),-binary64(principals[i]['bits']),'add',model)))
            counts['singleton_interest/payment_minus_cumprinc_'+model]+=interests[i]==dict(type='number',bits=predicted)
    span_predictions=[]
    for r in spans:
        point_rows=[lookup[key(r)+(k,)] for k in range(r['start'],r['end']+1)]
        predictions={}
        for source in ('cumprinc','public_ppmt'):
            values=[binary64(p[source]['bits']) for p in point_rows]
            for model in ('strict','stored_pc64','resident_pc64'):
                total=0
                for v in values:total=operation(total,v,'add',model)
                predictions[source+'/'+model]=bits(float(round_model(total,'strict')))
        span_predictions.append(dict(r,predictions=predictions))
    write(str(output)+'.span-predictions.json',dict(format='dnacalc-cumulative-span-predictions/1',complete=True,
        records=span_predictions,singletons=point_predictions,receipts=service.receipts,source_bindings=SOURCE_BINDINGS))
    span_answers=ask('CUMPRINC',[cumulative_args(r,r['start'],r['end']) for r in spans])
    for row,answer in zip(span_predictions,span_answers):
        row['answer']=answer;counts['span_rows']+=1
        for model,prediction in row['predictions'].items():
            counts['span/'+model]+=answer==dict(type='number',bits=prediction)
    result=dict(format='dnacalc-cumulative-public-route-pairs/1',complete=True,counts=counts,
        singletons=point_predictions,spans=span_predictions,receipts=service.receipts,source_bindings=SOURCE_BINDINGS,
        questions_sha256=sha(str(output)+'.questions.json'),
        singleton_predictions_sha256=sha(str(output)+'.singleton-predictions.json'),
        span_predictions_sha256=sha(str(output)+'.span-predictions.json'),
        limits='Original PMT/IPMT/PPMT leaves precede independent cumulative singleton answers, '
        'which precede the short-interval answers. Fixed scalar routes and three ascending-sum '
        'precisions are tested separately. A contradiction excludes only that exact public '
        'composition, not an unpublished component or arbitrary cumulative algorithm. No '
        'principal, summation selector, tolerance or coefficient is fitted or installed.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def cumulative_fv(checkout,study,output,cache,offline=False):
    """Condition end-payment principal expressions on independent public FV."""
    from geomean import verified_study_answers
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete cumulative pairs')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    frozen={**SOURCE_BINDINGS,**bindings,str(source.resolve()):sha(source)}
    points=[p for p in doc['singletons'] if p['timing']==bits(0.)]
    spans=[r for r in doc['spans'] if r['timing']==bits(0.)]
    questions=[[p['rate'],bits(value(p['period'])-1.),p['payment']['bits'],p['pv'],bits(0.)] for p in points]
    write(str(output)+'.questions.json',dict(format='dnacalc-cumulative-fv-leaf-questions/1',complete=True,
        rows=questions,source_bindings=frozen,
        hypothesis='End-payment principal is PMT-FV(rate,k-1,PMT,pv,0)*rate, accumulated in ascending period order.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('FV',[[client.number(x) for x in row] for row in questions])['rows']
    after_questions=[[p['rate'],p['period'],p['payment']['bits'],p['pv'],bits(0.)] for p in points]
    write(str(output)+'.balance-questions.json',dict(format='dnacalc-cumulative-fv-balance-questions/1',complete=True,
        rows=after_questions,source_bindings=frozen,
        hypothesis='Principal over an interval is public FV before the first period minus public FV after the last period.'))
    after_answers=service.ask('FV',[[client.number(x) for x in row] for row in after_questions])['rows']
    # Fixed simple source forms, ranked by code shape, not observed matches.
    models=dict(compiled_expression=('resident_pc64','resident_pc64','stored_pc64'),
                stored_assignments=('stored_pc64','stored_pc64','stored_pc64'),
                strict=('strict','strict','strict'),
                stored_interest=('stored_pc64','resident_pc64','stored_pc64'))
    counts=Counter();records=[];lookup={}
    def key(p):return tuple(p[k] for k in ('rate','periods','pv','timing'))
    def term(p,fv,model):
        multiply,subtract,_=models[model]
        interest=operation(binary64(fv['bits']),binary64(p['rate']),'mul',multiply)
        return operation(binary64(p['payment']['bits']),-interest,'add',subtract)
    for p,fv,after_fv in zip(points,answers,after_answers):
        if fv['type']!='number' or after_fv['type']!='number':raise ValueError('finite FV control must be numeric')
        lookup[key(p)+(int(value(p['period'])),)]=(p,fv,after_fv)
        predictions={name:bits(float(round_model(term(p,fv,name),'stored_pc64'))) for name in models}
        for model in ('strict','stored_pc64'):
            predictions['balance_'+model]=bits(float(operation(binary64(fv['bits']),-binary64(after_fv['bits']),'add',model)))
        records.append(dict(p,public_fv=fv,public_fv_after=after_fv,predictions=predictions));counts['singletons']+=1
        for name,b in predictions.items():counts['singleton/'+name]+=p['cumprinc']==dict(type='number',bits=b)
    interval_records=[]
    for row in spans:
        predictions={}
        for name,(_,_,summation) in models.items():
            total=0
            for k in range(row['start'],row['end']+1):
                p,fv,_=lookup[key(row)+(k,)]
                total=operation(total,term(p,fv,name),'add',summation)
            predictions[name]=bits(float(round_model(total,'strict')))
        _,first_fv,_=lookup[key(row)+(row['start'],)]
        _,_,last_fv=lookup[key(row)+(row['end'],)]
        for model in ('strict','stored_pc64'):
            predictions['balance_'+model]=bits(float(operation(binary64(first_fv['bits']),-binary64(last_fv['bits']),'add',model)))
        interval_records.append(dict(row,public_fv_predictions=predictions));counts['spans']+=1
        for name,b in predictions.items():counts['span/'+name]+=row['answer']==dict(type='number',bits=b)
    result=dict(format='dnacalc-cumulative-public-fv-route/1',complete=True,counts=counts,models=models,
        singletons=records,spans=interval_records,receipts=doc['receipts']+service.receipts,
        source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Discovery comparison on the earlier end-payment cumulative cohort. Public FV leaves '
        'are newly observed at the original public PMT amount. The four literal expression/store '
        'forms and direct FV endpoint differences are fixed before those FV answers; no logarithm, private principal or selector is '
        'fitted. A complete survivor would still require fresh discriminants and validation. '
        'Beginning-payment semantics are outside this isolated route test.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def cumulative_denominators(checkout,output,cache,offline=False):
    """Pin a possible first-principal denominator without a growth producer."""
    rates=(.0025,.03,.125,2.**-10);periods=(12.,60.,360.)
    balances=[.125,1.,1.1,1000.,125000.,2.**20]
    balances += [value(f'{int(bits(1.),16)+k:016x}') for k in (-2,-1,1,2)]
    groups=[(r,n) for r in rates for n in periods]
    rows=[[r,n,pv,1.,1.,0.] for r,n in groups for pv in balances]
    write(str(output)+'.questions.json',dict(format='dnacalc-cumulative-denominator-questions/1',complete=True,
        rows=[[bits(x) for x in r] for r in rows],source_bindings=SOURCE_BINDINGS,
        hypothesis='One shared binary64 positive growth-minus-one denominator per loan; quotient-first and product-first principal expressions, each with coherent strict or stored-PC64 arithmetic.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('CUMPRINC',[[client.number(x) for x in row] for row in rows])['rows']
    fv=service.ask('FV',[[client.number(r),client.number(n),client.number(0.),client.number(1.),client.number(0.)] for r,n in groups])['rows']
    powers=service.ask('POWER',[[client.number(1.+r),client.number(n)] for r,n in groups])['rows']
    records=[];counts=Counter()
    for i,(rate,n) in enumerate(groups):
        rr=binary64(bits(rate));observed=answers[i*len(balances):(i+1)*len(balances)]
        candidates=[]
        for order in ('quotient_first','product_first'):
            for model in ('strict','stored_pc64'):
                domain=None;contradicted=False;first_contradiction=None
                for j,(pv,answer) in enumerate(zip(balances,observed)):
                    if answer['type']!='number':raise ValueError('first principal must be numeric')
                    target=-binary64(answer['bits']);pp=binary64(bits(pv))
                    if order=='quotient_first':
                        q=inverse(Span(target,target),rr,'mul',True,model)
                        constraint=inverse(q,pp,'div',False,model) if q is not None else None
                    else:
                        numerator=operation(pp,rr,'mul',model)
                        constraint=inverse(Span(target,target),numerator,'div',False,model)
                    if not contradicted:
                        domain=constraint if domain is None else (None if constraint is None else domain.meet(constraint))
                        if domain is None:contradicted=True;first_contradiction=j
                coords=B.coordinates(domain) if domain else None
                total=0 if coords is None else coords[1]-coords[0]+1
                name=order+'/'+model;counts[name+'/groups']+=1;counts[name+'/nonempty']+=total>0
                counts[name+'/unique']+=total==1
                values=[] if coords is None else [bits(float(B.value(k))) for k in range(coords[0],min(coords[1]+1,coords[0]+8))]
                public={}
                for leaf,growth in (('FV',fv[i]),('POWER',powers[i])):
                    if growth['type']!='number':raise ValueError('ordinary growth must be numeric')
                    g=binary64(growth['bits'])*(-1 if leaf=='FV' else 1)
                    denominator=operation(g,-1,'add',model)
                    public[leaf]=dict(denominator=bits(float(denominator)),possible=domain is not None and domain.contains(denominator))
                    counts[name+'/'+leaf+'_possible']+=public[leaf]['possible']
                candidates.append(dict(order=order,model=model,domain=domain.json() if domain else None,
                    binary64_count=total,values=values,first_contradiction=first_contradiction,public_growth=public))
        records.append(dict(rate=bits(rate),periods=bits(n),balances=[bits(v) for v in balances],
                            answers=observed,public_fv=fv[i],public_power=powers[i],candidates=candidates))
    result=dict(format='dnacalc-cumulative-first-principal-denominators/1',complete=True,counts=counts,
        records=records,receipts=service.receipts,source_bindings=SOURCE_BINDINGS,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Necessary common-denominator constraints from ten independent PV scales/adjacent balances. '
        'A surviving denominator is not a growth producer and is never installed. The public FV/POWER '
        'leaves are independent controls, not assumed private helpers. No finite program is selected '
        'by how many rows it matches; every empty intersection is retained.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def cumulative_first_residual(study,output,cache):
    """Invert the first-period cancellation without prescribing an exponential.

    The unknown em is shared across balances, not chosen per answer. Every
    fixed source form is tested independently; no best-match form is selected.
    """
    from geomean import verified_study_answers
    from exact import EXTENDED64
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete first-principal observations')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    models=dict(stored_assignments=('stored_pc64',)*3,
                single_expression=('resident_pc64',)*3,
                stored_quotient=('stored_pc64','resident_pc64','resident_pc64'),
                stored_quotient_and_interest=('stored_pc64','resident_pc64','stored_pc64'),
                stored_payment=('resident_pc64','stored_pc64','resident_pc64'))
    records=[];counts={name:Counter() for name in models}
    for group in doc['records']:
        rate=binary64(group['rate']);candidates=[]
        for name,(division,product,interest_model) in models.items():
            domain=None;contradicted=False;constraints=[]
            for pv,answer in zip(group['balances'],group['answers']):
                if answer['type']!='number':raise ValueError('nonnumeric original principal')
                p=binary64(pv);target=binary64(answer['bits'])
                interest=operation(p,rate,'mul',interest_model)
                payment=inverse(Span(target,target),interest,'add',True,'stored_pc64')
                quotient=inverse(payment,rate,'mul',True,product) if payment is not None else None
                needed=inverse(quotient,p,'div',False,division) if quotient is not None else None
                constraints.append(needed.json() if needed is not None else None)
                if not contradicted:
                    domain=needed if domain is None else (None if needed is None else domain.meet(needed))
                    if domain is None:contradicted=True
            binary=B.coordinates(domain) if domain is not None else None
            extended=EXTENDED64.coordinates(domain) if domain is not None else None
            counts[name]['groups']+=1
            counts[name]['binary64_nonempty']+=binary is not None
            counts[name]['resident64_nonempty']+=extended is not None
            counts[name]['real_nonempty']+=domain is not None
            candidates.append(dict(name=name,division=division,payment_product=product,
                interest_product=interest_model,required_shared_em=domain.json() if domain is not None else None,
                binary64_count=0 if binary is None else binary[1]-binary[0]+1,
                resident64_count=0 if extended is None else extended[1]-extended[0]+1,
                constraints=constraints))
        records.append(dict(rate=group['rate'],periods=group['periods'],balances=group['balances'],
                            answers=group['answers'],candidates=candidates))
    result=dict(format='dnacalc-cumulative-first-residual-preimages/1',complete=True,counts=counts,
        source_bindings={**SOURCE_BINDINGS,**bindings,str(source.resolve()):sha(source)},records=records,
        original_receipts=doc['receipts'],
        limits='Conditional first-principal source forms: q=pv/em, payment=q*rate, '
        'interest=pv*rate, result=stored-PC64(payment+interest). A shared em is constrained '
        'by all ten balances before lattice projection. Both binary64 and resident64 '
        'possibilities are reported. Nonempty cells are not a growth producer or an '
        'implementation; an empty intersection excludes only its declared store form. '
        'No new Excel reading, fitted helper, or product-code change.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def loss_controls(checkout,requests,output,cache,executable,offline=False):
    """Locate downstream losses using original public payment/component calls."""
    source=Path(requests);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete original loss list')
    frozen={**SOURCE_BINDINGS,str(source.resolve()):sha(source),str(Path(executable).resolve()):sha(executable),
            **doc.get('source_bindings',{})}
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    payments={};components={};records=[]
    for original in doc['records']:
        name=original['function'];args=original['args']
        if name in ('IPMT','PPMT') and len(args)==6:
            pmt=[args[i] for i in (0,2,3,4,5)];parts=[]
        elif name in ('CUMPRINC','CUMIPMT') and len(args)==6:
            pmt=[args[0],args[1],args[2],bits(0.),args[5]]
            start,end=(int(value(args[i])) for i in (3,4))
            if not 1<=start<=end<=10000:raise ValueError('loss interval outside bounded controls')
            parts=[[args[0],bits(k),args[1],args[2],bits(0.),args[5]] for k in range(start,end+1)]
        else:raise ValueError('unsupported downstream loss shape')
        payments[tuple(pmt)]=None
        for p in parts:components[tuple(p)]=None
        records.append(dict(original,payment_args=pmt,component_args=parts))
    # These are independent leaves, never re-asked parent outcomes.
    for name,groups in [('PMT',payments),('PPMT',components)]:
        keys=list(groups)
        if not keys:continue
        answers=service.ask(name,[[client.number(x) for x in row] for row in keys])['rows']
        for key,answer in zip(keys,answers):groups[key]=dict(answer=answer)
        inputs=Path(str(output)+'.'+name+'.inputs.json');actuals=Path(str(output)+'.'+name+'.actuals.jsonl')
        witnesses=[dict(id='control-'+str(i),args=['0x'+b for b in row]) for i,row in enumerate(keys)]
        write(inputs,dict(function=name,witnesses=witnesses))
        # ubs:ignore — fixed argv to the bound local application binary, no shell, bounded execution.
        subprocess.run([str(executable),'invoke-inputs',str(inputs),str(actuals)],check=True,timeout=300)
        computed=[json.loads(s) for s in actuals.read_text().splitlines()]
        if len(computed)!=len(witnesses):raise ValueError('incomplete application controls')
        for key,expected,actual in zip(keys,witnesses,computed):
            if any(expected[k]!=actual[k] for k in ('id','args')):raise ValueError('application control identity changed')
            groups[key]['calculated']=actual['calc_func_prediction']
    counts=Counter()
    for row in records:
        p=payments[tuple(row['payment_args'])];a=p['answer']
        payment_correct=a['type']=='number' and p['calculated']=='0x'+a['bits']
        row['payment']=p;row['payment_matches_original_public_call']=payment_correct
        counts['rows']+=1;counts['correct_payment']+=payment_correct
        if row['component_args']:
            ps=[components[tuple(k)] for k in row['component_args']]
            all_correct=all(p['answer']['type']=='number' and p['calculated']=='0x'+p['answer']['bits'] for p in ps)
            total=0
            for p in ps:
                if p['answer']['type']!='number':raise ValueError('public principal component is not numeric')
                total=operation(total,binary64(p['answer']['bits']),'add','strict')
            published=bits(float(total));row.update(components=ps,public_ppmt_strict_sum=published,
                all_components_match_original_public_calls=all_correct)
            counts['correct_all_components']+=all_correct
            classification=('published_PPMT_sum_differs_from_original_cumulative_route' if all_correct
                else 'published_PPMT_component_discrepancy_in_cumulative_route')
        else:
            classification='downstream_IPMT_or_PPMT_discrepancy_despite_correct_public_payment'
        if not payment_correct:classification='existing_public_payment_helper_discrepancy'
        row['classification']=classification;counts[classification]+=1
    result=dict(format='dnacalc-annuity-downstream-loss-controls/1',complete=True,counts=counts,
        records=records,receipts=service.receipts,source_bindings=frozen,
        limits='The original lost parent answers are never replaced or re-asked. Independent public '
        'PMT and, for cumulative rows, every public PPMT component locate whether a payment helper, '
        'component or cumulative composition still differs. Correct payment alone identifies no '
        'private downstream algorithm. These observations classify remaining work, not gains.')
    write(output,result);print(json.dumps(counts),flush=True);return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout')
    p.add_argument('output')
    p.add_argument('--cache',required=True)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--discriminants')
    p.add_argument('--negative',action='store_true')
    p.add_argument('--compose')
    p.add_argument('--edges',action='store_true')
    p.add_argument('--products')
    p.add_argument('--timing-edges',action='store_true')
    p.add_argument('--cumulative-pairs',action='store_true')
    p.add_argument('--cumulative-fv')
    p.add_argument('--cumulative-denominators',action='store_true')
    p.add_argument('--loss-controls',nargs=2,metavar=('REQUESTS','EXECUTABLE'))
    p.add_argument('--first-residual',metavar='STUDY')
    a=p.parse_args()
    if a.first_residual:
        cumulative_first_residual(a.first_residual,a.output,a.cache)
    elif a.loss_controls:
        loss_controls(a.checkout,a.loss_controls[0],a.output,a.cache,a.loss_controls[1],a.offline)
    elif a.cumulative_denominators:
        cumulative_denominators(a.checkout,a.output,a.cache,a.offline)
    elif a.cumulative_fv:
        cumulative_fv(a.checkout,a.cumulative_fv,a.output,a.cache,a.offline)
    elif a.cumulative_pairs:
        cumulative_pairs(a.checkout,a.output,a.cache,a.offline)
    elif a.timing_edges:
        timing_edges(a.checkout,a.output,a.cache,a.offline)
    elif a.products:
        products(a.checkout,a.products,a.output,a.cache,a.offline)
    elif a.edges:
        edges(a.checkout,a.output,a.cache,a.offline)
    elif a.compose:
        compose(a.checkout,a.compose,a.output,a.cache,a.offline)
    elif a.discriminants:
        discriminants(a.checkout,a.discriminants,a.output,a.cache,a.offline)
    else:
        run(a.checkout,a.output,a.cache,a.offline,a.negative)
