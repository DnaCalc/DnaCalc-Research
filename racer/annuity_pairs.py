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


def due_parts(checkout, output, cache, cohort="discovery", offline=False):
    """Constrain a shared stored period component with paired public calls.

    Ordinary annuity algebra gives the beginning-payment coefficient
    1/(1+1/r). End-payment publications constrain the entire latent component,
    so no logarithm, exponential or private quotient is assumed here.
    """
    import itertools
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    rows=[]
    if cohort not in ("discovery","holdout","validation35"):raise ValueError('unknown cohort')
    holdout=cohort=="holdout"
    rates=(2.**-20,.00375,.0625,.2,.375,math.nextafter(.5,0.)) if holdout else (2.**-10,.05,.1,.125,.25,.5)
    periods=(3.,17.,128.5) if holdout else (4.,12.,360.)
    balances=(math.nextafter(1.,0.),1.1,math.ldexp(1000.,-20),math.ldexp(1000.,20)) if holdout else (1.,math.nextafter(1000.,math.inf),math.ldexp(1000.,-40),math.ldexp(1000.,40))
    if cohort=='validation35':
        rates=(2.**-12,.007,.03125,.1875,.4,.8)
        periods=(5.,29.,96.25)
        balances=(1.25,math.nextafter(300.,math.inf),math.ldexp(750.,-10),math.ldexp(750.,30))
    for r,n,pv in itertools.product(rates,periods,balances):
        for per in sorted(set((2.,2.5 if holdout else n/2.,n-1.,n))):
            for fv in (0.,(-.3 if holdout else .7)*pv):
                rows.append([r,per,n,pv,fv])
    questions=[[client.number(x) for x in row+[timing]] for row in rows for timing in (0.,1.)]
    answers={f:service.ask_rows(f,questions) for f in ('PPMT','IPMT')}
    records=[];counts=Counter()
    for f, observed in answers.items():
        for i,row in enumerate(rows):
            end,due=observed[2*i:2*i+2];r=binary64(bits(row[0]));names={}
            if end['type']!='number' or due['type']!='number':
                raise ValueError('finite positive ordinary annuity must calculate')
            e,d=map(lambda a:binary64(a['bits']),(end,due))
            domain=inverse(Span(e,e),r,'mul',True,'stored_pc64')
            coords=B.coordinates(domain) if domain else None
            if coords is not None and coords[1]-coords[0]>8:
                raise ValueError('shared component is not sufficiently isolated')
            latent=[] if coords is None else [B.value(k) for k in range(coords[0],coords[1]+1)]
            factor=operation(1,operation(1,operation(1,r,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64')
            names['shared_component_due_coefficient']=[operation(t,factor,'mul','stored_pc64') for t in latent]
            names['outer_timing_division']=[operation(e,operation(1,r,'add','stored_pc64'),'div','stored_pc64')]
            names['shared_component_timing_then_rate']=[operation(operation(t,operation(1,r,'add','stored_pc64'),'div','stored_pc64'),r,'mul','stored_pc64') for t in latent]
            possible={k:d in vv for k,vv in names.items()}
            for k,v in possible.items():counts[f+'/'+('fv0' if row[4]==0 else 'fv')+'/'+k+'/possible']+=v
            counts[f+'/rows']+=1
            records.append(dict(function=f,args=list(map(bits,row)),end=end,due=due,
                latent_binary64_domain=domain.json() if domain else None,possible=possible,
                predictions={k:list(map(lambda x:bits(float(x)),vv)) for k,vv in names.items()}))
    write(output,dict(complete=True,cohort=cohort,calls=len(questions)*2,counts=counts,records=records,
        receipts=service.receipts,source_sha256=sha(Path(__file__)),
        limits='Conditional shared stored component only. Every end-call preimage is retained in full. Nonempty domains do not identify a complete period kernel; no observed helper is installed.'))
    print(json.dumps(counts),flush=True)


def due_quotients(checkout, source, output, cache, offline=False):
    """Read PMT independently for each period-component loan."""
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    document=json.loads(Path(source).read_text());loans={}
    for row in document['records']:
        a=row['args'];key=tuple(a[i]for i in (0,2,3,4));loans.setdefault(key,len(loans))
    records=[];questions=[]
    for key,index in loans.items():
        for timing in (0.,1.):questions.append([client.number(value(x))for x in key]+[client.number(timing)])
        records.append(dict(id=index,args=list(key)))
    answers=service.ask_rows('PMT',questions)
    for i,row in enumerate(records):row['end'],row['due']=answers[2*i:2*i+2]
    write(output,dict(complete=True,loans=records,receipts=service.receipts,
        source_sha256={str(Path(source).resolve()):sha(source),str(Path(__file__).resolve()):sha(__file__)},
        limits='Independent public PMT observations only. No private quotient is recovered or installed by these readings.'))
    print(json.dumps(dict(loans=len(records),calls=len(answers))))


def due_quotient_domains(source, components, output):
    """Inverse PMT publications before inspecting a proposed period quotient."""
    doc=json.loads(Path(source).read_text());observed=json.loads(Path(components).read_text())
    def qvalue(s):return binary64(s)
    groups={};current={}
    for row in observed:
        a=row['original']['args'];key=tuple(a[i].removeprefix('0x')for i in (0,2,3,4))
        if row['quotient']is not None:
            if key in current and current[key]!=row['quotient']:raise ValueError('period-dependent current quotient')
            current[key]=row['quotient']
    for row in doc['loans']:
        a=row['args'];r,n,p,f=map(qvalue,a)
        if f!=0:continue
        y=qvalue(row['end']['bits']);qd=inverse(Span(y,y),r,'mul',True,'stored_pc64')
        em=inverse(qd,p,'div',False,'stored_pc64');key=tuple(a[:2])
        groups[key]=em if key not in groups else groups[key].meet(em)
    ems={}
    for k,d in groups.items():
        coords=B.coordinates(d)if d else None
        if coords is None or coords[1]-coords[0]>8:raise ValueError('unresolved independent em domain')
        ems[k]=[B.value(i)for i in range(coords[0],coords[1]+1)]
    records=[];counts=Counter()
    for row in doc['loans']:
        a=row['args'];r,n,p,f=map(qvalue,a);factor=operation(1,operation(1,operation(1,r,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64')
        end,due=map(lambda x:qvalue(row[x]['bits']),('end','due'))
        d0=inverse(Span(end,end),r,'mul',True,'stored_pc64');d1=inverse(Span(due,due),factor,'mul',True,'stored_pc64');joint=d0.meet(d1)
        predictions=[]
        for em in ems[tuple(a[:2])]:
            balance=operation(p,f,'add','stored_pc64');newq=operation(operation(balance,em,'div','stored_pc64'),f,'add','stored_pc64')
            predictions.append(dict(em=bits(float(em)),quotient=bits(float(newq)),admitted=joint.contains(newq)if joint else False))
        old=current.get(tuple(a));oldok=joint.contains(qvalue(old))if joint and old else False
        counts['loans']+=1;counts['unique_independent_em']+=len(ems[tuple(a[:2])])==1;counts['new_fold_admitted']+=any(x['admitted']for x in predictions);counts['current_quotient_admitted']+=oldok
        records.append(dict(original=row,joint_quotient_domain=joint.json()if joint else None,current_quotient=old,current_admitted=oldok,predictions=predictions))
    write(output,dict(complete=True,counts=counts,records=records,
        source_sha256={str(Path(source).resolve()):sha(source),str(Path(components).resolve()):sha(components),str(Path(__file__).resolve()):sha(__file__)},
        limits='Binary64 em cells come only from independent zero-future PMT calls. The fixed common balance fold is evaluated at every remaining em. Neither a nonempty domain nor shared mathematical quotient proves a whole period kernel.'))
    print(json.dumps(counts))


def period_boundaries(checkout, output, cache, cohort="boundary", offline=False):
    """Read caller boundary cases independently of the period-component premise."""
    import itertools
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline);cases=[]
    if cohort not in ('boundary','higher','higher_holdout','guard','guard_holdout'):raise ValueError('unknown cohort')
    rates=(.625,1.5,4.,16.)if cohort=='higher_holdout'else(.75,1.,2.,10.)if cohort=='higher'else(2.**-10,.05,.375,1.,-.125)
    periods=(5.,17.,50.)if cohort=='higher_holdout'else(3.,12.,60.)if cohort=='higher'else(3.,12.,360.)
    if cohort not in ('guard','guard_holdout'):
        for r,n,pv in itertools.product(rates,periods,(.125,math.nextafter(1.,0.),12345.678)if cohort=='higher_holdout'else(1.,math.nextafter(1000.,math.inf),200000.)):
            points=(1.,2.5,n-.25,n)if cohort=='higher_holdout'else(1.,2.5,n/2.,n-1.,n)if cohort=='higher'else(1.,math.nextafter(1.,math.inf),1.5,math.nextafter(n,0.),n)
            for fv,per,timing in itertools.product((0.,-.3*pv,.7*pv,-pv),points,(0.,1.)):
                cases.append(dict(args=[bits(x)for x in(r,per,n,pv,fv,timing)]))
    else:
        import struct
        seen=set()
        def add(per,n):
            args=[bits(x)for x in(.125,per,n,1000.,0.,0.)]
            if tuple(args)not in seen:seen.add(tuple(args));cases.append(dict(args=args))
        if cohort=='guard_holdout':
            small=(2.**-100,2.**-60,2.**-54,2.**-53,2.**-52,.125,.5,.875,1.,1.25,1.5)
            for n in small:
                for per in(math.nextafter(1.,0.),1.,math.nextafter(1.,math.inf),1.125,1.5,
                        n+.5,n+1.,math.nextafter(n+1.,0.),math.nextafter(n+1.,math.inf)):
                    for timing in(0.,1.):
                        args=[bits(x)for x in(.125,per,n,1000.,0.,timing)]
                        if tuple(args)not in seen:seen.add(tuple(args));cases.append(dict(args=args))
            for exponent in(53,54,63,64,65,66,100,1023):
                n=2.**exponent
                for per in(math.nextafter(n,0.),n,math.nextafter(n,math.inf)):
                    for timing in(0.,1.):
                        args=[bits(x)for x in(.125,per,n,1000.,0.,timing)]
                        if tuple(args)not in seen:seen.add(tuple(args));cases.append(dict(args=args))
        for base in (()if cohort=='guard_holdout'else(2.,3.,8.,17.,128.,378.,1024.,2.**20,2.**32,2.**52,2.**53)):
            for n in ((base,base+.5)if base<2.**52 else(base,)):
                raw=int(bits(n),16)
                for delta in range(-4,5):add(value(f'{raw+delta:016x}'),n)
                for delta in (-1.,-.5,0.,.125,.5,.875,1.,1.125):add(n+delta,n)
            raw=int(bits(base),16)
            for delta in range(-4,5):add(base,value(f'{raw+delta:016x}'))
            for delta in (-1.,-.875,-.5,-.125,0.,.5,1.):add(base,base+delta)
    write(str(output)+'.questions.json',dict(cases=cases,source_sha256=sha(Path(__file__)),
        selection=('Fixed period/nper adjacent windows and fractional/exponent boundaries. No tolerance fitted.'if cohort in('guard','guard_holdout')else'Fixed higher positive rates and legal period points, first and opposing-balance exceptions included.'if cohort in('higher','higher_holdout')else'Fixed rates, signs, caller first-period and adjacent-double boundaries, final periods and opposite balances. No period-component premise imposed on first beginning payments or fallback branches.')))
    answers={f:service.ask_rows(f,[[client.number(value(x))for x in c['args']]for c in cases])for f in('IPMT','PPMT')}
    write(output,dict(complete=True,cohort=cohort,cases=cases,answers=answers,receipts=service.receipts,
        source_sha256=sha(Path(__file__)),limits='Caller boundaries and numerical fallback readings only. No expected answer is inferred or substituted for a recorded result.'))
    print(json.dumps(dict(cases=len(cases),calls=sum(map(len,answers.values())))))


if __name__ == '__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='period-boundaries':
        period_boundaries(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='due-quotient-domains':
        due_quotient_domains(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='due-quotients':
        due_quotients(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='due-parts':
        due_parts(*sys.argv[2:]);raise SystemExit(0)


def run(checkout, output, cache, offline=False, negative=False, final_model=None, discount=False):
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    # Fixed ordinary rates, periods and scaled/adjacent balances; no bank
    # failures or current application predictions select these questions.
    rates = (.0025, .0075, .03, .1, .125, .25, .5, 2.**-10)
    periods = (1., 12., 360., 1024.)
    if negative:
        rates=tuple(-r for r in rates)
        periods=(1.,12.,60.,120.)
        if discount:
            periods=(-.5,-1.,-12.,-60.)
    final_model=final_model or ('strict' if negative else 'stored_pc64')
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
        branch=('negative_rate_discount' if discount else 'negative_rate_growth') if negative else 'positive_rate_discount', final_model=final_model,
        groups=groups, factor_groups=factor_groups, receipts=service.receipts, source_bindings=SOURCE_BINDINGS,
        limits='Finite PMT at the declared rate sign and final-product word only. One shared binary64 em per rate/period is constrained by ten independent type-zero balances. Timing and nonzero-future-value stages are tested separately against all remaining em values. Nonempty constraints do not identify the logarithm or expm1; no parameter or selector is fitted. Opposing balances retain their separately measured cancellation route.')
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


def products(checkout,prior,output,cache,offline=False,timing=False):
    """Mine final-product ties with a denominator pinned by earlier public calls."""
    data=json.loads(Path(prior).read_text())
    if data.get('complete') is not True:
        raise ValueError('incomplete denominator study')
    rng=random.Random(2901);rows=[];excluded=[];factor_constraints={}
    if timing:
        for group in data['groups']:
            if len(group['em_values'])!=1:raise ValueError('unresolved denominator')
            em=binary64(group['em_values'][0])
            state=factor_constraints.setdefault(group['rate'],{word:dict(domain=None,dead=False) for word in ('strict','stored_pc64')})
            for observation in group['timing']:
                q=operation(binary64(observation['pv']),em,'div','stored_pc64')
                y=binary64(observation['due']['bits'])
                for word,item in state.items():
                    needed=inverse(Span(y,y),q,'mul',True,word)
                    if not item['dead']:
                        item['domain']=needed if item['domain'] is None else item['domain'].meet(needed) if needed else None
                        item['dead']=item['domain'] is None
        for raw,state in factor_constraints.items():
            r=binary64(raw)
            fixed=operation(1,operation(1,operation(1,r,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64')
            for word,item in state.items():
                domain=item.pop('domain');coords=B.coordinates(domain)if domain else None
                if coords is None or coords[0]!=coords[1] or B.value(coords[0])!=fixed:
                    raise ValueError('final words do not independently constrain the same fixed coefficient')
                item.update(domain=domain.json(),unique_coefficient=bits(float(fixed)))
    for group in data['groups']:
        if len(group['em_values'])!=1:
            raise ValueError('unresolved denominator')
        rate=binary64(group['rate']);em=binary64(group['em_values'][0])
        factor=operation(1,operation(1,operation(1,rate,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64') if timing else rate
        magnitude=abs(factor)
        if magnitude.numerator.bit_count()==1 and magnitude.denominator.bit_count()==1:
            excluded.append(dict(rate=group['rate'],periods=group['periods'],reason='Power-of-two multiplication is exact on these normal operands.'))
            continue
        selected=0
        for trial in range(100000):
            pv=f'{((1023+rng.randint(-10,10))<<52)|rng.getrandbits(52):016x}'
            q=operation(binary64(pv),em,'div','stored_pc64')
            single=operation(q,factor,'mul','strict')
            double=operation(q,factor,'mul','stored_pc64')
            if single==double:
                continue
            rows.append(dict(args=[group['rate'],group['periods'],pv,'0000000000000000',bits(1. if timing else 0.)],
                em=group['em_values'][0],factor=bits(float(factor)),strict=bits(float(single)),stored_pc64=bits(float(double))))
            selected+=1
            if selected==2:
                break
        if selected!=2:
            if not timing:
                raise ValueError('no complete final-product discriminant set')
            excluded.append(dict(rate=group['rate'],periods=group['periods'],selected=selected,trials=100000,
                reason='The fixed bounded primitive search did not find two product-word discriminants. No mathematical impossibility is claimed.'))
    write(str(output)+'.questions.json',dict(format='dnacalc-annuity-product-questions/1',complete=True,
        rows=rows,excluded=excluded,prior_sha256=sha(prior),source_bindings=SOURCE_BINDINGS,
        timing=timing,factor_constraints=factor_constraints,
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
        timing=timing,factor_constraints=factor_constraints,
        questions_sha256=sha(str(output)+'.questions.json'),source_bindings=SOURCE_BINDINGS,
        limits='One final multiplication, fv0, at the declared timing, with independently pinned binary64 em and a stored-PC64 quotient. '
               'For beginning payments both final-word models independently constrain the same unique binary64 coefficient from earlier public calls before new ties are selected. '
               'Neither results refute the conditional earlier-stage model; no result changes a denominator or selector. Power-of-two factors cannot distinguish the two final roundings here.')
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


def cumulative_power_denominators(checkout,output,cache,offline=False):
    """Question geometric principal numerators with any common loan denominator."""
    import itertools
    from exact import EXTENDED64 as E
    rates=sorted({f'{int(bits(center),16)-4+i:016x}'for center in (.0025,.03125,2.**-10)for i in range(8)})
    spans=((1,1),(2,2),(4,4),(1,2),(1,4),(2,5),(1,8))
    cases=[dict(rate=r,periods=n,pv=bits(pv),start=s,end=e)
           for r,n,pv,(s,e)in itertools.product(rates,(8,16,36),(1.,1.1,128.,10000.),spans)]
    profiles=list(itertools.product(('difference','factored'),('strict','stored_pc64','resident_pc64'),
                                   ('strict','stored_pc64','resident_pc64'),('strict','stored_pc64')))
    def label(p):return '/'.join(p)
    bases={}
    for r in rates:
        base=operation(binary64(bits(1.)),binary64(r),'add','strict')
        if base!=operation(binary64(bits(1.)),binary64(r),'add','stored_pc64'):
            raise ValueError('stored base word differs in this fixed cohort')
        bases[r]=bits(float(base))
    powers=sorted({(bases[c['rate']],k)for c in cases for k in(c['start']-1,c['end'],c['end']-c['start']+1,c['periods'])})
    frozen=dict(SOURCE_BINDINGS)
    write(str(output)+'.questions.json',dict(format='dnacalc-cumulative-power-denominator-questions/1',complete=True,
        cases=cases,powers=powers,profiles=profiles,source_bindings=frozen,
        premise='End-payment principal over s..e is mathematically -PV*(q^e-q^(s-1))/(q^n-1), q=1+r. '
                'Compare direct/factored numerators using independently published POWER leaves. '
                'Relax the denominator to any one resident positive value per fixed rate/n, shared across PV and intervals. '
                'No PMT/FV call or private payment is assumed. Stored q is a declared premise, not an identified private base.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    leaves=service.ask('POWER',[[client.number(q),client.number(k)]for q,k in powers])['rows']
    if any(r['type']!='number'for r in leaves):raise ValueError('finite positive power leaves required')
    lookup={key:binary64(row['bits'])for key,row in zip(powers,leaves)}
    args=[[c['rate'],bits(c['periods']),c['pv'],bits(c['start']),bits(c['end']),bits(0.)]for c in cases]
    answers=service.ask('CUMPRINC',[[client.number(v)for v in row]for row in args])['rows']
    states={};records=[]
    for index,(case,answer)in enumerate(zip(cases,answers)):
        if answer['type']!='number'or binary64(answer['bits'])>=0:
            raise ValueError('negative principal reading required')
        r,n=case['rate'],case['periods'];q=bases[r]
        first,last=lookup[(q,case['start']-1)],lookup[(q,case['end'])]
        count=lookup[(q,case['end']-case['start']+1)]
        target=-binary64(answer['bits']);pv=binary64(case['pv'])
        for profile in profiles:
            form,dword,mword,fword=profile;key=(r,n,label(profile))
            state=states.setdefault(key,dict(domain=None,applied=0,contradiction=None))
            if state['contradiction']is not None:continue
            delta=(operation(last,-first,'add',dword)if form=='difference'else
                   operation(first,operation(count,-binary64(bits(1.)),'add',dword),'mul',dword))
            numerator=operation(pv,delta,'mul',mword)
            needed=inverse(Span(target,target),numerator,'div',False,fword)
            state['domain']=needed if state['applied']==0 else state['domain'].meet(needed)if state['domain']and needed else None
            state['applied']+=1
            if state['domain']is None or E.coordinates(state['domain'])is None:
                state['contradiction']=dict(index=index,case=case,answer=answer,
                    numerator=str(numerator),last_denominator_cell=needed.json()if needed else None)
        records.append(dict(case=case,args=args[index],answer=answer))
    models=[];feasible={label(p):True for p in profiles};loan_feasible=Counter()
    for (r,n,name),state in states.items():
        domain=state['domain'];coords=E.coordinates(domain)if domain else None
        feasible[name]&=coords is not None
        if coords:loan_feasible[(r,n)]+=1
        candidates={word:operation(lookup[(bases[r],n)],-binary64(bits(1.)),'add',word)
                    for word in('strict','stored_pc64','resident_pc64')}
        models.append(dict(rate=r,periods=n,profile=name,applied_constraints=state['applied'],
            contradiction=state['contradiction'],domain=domain.json()if domain else None,
            resident_count=coords[1]-coords[0]+1 if coords else 0,
            contains_public_power_denominator={word:bool(domain and domain.contains(v))for word,v in candidates.items()}))
    result=dict(format='dnacalc-cumulative-power-denominators/1',complete=True,records=records,models=models,
        counts=dict(rows=len(records),loans=len(rates)*3,loans_with_every_profile_refuted=len(rates)*3-len(loan_feasible)),
        complete_profiles=[name for name,ok in feasible.items()if ok],power_leaves=[dict(base=q,exponent=k,answer=a)for(q,k),a in zip(powers,leaves)],
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Necessary exact denominator constraints under declared published POWER and stored-base premises. '
               'An empty resident intersection refutes that whole numerator/word program for the loan; '
               'a surviving interval identifies no private denominator and never scores as an application pass. '
               'Private resident bases, unpublished powers, different numerators and beginning-payment routes remain open.')
    write(output,result);print(json.dumps(dict(counts=result['counts'],complete_profiles=result['complete_profiles'])),flush=True)
    return result


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
    payments={};components={};future_leaves={};records=[]
    for original in doc['records']:
        name=original['function'];args=original['args']
        if name=='FV' and len(args)==5:
            r,n=args[:2]
            leaves=[[r,n,bits(0.),bits(-1.),bits(0.)],
                    [r,n,bits(-1.),bits(0.),bits(0.)],
                    [r,n,bits(-1.),bits(0.),bits(1.)]]
            for p in leaves:future_leaves[tuple(p)]=None
            records.append(dict(original,future_leaf_args=leaves))
            continue
        elif name in ('IPMT','PPMT') and len(args)==6:
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
    for name,groups in [('PMT',payments),('PPMT',components),('FV',future_leaves)]:
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
        if row['function']=='FV':
            args=row['args'];leaves=[future_leaves[tuple(k)] for k in row['future_leaf_args']]
            row['future_leaves']=leaves;counts['rows']+=1
            if all(p['answer']['type']=='number' for p in leaves):
                growth=binary64(leaves[0]['answer']['bits'])
                term=binary64(leaves[2 if value(args[4])!=0 else 1]['answer']['bits'])
                present=operation(binary64(args[3]),growth,'mul','stored_pc64')
                payment=operation(binary64(args[2]),term,'mul','stored_pc64')
                prediction='0x'+bits(float(-operation(present,payment,'add','stored_pc64')))
                row['public_leaf_composition']=prediction
                row['public_leaf_composition_matches_original']=prediction==row['expected']
                failed=[i for i,p in enumerate(leaves)if p['calculated']!='0x'+p['answer']['bits']]
                row['mismatched_leaf_indices']=failed
                classification=('FV_growth_or_annuity_production_discrepancy' if failed
                    else 'FV_composition_discrepancy_despite_correct_public_leaves')
                if prediction!=row['expected']:classification='FV_public_leaf_composition_counterexample'
            else:classification='FV_leaf_domain_requires_separate_review'
            row['classification']=classification;counts[classification]+=1
            continue
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
        'private downstream algorithm. FV unit growth/end/beginning leaves separately test '
        'its declared stored two-product/addition composition. These observations classify '
        'remaining work, not gains.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def growth_reciprocal_edges(checkout, output, cache, offline=False):
    """Pair payment/parts across exact power-of-two growth overflow edges."""
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    cases=[]
    for rate,center in ((-.5,1024.),(-.75,512.),(value('bfefffffffffffff'),1024./53.)):
        first=int(bits(center),16)
        cases.extend(dict(rate=bits(rate),periods=f'{first+k:016x}')for k in range(-3,4))
    questions={name:[]for name in ('PMT','IPMT','PPMT')}
    rows={name:[]for name in questions}
    for i,c in enumerate(cases):
        for timing in (0.,1.):
            args=[c['rate'],c['periods'],bits(3.7),bits(.9),bits(timing)]
            questions['PMT'].append([client.number(a)for a in args]);rows['PMT'].append(dict(case=i,timing=timing,args=args))
            for period in (1.,2.):
                aa=[c['rate'],bits(period),c['periods'],bits(3.7),bits(.9),bits(timing)]
                for name in ('IPMT','PPMT'):
                    questions[name].append([client.number(a)for a in aa]);rows[name].append(dict(case=i,timing=timing,period=period,args=aa))
    write(str(output)+'.questions.json',dict(format='dnacalc-annuity-growth-reciprocal-questions/1',complete=True,
        cases=cases,questions=questions,rows=rows,source_bindings=SOURCE_BINDINGS,
        selection='Seven neighboring periods at the exact mathematical 2^1024 growth boundary '
                  'for three dyadic growth bases. Both timing values and first/second parts; '
                  'first beginning IPMT zero is the only proposed exception to payment refusals.'))
    observed={name:service.ask(name,q)['rows']for name,q in questions.items()}
    payment={(r['case'],r['timing']):a for r,a in zip(rows['PMT'],observed['PMT'])}
    records=[];counts=Counter()
    for name in questions:
        for row,answer in zip(rows[name],observed[name]):
            expected='not_asserted'
            if name=='IPMT'and row['period']==1 and row['timing']==1:
                expected='zero'
            elif name!='PMT'and payment[(row['case'],row['timing'])]==dict(type='error',code='#NUM!'):
                expected='error:Num'
            seen='zero'if answer['type']=='number'and int(answer['bits'],16)&((1<<63)-1)==0 else 'error:Num'if answer==dict(type='error',code='#NUM!')else 'other'
            counts[name+'/calls']+=1;counts[name+'/num']+=seen=='error:Num'
            counts[name+'/constrained']+=expected!='not_asserted'
            counts[name+'/contradictions']+=expected!='not_asserted'and expected!=seen
            records.append(dict(function=name,**row,answer=answer,expected_admission=expected,observed_admission=seen))
    result=dict(format='dnacalc-annuity-growth-reciprocal-edges/1',complete=True,records=records,
        counts=counts,receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),source_bindings=SOURCE_BINDINGS,
        limits='Adjacent period observations and paired public-domain relations only. No exact '
               'private overflow threshold or exponential producer is inferred. Numerical values '
               'remain controls, not an identified payment graph. A shared PMT refusal must '
               'preserve first-period beginning IPMT zero after finite/period-range validation.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def cumulative_component_constraints(checkout, study, prior_cache, output, cache, offline=False):
    """Test unrounded singleton components and paired cumulative returns.

    Hidden components are relaxed to complete exact publication cells. Their
    sum ranges overapproximate every permitted component, so exclusion is a
    necessary refutation, never a fitted private principal value.
    """
    from geomean import verified_study_answers
    from exact import EXTENDED64
    source = Path(study)
    doc = json.loads(source.read_text())
    _, bindings, _ = verified_study_answers(doc,Path(prior_cache))
    frozen = {**SOURCE_BINDINGS, **bindings, str(source.resolve()):sha(source)}
    key = lambda row: tuple(row[k] for k in ('rate','periods','pv','timing'))
    points = {key(row)+(int(value(row['period'])),):row for row in doc['singletons']}
    profiles = [dict(component=component,accumulator=model)
                for component in ('published_binary64','resident_before_publication','real_before_pc64_publication')
                for model in ('strict','stored_pc64','resident_pc64')]
    counts = [Counter() for _ in profiles]
    rows = []
    questions = []
    for row in doc['spans']:
        singletons = [points[key(row)+(k,)] for k in range(row['start'],row['end']+1)]
        answer = binary64(row['answer']['bits'])
        predictions = []
        for index, profile in enumerate(profiles):
            total = Span(0,0)
            for point in singletons:
                observed = binary64(point['cumprinc']['bits'])
                component = Span(observed,observed)
                if profile['component'] != 'published_binary64':
                    component = B.pullback(component)
                    if profile['component'] == 'real_before_pc64_publication':
                        component = EXTENDED64.pullback(component)
                # Closed rounded endpoint enclosure is a deliberate relaxation.
                low = round_model(total.low+component.low,profile['accumulator'])
                high = round_model(total.high+component.high,profile['accumulator'])
                total = Span(low,high)
            published = Span(B.round(total.low),B.round(total.high))
            possible = published.contains(answer)
            counts[index]['spans'] += 1
            counts[index]['possible'] += possible
            counts[index]['excluded'] += not possible
            predictions.append(dict(possible=possible,published_enclosure=published.json()))
        args = [row['rate'],row['periods'],row['pv'],bits(row['start']),bits(row['end']),row['timing']]
        questions.append([dict(type='number',bits=a) for a in args])
        rows.append(dict(args=args,answer=row['answer'],payment=singletons[0]['payment'],
            singleton_principals=[point['cumprinc'] for point in singletons],component_models=predictions))
    write(str(output)+'.questions.json',dict(format='dnacalc-cumulative-component-questions/1',complete=True,
        questions=questions,profiles=profiles,rows=rows,source_bindings=frozen,
        selection='The648 retained short intervals are unchanged. New paired CUMIPMT answers '
                  'test count*PMT-CUMPRINC separately from all hidden-component sum bounds.'))
    client,_ = load_service(checkout)
    service = CachedService(client.Client(timeout=300),cache,offline)
    answers = service.ask('CUMIPMT',questions)['rows']
    paired_counts = Counter()
    for row, answer in zip(rows,answers):
        if answer['type'] != 'number':
            raise ValueError('ordinary cumulative interest must be numeric')
        length = binary64(row['args'][4])-binary64(row['args'][3])+1
        payment = binary64(row['payment']['bits'])
        principal = binary64(row['answer']['bits'])
        row['interest'] = answer
        row['paired_predictions'] = {}
        for multiply in ('strict','stored_pc64','resident_pc64'):
            for subtract in ('strict','stored_pc64'):
                product = operation(length,payment,'mul',multiply)
                predicted = bits(float(operation(product,-principal,'add',subtract)))
                name = multiply+'/'+subtract
                row['paired_predictions'][name] = predicted
                paired_counts[name+'/matched'] += answer == dict(type='number',bits=predicted)
                paired_counts[name+'/rows'] += 1
    result = dict(format='dnacalc-cumulative-component-constraints/1',complete=True,
        profiles=profiles,counts=counts,paired_counts=paired_counts,records=rows,
        source_bindings=frozen,receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Every inherited singleton and principal interval answer is unchanged. '
               'Hidden components are constrained only by the original singleton publication, '
               'then all rounded-sum endpoint bounds are conservatively closed. '
               'Exclusion rules out even this relaxation of independent, interval-invariant components. '
               'A surviving interval identifies no component or route. Paired count*PMT-CUMPRINC '
               'compositions are fixed before new CUMIPMT answers; none is selected by counts.')
    write(output,result)
    print(json.dumps(dict(component_counts=counts,paired_counts=paired_counts)),flush=True)
    return result


def future_value_stages(checkout,output,cache,offline=False):
    """Separate FV's public unit leaves, products and final balance sum."""
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    groups=[dict(rate=bits(r),periods=bits(n))
            for r in (.0025,.03,.125,1.,-.0025,-.1,-.5,value('bfeffffffffffffe'))
            for n in (1.,2.,12.,60.,-.5,-1.,-12.)]
    leaf_questions=[]
    for g in groups:
        for payment,present,timing in ((0.,-1.,0.),(-1.,0.,0.),(-1.,0.,1.)):
            leaf_questions.append([client.number(v) for v in (g['rate'],g['periods'],bits(payment),bits(present),bits(timing))])
    write(str(output)+'.leaf-questions.json',dict(format='dnacalc-fv-unit-leaf-questions/1',complete=True,
        groups=groups,questions=leaf_questions,source_bindings=SOURCE_BINDINGS,
        selection='Eight fixed rates, seven fixed periods; growth and both payment-timing leaves are separate public unit calls, before product probes.'))
    leaves=service.ask('FV',leaf_questions)['rows']
    rng=random.Random(3107);probes=[];exclusions=[];questions=[]
    for i,g in enumerate(groups):
        g['leaves']=leaves[3*i:3*i+3]
        for leaf_index,leaf in enumerate(g['leaves']):
            if leaf['type']!='number':
                exclusions.append(dict(group=i,leaf=leaf_index,reason='Original unit leaf is nonnumeric; no arithmetic prediction.'))
                continue
            factor=binary64(leaf['bits'])
            balances=[bits(v) for v in (.125,1.,1.1,3.7,1000.,-1.1)]
            magnitude=abs(factor)
            exact_scale=not magnitude or (magnitude.numerator.bit_count()==1 and magnitude.denominator.bit_count()==1)
            found=0
            if not exact_scale:
                for trial in range(50000):
                    raw=f'{((1023+rng.randint(-10,10))<<52)|rng.getrandbits(52):016x}'
                    b=binary64(raw)
                    try:
                        strict=operation(b,factor,'mul','strict');stored=operation(b,factor,'mul','stored_pc64')
                    except ValueError:continue
                    if strict==stored:continue
                    balances.append(raw);found+=1
                    if found==2:break
                if found!=2:exclusions.append(dict(group=i,leaf=leaf_index,reason='Bounded50000-draw primitive search did not find two ties; no impossibility claim.',selected=found))
            for balance in balances:
                args=[g['rate'],g['periods'],bits(0.) if leaf_index==0 else balance,
                      balance if leaf_index==0 else bits(0.),bits(1.) if leaf_index==2 else bits(0.)]
                predictions={word:bits(float(-operation(binary64(balance),factor,'mul',word)))
                             for word in ('strict','stored_pc64')}
                probes.append(dict(group=i,leaf=leaf_index,args=args,balance=balance,factor=leaf,predictions=predictions))
                questions.append([client.number(a) for a in args])
    write(str(output)+'.product-questions.json',dict(format='dnacalc-fv-product-questions/1',complete=True,
        probes=probes,exclusions=exclusions,receipts=service.receipts,source_bindings=SOURCE_BINDINGS,
        selection='Six fixed balances and the first two strict/PC64 product ties per non-power-of-two public unit factor, seed3107; primitive-only search before consumer answers.'))
    answers=service.ask('FV',questions)['rows'];counts=Counter();domains={}
    for row,answer in zip(probes,answers):
        row['answer']=answer
        for word,prediction in row['predictions'].items():
            counts[word+'/calls']+=1;counts[word+'/matched']+=answer==dict(type='number',bits=prediction)
            key=(row['group'],row['leaf'],word)
            state=domains.setdefault(key,dict(dead=False,domain=None))
            if answer['type']!='number':state['dead']=True;continue
            needed=inverse(Span(-binary64(answer['bits']),-binary64(answer['bits'])),binary64(row['balance']),'mul',True,word)
            if not state['dead']:
                state['domain']=needed if state['domain'] is None else state['domain'].meet(needed) if needed else None
                state['dead']=state['domain'] is None
    inferred=[]
    for (group,leaf,word),state in domains.items():
        interval=state['domain'] if not state['dead'] else None
        coords=B.coordinates(interval) if interval else None
        inferred.append(dict(group=group,leaf=leaf,word=word,domain=interval.json() if interval else None,
            binary64_count=0 if coords is None else coords[1]-coords[0]+1))
    result=dict(format='dnacalc-fv-isolated-products/1',complete=True,groups=groups,records=probes,
        counts=counts,inferred_factor_domains=inferred,exclusions=exclusions,receipts=service.receipts,source_bindings=SOURCE_BINDINGS,
        questions_sha256=sha(str(output)+'.product-questions.json'),
        limits='Public unit factors are observed controls, not assumed private returns. Product '
               'predictions test their binary64 reuse. Necessary factor intervals allow any '
               'common underlying factor across each fixed group and final-word model. Empty '
               'domains refute only that product composition. Nonnumeric leaves remain counted '
               'exclusions, not arithmetic passes. No factor, tree or domain threshold is fitted.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def future_value_coupled(checkout,study,output,cache,offline=False):
    """Check the isolated FV stores when both balances are nonzero."""
    import itertools
    from fractions import Fraction as Q
    from exact import exponent_of,two
    from geomean import verified_study_answers
    source=Path(study);doc=json.loads(source.read_text())
    _,bindings,_=verified_study_answers(doc,Path(cache))
    frozen={**SOURCE_BINDINGS,**bindings,str(source.resolve()):sha(source)}
    profiles=[dict(present_product=a,payment_product=b,addition=c)
        for a,b,c in itertools.product(('strict','stored_pc64','resident_pc64'),
            ('strict','stored_pc64','resident_pc64'),('strict','stored_pc64'))]
    rows=[];missing=[]
    def add(gi,timing,present,payment,role):
        g=doc['groups'][gi];f,t=(binary64(g['leaves'][j]['bits']) for j in (0,1+timing))
        p,m=binary64(present),binary64(payment);predictions=[]
        for profile in profiles:
            left=operation(p,f,'mul',profile['present_product'])
            right=operation(m,t,'mul',profile['payment_product'])
            value=-operation(left,right,'add',profile['addition'])
            predictions.append(bits(float(value)))
        rows.append(dict(group=gi,role=role,args=[g['rate'],g['periods'],payment,present,bits(float(timing))],predictions=predictions))
    for gi,g in enumerate(doc['groups']):
        for timing in (0,1):
            if any(g['leaves'][j]['type']!='number' for j in (0,1+timing)):
                missing.append(dict(group=gi,timing=timing,reason='Nonnumeric original leaf.'));continue
            f,t=(binary64(g['leaves'][j]['bits']) for j in (0,1+timing))
            for present,payment in ((1.1,.7),(-1.1,.7),(.125,1000.),(-3.7,1000.)):
                add(gi,timing,bits(present),bits(payment),'fixed_nonzero_pair')
            p=bits(1.1);left=operation(binary64(p),f,'mul','stored_pc64')
            if not left or not t:continue
            unit=two(exponent_of(abs(left))-52)
            for sign in (1,-1):
                ideal=Q(sign)*unit/(2*t);center=int(bits(float(ideal)),16);selected=None
                for offset in range(-64,65):
                    payment=f'{center+offset:016x}';right=operation(binary64(payment),t,'mul','stored_pc64')
                    a=operation(left,right,'add','strict');b=operation(left,right,'add','stored_pc64')
                    if a!=b:selected=payment;break
                if selected is None:
                    missing.append(dict(group=gi,timing=timing,sign=sign,
                        reason='No addition-word tie in the fixed129-adjacent-payment window.'))
                else:add(gi,timing,p,selected,'first_primitive_addition_word_tie')
    write(str(output)+'.questions.json',dict(format='dnacalc-fv-coupled-questions/1',complete=True,
        profiles=profiles,rows=rows,excluded=missing,source_bindings=frozen,
        selection='Four fixed nonzero PV/payment pairs and the first final-addition strict/PC64 tie in each signed129-double payment window around half an ULP of the independently known PV product. Before fresh answers.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('FV',[[client.number(a) for a in row['args']] for row in rows])['rows']
    counts=[Counter() for _ in profiles]
    for row,answer in zip(rows,answers):
        row['answer']=answer
        for i,prediction in enumerate(row['predictions']):
            counts[i]['calls']+=1;counts[i]['matched']+=answer==dict(type='number',bits=prediction)
            counts[i]['failed']+=answer!=dict(type='number',bits=prediction)
    result=dict(format='dnacalc-fv-coupled-stores/1',complete=True,profiles=profiles,
        counts=counts,records=rows,excluded=missing,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='The original independent unit factors and isolated product evidence are '
               'unchanged. These18 literal two-product/addition programs check their '
               'composition on new nonzero pairs, with final-addition ties selected using '
               'primitive arithmetic only. Missing ties are bounded exclusions. No model '
               'is ranked or selected by count; all original errors remain failures.')
    write(output,result);print(json.dumps(list(zip(profiles,counts))),flush=True);return result


def future_value_timing(checkout,study,output,cache,offline=False):
    """Keep FV's timing-factor production separate from its final stores."""
    import itertools
    from geomean import verified_study_answers
    source=Path(study);doc=json.loads(source.read_text())
    _,bindings,_=verified_study_answers(doc,Path(cache))
    frozen={**SOURCE_BINDINGS,**bindings,str(source.resolve()):sha(source)}
    profiles=[dict(form=form,timing_factor=tf,inner=inner,finish=finish,difference='strict')
        for form,tf,inner,finish in itertools.product(
            ('multiply_divide','divide_multiply','coefficient_multiply'),
            ('strict','resident_pc64'),('strict','stored_pc64','resident_pc64'),('strict','stored_pc64'))]
    first={};alive=set(range(len(profiles)))
    for gi,g in enumerate(doc['groups']):
        r=binary64(g['rate']);u=operation(binary64(g['leaves'][0]['bits']),-1,'add','strict')
        for timing in (0,1):
            target=g['leaves'][1+timing]
            for i in sorted(alive):
                p=profiles[i];tf=operation(1,r,'add',p['timing_factor']) if timing else 1
                if p['form']=='multiply_divide':v=operation(operation(tf,u,'mul',p['inner']),r,'div',p['finish'])
                elif p['form']=='divide_multiply':v=operation(tf,operation(u,r,'div',p['inner']),'mul',p['finish'])
                else:v=operation(operation(tf,r,'div',p['inner']),u,'mul',p['finish'])
                prediction=bits(float(v))
                if target!=dict(type='number',bits=prediction):
                    alive.remove(i);first[i]=dict(group=gi,timing=timing,target=target,prediction=prediction)
    rng=random.Random(3119);selected=[];excluded=[]
    for stage in ('end_division','timing_product'):
        for sign in (0,1):
            for exponent in (-20,-10,-4,-1):
                found=0
                for trial in range(100000):
                    raw=f'{(sign<<63)|((1023+exponent)<<52)|rng.getrandbits(52):016x}'
                    r=binary64(raw);tf=operation(1,r,'add','strict');u=tf-1
                    strict=operation(u,r,'div','strict');stored=operation(u,r,'div','stored_pc64')
                    discriminates=(strict!=stored if stage=='end_division' else
                        operation(strict,tf,'mul','strict')!=operation(strict,tf,'mul','stored_pc64'))
                    if not discriminates:continue
                    selected.append(dict(rate=raw,stage=stage,sign=sign,exponent=exponent,trial=trial))
                    found+=1
                    if found==2:break
                if found!=2:excluded.append(dict(stage=stage,sign=sign,exponent=exponent,selected=found,trials=100000,
                    reason='The fixed primitive search was incomplete; no impossibility claim.'))
    selected.extend(dict(rate=bits(sign*2.**exponent),stage='fixed_small_rate_control')
                    for sign in (1,-1)for exponent in (-53,-52,-40))
    questions=[]
    for row in selected:
        r=row['rate']
        questions.extend([[client_value for client_value in (r,bits(1.),bits(pmt),bits(pv),bits(timing))]
                          for pmt,pv,timing in ((0.,-1.,0.),(-1.,0.,0.),(-1.,0.,1.),(.7,1.1,0.),(.7,1.1,1.))])
    write(str(output)+'.questions.json',dict(format='dnacalc-fv-timing-questions/1',complete=True,
        selected=selected,questions=questions,excluded=excluded,profiles=profiles,
        earlier_survivors=sorted(alive),earlier_first_contradictions=first,source_bindings=frozen,
        selection='Thirty-six literal source forms are constrained by all168 original unit leaves. Seed3119 selects first two n=1 primitive end-division or timing-product ties per signed rate binade, before new readings. Six fixed small-rate controls are separate. No model selected by count.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('FV',[[client.number(v) for v in q] for q in questions])['rows']
    counts=Counter();records=[]
    for i,row in enumerate(selected):
        observed=answers[5*i:5*i+5];r=binary64(row['rate']);tf=operation(1,r,'add','strict')
        counts['groups']+=1
        predictions={}
        if all(a['type']=='number' for a in observed[:3]):
            f,t,d=(binary64(a['bits'])for a in observed[:3])
            counts['growth_equals_stored_one_plus_rate']+=f==tf
            for word in ('strict','stored_pc64'):
                end=operation(f-1,r,'div',word);due=operation(t,tf,'mul',word)
                predictions[word]=dict(end=bits(float(end)),due_from_observed_end=bits(float(due)))
                counts[word+'/end_matches']+=t==end
                counts[word+'/due_from_observed_end_matches']+=d==due
        records.append(dict(row,questions=questions[5*i:5*i+5],answers=observed,predictions=predictions))
    result=dict(format='dnacalc-fv-timing-production/1',complete=True,profiles=profiles,
        earlier_survivors=sorted(alive),earlier_first_contradictions=first,
        records=records,counts=counts,excluded=excluded,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='A separate timing-helper study, not an application change. Earlier source forms '
               'use the independently observed unit growth and an explicitly stored binary64 '
               'difference growth-minus-one. Fresh n=1 hypotheses check their '
               'growth and end-term premises against separate public FV leaves, then compare '
               'the beginning term to those observed leaves. Fixed small-rate controls can '
               'refute the generic quotient formula and are never dropped. Bounded searches '
               'do not establish impossibility or an approximation cutoff. All original '
               'mixed calls remain available for the registered application replay.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def future_value_growth(checkout,timing,controls,executable,output,cache,offline=False):
    """Isolate FV's base addition before the existing public POWER route."""
    from geomean import verified_study_answers
    paths=list(map(Path,(timing,controls,executable)))
    documents=[json.loads(p.read_text())for p in paths[:2]]
    frozen={**SOURCE_BINDINGS,**{str(p.resolve()):sha(p)for p in paths}}
    for document in documents:
        _,bindings,_=verified_study_answers(document,Path(cache));frozen.update(bindings)
    base_counts=Counter()
    for row in documents[0]['records']:
        r=binary64(row['rate']);answer=row['answers'][0]
        for word in ('strict','stored_pc64'):
            predicted=bits(float(operation(1,r,'add',word)))
            base_counts[word+'/calls']+=1
            base_counts[word+'/matches']+=answer==dict(type='number',bits=predicted)
    questions=[];rows=[]
    for row in documents[1]['records']:
        r,n=row['args'][:2];growth=row['future_leaves'][0]['answer']
        for word in ('strict','stored_pc64'):
            base=bits(float(operation(1,binary64(r),'add',word)))
            rows.append(dict(parent_id=row['id'],rate=r,periods=n,word=word,base=base,observed_FV_growth=growth))
            questions.append([base,n])
    write(str(output)+'.questions.json',dict(format='dnacalc-fv-growth-questions/1',complete=True,
        rows=rows,base_counts=base_counts,source_bindings=frozen,
        hypothesis='A stored-PC64 one-plus-rate base, then the unchanged public POWER kernel. The32 earlier n=1 growth observations test base formation directly. All14 remaining FV parents provide independent public growth leaves; parent outcomes are never re-asked.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('POWER',[[client.number(v)for v in q]for q in questions])['rows']
    inputs=Path(str(output)+'.inputs.json');actuals=Path(str(output)+'.actuals.jsonl')
    witnesses=[dict(id='power-'+str(i),args=['0x'+v for v in q])for i,q in enumerate(questions)]
    write(inputs,dict(function='POWER',witnesses=witnesses))
    # ubs:ignore — fixed argv to the bound local executable, no shell, original numeric inputs.
    subprocess.run([str(executable),'invoke-inputs',str(inputs),str(actuals)],check=True,timeout=300)
    calculated=[json.loads(s)for s in actuals.read_text().splitlines()]
    if len(calculated)!=len(rows):raise ValueError('incomplete registered POWER controls')
    counts=Counter()
    for row,answer,request,actual in zip(rows,answers,witnesses,calculated):
        if any(request[k]!=actual[k]for k in ('id','args')):raise ValueError('POWER control identity changed')
        row['public_POWER']=answer;row['registered_POWER']=actual['calc_func_prediction']
        expected=row['observed_FV_growth'];prefix=row['word']
        counts[prefix+'/calls']+=1
        counts[prefix+'/public_POWER_matches_FV_growth']+=answer==expected
        counts[prefix+'/registered_POWER_matches_FV_growth']+=expected['type']=='number'and actual['calc_func_prediction']=='0x'+expected['bits']
    result=dict(format='dnacalc-fv-growth-base-production/1',complete=True,base_counts=base_counts,
        counts=counts,records=rows,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='One base-addition stage and the already existing POWER producer. '
               'Earlier n=1 growth calls directly expose the base; residual growth and '
               'POWER are independent public leaves. Registered POWER is a calculated '
               'control, never Excel evidence. This identifies no replacement power '
               'algorithm, coefficient or cutoff. Every remaining parent stays a failure '
               'until the actual FV route is changed and remeasured.')
    write(output,result);print(json.dumps(dict(base_counts=base_counts,power_counts=counts)),flush=True);return result


def cumulative_aggregate_constraints(study, prior, cache, prior_cache, output):
    """Can one unpublished interest aggregate finish both cumulative calls?"""
    from exact import BINARY64 as B, EXTENDED64 as E, inverse, pull_rounding
    from geomean import verified_study_answers
    source, before = Path(study), Path(prior)
    doc, old = json.loads(source.read_text()), json.loads(before.read_text())
    answers, bindings, _ = verified_study_answers(doc, Path(cache))
    inherited, old_bindings, _ = verified_study_answers(old, Path(prior_cache))
    answers.update(inherited)
    frozen = {**SOURCE_BINDINGS, **bindings, **old_bindings,
              str(source.resolve()):sha(source), str(before.resolve()):sha(before)}
    profiles = [dict(payment_product=product, principal_return=word, interest=hidden)
                for product in ('strict','stored_pc64','resident_pc64')
                for word in ('strict','stored_pc64')
                for hidden in ('binary64','resident_pc64','real_before_pc64')]
    counts = [Counter() for _ in profiles]
    records = []
    for row in doc['records']:
        args = row['args']
        payment_args = [args[0],args[1],args[2],bits(0.),args[5]]
        for name, request, expected in (
                ('PMT',payment_args,row['payment']), ('CUMPRINC',args,row['answer']),
                ('CUMIPMT',args,row['interest'])):
            key = name,json.dumps([dict(type='number',bits=x) for x in request],sort_keys=True)
            if answers.get(key) != expected:
                raise ValueError('aggregate constraint lost an original public reading')
        n = binary64(args[4])-binary64(args[3])+1
        payment = binary64(row['payment']['bits'])
        principal = binary64(row['answer']['bits'])
        interest = binary64(row['interest']['bits'])
        domains = []
        for i, profile in enumerate(profiles):
            product = operation(n,payment,'mul',profile['payment_product'])
            negative = inverse(Span(principal,principal),product,'add',True,profile['principal_return'])
            required = negative.linear(-1) if negative else None
            observed = (pull_rounding(Span(interest,interest),'stored_pc64')
                        if profile['interest']=='real_before_pc64' else B.pullback(Span(interest,interest)))
            domain = observed.meet(required) if required else None
            lattice = B if profile['interest']=='binary64' else E
            possible = domain is not None and (profile['interest']=='real_before_pc64' or lattice.coordinates(domain) is not None)
            counts[i]['spans'] += 1
            counts[i]['possible'] += possible
            counts[i]['excluded'] += not possible
            domains.append(dict(possible=possible,required_interest=domain.json() if domain else None))
        records.append(dict(args=args,payment=row['payment'],principal=row['answer'],
                            interest=row['interest'],models=domains))
    result = dict(format='dnacalc-cumulative-aggregate-constraints/1',complete=True,
        profiles=profiles,counts=counts,records=records,source_bindings=frozen,
        common_exclusions=sum(not any(p['possible'] for p in row['models']) for row in records),
        limits='Fixed final aggregate expression principal = count*public_PMT - hidden_interest. '
               'The same hidden interest must publish as the original CUMIPMT. Three payment '
               'product precisions, two principal return words and three hidden-interest domains '
               'are stated independently of matches. Empty cells refute only this final aggregate '
               'with the public payment; they do not refute a private payment, interval-specific '
               'loops or other cumulative routes. A possible cell identifies no hidden value '
               'and never counts as a calculated application pass. No new Excel answers.')
    write(output,result)
    print(json.dumps(dict(counts=counts,common_exclusions=result['common_exclusions'])),flush=True)
    return result


def reciprocal_domain(checkout, output, cache, offline=False, transfer=False):
    """Discriminate the PC64 reciprocal pole before cancellation shortcuts.

    The reciprocal is rounded independently with rational arithmetic. The
    first beginning-payment interest is a separate public-domain question.
    No numerical payment result is predicted by this admission experiment.
    """
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    first = 0xbff0000000000000
    rates = [f'{first-i:016x}' for i in range(17)]
    rates += [bits(r) for r in (-.75, -.5, -1./3., -.25, .125)]
    if transfer:
        rates += [bits(r) for r in (-2., -1.25, 0.)]
    questions = {name: [] for name in ('PMT', 'IPMT', 'PPMT', 'FV')}
    metadata = {name: [] for name in questions}
    for raw in rates:
        r = binary64(raw)
        reciprocal = operation(binary64(bits(1.)), r, 'div', 'stored_pc64') if r else None
        coefficient_denominator = operation(binary64(bits(1.)), reciprocal, 'add', 'stored_pc64') if r else None
        pole = coefficient_denominator == 0
        for n in (.5, 1., 2., 17., 97.):
            for pv, fv in ((3.7, 0.), (3.7, -3.7), (3.7, .9), (0., 3.7)):
                for timing in (0., 1.):
                    args = [raw, bits(n), bits(pv), bits(fv), bits(timing)]
                    questions['PMT'].append([client.number(a) for a in args])
                    metadata['PMT'].append(dict(args=args, reciprocal_pole=pole,
                        expected_admission='error:Num' if r <= -1 or (pole and (not transfer or timing == 1)) else 'not_asserted'))
        for n in (2., 17., 97.):
            for per in (1., 2.):
                for timing in (0., 1.):
                    args = [raw, bits(per), bits(n), bits(3.7), bits(.9), bits(timing)]
                    for name in ('IPMT', 'PPMT'):
                        questions[name].append([client.number(a) for a in args])
                        exception = name == 'IPMT' and per == 1 and timing == 1 and (transfer or r > -1)
                        rejected = r <= -1 or (pole and (not transfer or timing == 1))
                        expected = 'zero' if transfer and exception else ('zero' if exception else 'error:Num') if rejected else 'not_asserted'
                        metadata[name].append(dict(args=args, reciprocal_pole=pole,
                            first_interest_exception=exception, expected_admission=expected))
            args = [raw, bits(n), bits(.9), bits(3.7), bits(1.)]
            questions['FV'].append([client.number(a) for a in args])
            metadata['FV'].append(dict(args=args, reciprocal_pole=pole,
                expected_admission='not_asserted'))
    manifest = dict(format='dnacalc-annuity-reciprocal-domain-questions/1', complete=True,
        questions=questions, metadata=metadata, source_bindings=SOURCE_BINDINGS,
        phase=2 if transfer else 1,
        selection='Seventeen adjacent rates at -1 and five fixed ordinary controls; fixed periods, balances and timing. '
                  'The PC64 reciprocal-pole hypothesis precedes every fresh answer. '
                  'Phase two separates beginning timing after phase one counterexamples and adds -2, -1.25, zero. '
                  'First-period beginning-payment IPMT zero is tested before the rate-domain checks.')
    write(str(output)+'.questions.json', manifest)
    counts = {name: Counter() for name in questions}
    records = []
    for name, rows in questions.items():
        answers = service.ask(name, rows)['rows']
        for m, answer in zip(metadata[name], answers):
            expected = m['expected_admission']
            observed = ('zero' if answer['type'] == 'number' and int(answer['bits'], 16) & ((1<<63)-1) == 0
                        else 'error:Num' if answer == dict(type='error', code='#NUM!') else 'other')
            counts[name]['calls'] += 1
            counts[name]['constrained'] += expected != 'not_asserted'
            counts[name]['contradictions'] += expected != 'not_asserted' and expected != observed
            counts[name]['pole_calls'] += m['reciprocal_pole']
            records.append(dict(function=name, **m, answer=answer, observed_admission=observed))
    result = dict(format='dnacalc-annuity-reciprocal-domain/1', complete=True, records=records,
        phase=2 if transfer else 1,
        counts=counts, receipts=service.receipts, source_bindings=SOURCE_BINDINGS,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Admission evidence only, from exact input readbacks. '
               'The reciprocal and coefficient denominator use rational PC64 then binary64 rounding. '
               'Ordinary numerical values and FV are retained controls, not identified arithmetic. '
               'A PMT guard must not suppress the separately observed first-period IPMT zero. '
               'No row, coefficient or threshold is fitted to a match count.')
    write(output, result)
    print(json.dumps(counts), flush=True)
    return result


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout')
    p.add_argument('output')
    p.add_argument('--cache',required=True)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--discriminants')
    p.add_argument('--negative',action='store_true')
    p.add_argument('--negative-pc64',action='store_true')
    p.add_argument('--negative-discount',action='store_true')
    p.add_argument('--growth-reciprocal-edges',action='store_true')
    p.add_argument('--compose')
    p.add_argument('--edges',action='store_true')
    p.add_argument('--products')
    p.add_argument('--timing-products',action='store_true')
    p.add_argument('--timing-edges',action='store_true')
    p.add_argument('--cumulative-pairs',action='store_true')
    p.add_argument('--cumulative-fv')
    p.add_argument('--cumulative-denominators',action='store_true')
    p.add_argument('--cumulative-power-denominators',action='store_true')
    p.add_argument('--loss-controls',nargs=2,metavar=('REQUESTS','EXECUTABLE'))
    p.add_argument('--first-residual',metavar='STUDY')
    p.add_argument('--reciprocal-domain', action='store_true')
    p.add_argument('--reciprocal-transfer', action='store_true')
    p.add_argument('--cumulative-components',nargs=2,metavar=('STUDY','PRIOR_CACHE'))
    p.add_argument('--cumulative-aggregate',nargs=3,metavar=('STUDY','PRIOR','PRIOR_CACHE'))
    p.add_argument('--future-stages',action='store_true')
    p.add_argument('--future-coupled')
    p.add_argument('--future-timing')
    p.add_argument('--future-growth',nargs=3,metavar=('TIMING','CONTROLS','EXECUTABLE'))
    a=p.parse_args()
    if a.cumulative_power_denominators:
        cumulative_power_denominators(a.checkout,a.output,a.cache,a.offline)
    elif a.future_growth:
        future_value_growth(a.checkout,*a.future_growth,a.output,a.cache,a.offline)
    elif a.future_timing:
        future_value_timing(a.checkout,a.future_timing,a.output,a.cache,a.offline)
    elif a.future_coupled:
        future_value_coupled(a.checkout,a.future_coupled,a.output,a.cache,a.offline)
    elif a.future_stages:
        future_value_stages(a.checkout,a.output,a.cache,a.offline)
    elif a.cumulative_aggregate:
        cumulative_aggregate_constraints(a.cumulative_aggregate[0],a.cumulative_aggregate[1],a.cache,a.cumulative_aggregate[2],a.output)
    elif a.growth_reciprocal_edges:
        growth_reciprocal_edges(a.checkout,a.output,a.cache,a.offline)
    elif a.cumulative_components:
        cumulative_component_constraints(a.checkout,a.cumulative_components[0],a.cumulative_components[1],a.output,a.cache,a.offline)
    elif a.reciprocal_domain or a.reciprocal_transfer:
        reciprocal_domain(a.checkout,a.output,a.cache,a.offline,a.reciprocal_transfer)
    elif a.first_residual:
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
        products(a.checkout,a.products,a.output,a.cache,a.offline,a.timing_products)
    elif a.edges:
        edges(a.checkout,a.output,a.cache,a.offline)
    elif a.compose:
        compose(a.checkout,a.compose,a.output,a.cache,a.offline)
    elif a.discriminants:
        discriminants(a.checkout,a.discriminants,a.output,a.cache,a.offline)
    else:
        run(a.checkout,a.output,a.cache,a.offline,a.negative or a.negative_pc64 or a.negative_discount,
            'stored_pc64' if a.negative_pc64 or a.negative_discount else None,a.negative_discount)
