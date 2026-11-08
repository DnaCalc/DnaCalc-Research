"""Round40: executable retained-answer analysis and frozen, timed questions.

Original MIT research from Excel answers, IEEE754 and the existing mathematical
tables. No restricted implementation sources. R39 files remain immutable.
"""
import argparse
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import random
import time

from . import ets
from .ets import bits, value, read, grouped, sequential_total, initial_trend
from .fingerprint import CachedService, load_service, sha
from .question_designer import Census

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('/home/ubuntu/codex-runs/round40')
RETAINED = ROOT / 'evidence/round39/retained'


def save(name, document):
    path = OUT / (name + '.json')
    with path.open('x') as stream:
        json.dump(document, stream, separators=(',', ':'), allow_nan=False)
    return path


def frozen(name, document):
    document['producer_source'] = Path(__file__).read_text()
    document['source_sha256'] = sha(__file__)
    return save(name + '-frozen', document)


def service(offline=False):
    client, _ = load_service(ets.CHECKOUT)

    class TimedClient(client.Client):
        def ask(self, function, rows, **kwargs):
            start = time.monotonic()
            record = dict(function=function, calls=len(rows), result_index=kwargs.get('result_index'))
            try:
                reply = super().ask(function, rows, **kwargs)
                record.update(id=reply['id'], ok=reply.get('ok'), excel_build=reply.get('excel_build'))
                return reply
            except Exception as error:
                record.update(error=type(error).__name__, detail=str(error))
                raise
            finally:
                record['seconds'] = time.monotonic() - start
                with (OUT/'service-times.jsonl').open('a') as stream:
                    stream.write(json.dumps(record) + '\n')
                print('Excel wait', function, len(rows), round(record['seconds'], 3), flush=True)

    return client, CachedService(None if offline else TimedClient(timeout=300), OUT/'service-cache', offline)


def retained(name):
    path = RETAINED / (name + '.json.gz')
    return read(path)


def reanalysis():
    from .ets_seasonal import moving_initial
    grid = retained('step7-gamma-structured-fresh-frozen')['grid']
    rounded = [tuple(float(format(v, '.6g')) for v in row) for row in grid]
    representation = [dict(before=row, after=list(other)) for row, other in zip(grid, rounded) if list(other) != row]
    results = []
    for name in ['step7-initialization-fresh', 'step7-initialization-window-fresh', 'step7-initialization-partial-fresh']:
        document = retained(name + '-parameters')
        for i, (case, stat) in enumerate(zip(document['cases'], grouped(document))):
            y, period = case['y'], case['season']
            init = moving_initial(y, period)
            n = len(y)
            start = n - (n if n <= 10 else min(int(.4*n), 10))
            def score(par):
                a, b, g = par
                level, trend, seasons = init[0], init[1], init[2].copy()
                errors = [0.]
                for j, v in enumerate(y[1:], 1):
                    q = j % period
                    error = v - ((level + trend) + seasons[q])
                    errors.append(error)
                    if start and j >= start:
                        level += trend
                    else:
                        level += trend + a * error
                        trend += b * error
                        seasons[q] += g * error
                return sequential_total(v*v for v in errors)
            actual = [value(stat[str(k)]['bits']) for k in [1, 2, 3]]
            predicted = min(rounded, key=score)
            results.append(dict(source=name, case=i, predicted=predicted, actual=actual, matched=list(predicted)==actual))
    census = Census()
    # Reproduce the claim on exactly the original distinct reported parameter
    # observations, retaining post-fit changes rather than treating them as starts.
    for path in sorted(RETAINED.glob('step*-stat.json.gz')):
        document = read(path)
        if not all(k in document for k in ['cases', 'labels', 'answers']): continue
        for i, row in enumerate(grouped(document)):
            for k in [1, 2, 3]:
                a = row.get(str(k))
                if a and a.get('type') == 'number': census.observe([path.name, i, k], a['bits'])
    report = census.report()
    decimals = []
    for raw in report['counts']:
        x = value(json.loads(raw))
        decimals.append(dict(bits=bits(x), value=x, six_significant=float(format(x, '.6g'))==x,
                             nonseasonal_gamma=x==2.**-52))
    save('step1-reanalysis', dict(rounded_grid=rounded, representation_changes=representation,
         original=retained('step7-gamma-generalization-refutations')['matched'],
         rounded_matched=sum(row['matched'] for row in results), total=len(results), rows=results,
         census=report, parameters=decimals,
         limits='Six significant decimal digits describe the reported values, not an identified rounding instruction. The fixed 25-point structure still has counterexamples.'))
    print('rounded family', sum(row['matched'] for row in results), '/', len(results), 'representation changes', len(representation), flush=True)


def quantile_mapping():
    from .normal_quantile import CENTRAL_N, CENTRAL_D, tail
    rng = random.Random(400201)
    points = set()
    for center in [.5, .625, .75, .875, .9, .925, .9375, .0625, .075, .02425, .001, 2.**-26, 2.**-100, 2.**-1022]:
        raw = int(bits(center), 16)
        for offset in range(-64, 65):
            p = value(f'{raw+offset:016x}')
            if 2.**-1022 <= p < 1: points.add(p)
    points.update(rng.uniform(.075, .9375) for _ in range(512))
    points.update(math.exp(rng.uniform(-700, math.log(.075))) for _ in range(512))
    points = sorted(points)
    def as241(p):
        q=p-.5
        if abs(q) > .425: return tail(p)
        r=.180625-q*q
        def poly(coefficients):
            v=coefficients[0]
            for c in coefficients[1:]: v=v*r+c
            return v
        return q*poly(CENTRAL_N)/poly(CENTRAL_D)
    predictions = [bits(as241(p)) for p in points]
    frozen('step2-normal-map', dict(points=points, predictions=predictions,
           design='Fixed AS241 shape; adjacent input-bit regions, exponent edges, independent random central and log-uniform lower points.',
           helper_source=(ROOT/'racer/normal_quantile.py').read_text()))
    client, oracle = service()
    answers = oracle.ask_rows('NORM.S.INV', [[client.number(p)] for p in points], chunk_size=2000)
    save('step2-normal-map-answers', dict(points=points, answers=answers, receipts=oracle.receipts))
    print('normal map', len(points), 'AS241 counterexamples', sum(p!=a.get('bits') for p,a in zip(predictions,answers)), flush=True)


def quantile_refinement():
    """Fixed Newton/Halley trees with conditional public CDF/density leaves."""
    mapping=read(OUT/'step2-normal-map-frozen.json')
    # Small independent range sample of the already frozen input design.
    chosen=list(range(0,len(mapping['points']),5))
    points=[mapping['points'][i] for i in chosen]
    starts=[value(mapping['predictions'][i]) for i in chosen]
    frozen('step2-normal-refinement', dict(points=points, starts=[bits(x)for x in starts],
           structures=['Newton', 'Halley'], cdf_leaves=['public CDF', 'public ERF', 'public ERFC', 'host ERF', 'host ERFC'],
           limits='Conditional leaves identify or refute the consumer structure; complete proof needs an independent leaf implementation.'))
    c,s=service()
    cdf=s.ask_rows('NORM.S.DIST',[[c.number(x),c.boolean(True)]for x in starts],chunk_size=2000)
    phi=s.ask_rows('PHI',[[c.number(x)]for x in starts],chunk_size=2000)
    # Both common typed expressions for x/sqrt(2) are fixed in advance.
    arguments=[[abs(x)/math.sqrt(2), abs(x)*math.sqrt(.5)]for x in starts]
    erf=s.ask_rows('ERF.PRECISE',[[c.number(x)]for row in arguments for x in row],chunk_size=2000)
    erfc=s.ask_rows('ERFC.PRECISE',[[c.number(x)]for row in arguments for x in row],chunk_size=2000)
    answers=read(OUT/'step2-normal-map-answers.json')['answers']
    predictions={}
    for leaf in ['cdf','erf_div','erf_mul','erfc_div','erfc_mul','host_erf','host_erfc']:
        for method in ['Newton','Halley']:
            pred=[]
            for i,(p,x)in enumerate(zip(points,starts)):
                sign=math.copysign(1.,x)
                if leaf=='cdf':error=value(cdf[i]['bits'])-p
                elif leaf.startswith('erf_'):
                    v=value(erf[2*i+(leaf.endswith('mul'))]['bits'])
                    error=sign*.5*v-(p-.5)
                elif leaf.startswith('erfc_'):
                    v=value(erfc[2*i+(leaf.endswith('mul'))]['bits'])/2
                    error=v-p if x<0 else (1-p)-v
                elif leaf=='host_erf':error=.5*math.erf(x/math.sqrt(2))-(p-.5)
                else:
                    v=.5*math.erfc(abs(x)/math.sqrt(2));error=v-p if x<0 else (1-p)-v
                density=value(phi[i]['bits'])
                u=error/density if density else 0.
                result=x-u if method=='Newton'else x-u/(1+x*u/2)
                pred.append(bits(result))
            failures=[dict(index=chosen[i],p=points[i],initial=bits(starts[i]),prediction=v,answer=answers[chosen[i]])for i,v in enumerate(pred)if v!=answers[chosen[i]].get('bits')]
            predictions[leaf+'-'+method]=dict(selections=len(points),failures=len(failures),first=failures[:8])
    save('step2-normal-refinement-proof',dict(results=predictions,cdf=cdf,phi=phi,erf=erf,erfc=erfc,receipts=s.receipts))
    print('normal refinement refutations', {k:v['failures']for k,v in predictions.items()},flush=True)


def confidence_reanalysis():
    from .ets_confidence import variance
    from .ets import numeric_fit
    parameters=grouped(retained('step2-final-recurrence-parameters'))
    coupling=retained('step6-quantile-coupling-answers')
    normal=retained('step6-quantile-coupling-normal')
    z_by_conf=dict(zip(normal['confidences'],[value(a['bits'])for a in normal['answers']]))
    rows=[];counts=Counter()
    for (i,variant),answer in zip(coupling['labels'],coupling['answers']):
        case=coupling['cases'][i];y=case['y'];h,conf=variant
        a,b=[value(parameters[i][str(k)]['bits'])for k in [1,2]]
        _,_,_,errors,_=numeric_fit(y,a,b,case['trend'])
        sigma=sequential_total(e*e for e in errors)/len(y)
        prediction=bits(z_by_conf[conf]*math.sqrt(sigma*variance(a,b,h,39)))
        classification='matched'if prediction==answer.get('bits')else'error_sum_or_initial_state'if h==1 else'fit_or_variance'
        counts[classification]+=1
        rows.append(dict(case=i,n=len(y),alpha=a,beta=b,horizon=h,confidence=conf,
                         prediction=prediction,answer=answer,classification=classification))
    save('step2-confidence-reanalysis',dict(counts=dict(counts),rows=rows,
         conditional_inputs='Original independently measured Excel NORM.S.INV and STAT leaves; exact constructed initial trends.',
         limits='At horizon 1 the variance factor is exactly 1. Remaining losses are in the error sum or hidden initial-state arithmetic; no altered confidence expression is justified. This conditional check does not fix the product quantile.'))
    print('conditional confidence',dict(counts),flush=True)


def ask_series(c, s, cases, function, variants, chunk_size=400):
    """Existing ETS question construction with a throughput-sized batch."""
    rows=[];labels=[]
    for i,z in enumerate(cases):
        y=c.number_array([[v]for v in z['y']]);t=c.number_array([[float(j+1)]for j in range(len(z['y']))]);season=c.number(float(z.get('season',0)))
        for v in variants:
            if function=='FORECAST.ETS.STAT':args=[y,t,c.number(float(v)),season]
            elif function=='FORECAST.ETS':args=[c.number(len(z['y'])+v),y,t,season]
            elif function=='FORECAST.ETS.CONFINT':
                h,conf=v;args=[c.number(len(z['y'])+h),y,t,c.number(conf),season]
            else:args=[y,t]
            rows.append(args);labels.append([i,v])
    answers=s.ask_rows(function,rows,chunk_size=chunk_size)
    return dict(cases=cases,function=function,labels=labels,answers=answers,receipts=s.receipts)


def slope_mapping():
    rng=random.Random(400301);cases=[]
    for n in range(3,11):
        for base in [1., math.pi, 1/3]:
            raw=int(bits(base),16)
            for offset in range(-12,13):
                c=value(f'{raw+offset:016x}')
                cases.append(dict(y=[c]*n,season=0,kind='constant',offset=offset))
        for k in range(20):
            exponent=rng.randrange(-200,201)
            if k<10:y=[math.ldexp(rng.randrange(-2048,2049)/32,exponent)for _ in range(n)]
            else:
                raw=int(bits(math.ldexp(math.pi,exponent)),16)
                y=[value(f'{raw+rng.randrange(-16,17):016x}')for _ in range(n)]
            cases.append(dict(y=y,season=0,kind='short_general'if k<10 else'adjacent'))
    frozen('step3-slope-map',dict(cases=cases, horizons=[0.,1.,2.,5.],
           structures=['raw moments x=0', 'raw moments x=1', 'centered OLS', 'centered x only', 'endpoint', 'mean differences', 'LINEST']))
    c,s=service()
    for name,function,variants in [('parameters','FORECAST.ETS.STAT',[1,2,3]),('forecast','FORECAST.ETS',[0.,1.,2.,5.]),('confidence','FORECAST.ETS.CONFINT',[(1.,.5),(1.,.95)])]:
        save('step3-slope-map-'+name,ask_series(c,s,cases,function,variants))
    # The array-selection service is exercised on the same histories, selecting
    # the coefficient rather than trusting the old stale-loop diagnosis.
    args=[[c.number_array([[v]for v in z['y']]),c.number_array([[float(j+1)]for j in range(len(z['y']))]),c.boolean(True),c.boolean(False)]for z in cases]
    linest=s.ask_rows('LINEST',args,result_index=(1,1),chunk_size=400)
    save('step3-slope-map-linest',dict(cases=cases,answers=linest,receipts=s.receipts))
    print('slope map',len(cases),'histories',flush=True)


def fit_with_trend(y,a,b,trend):
    """The shared R39 recurrence and its measured post-fit reset expression."""
    level=y[0];seasons=[0.,0.];trends=[];errors=[0.]
    for j,v in enumerate(y[1:],1):
        q=j%2;error=v-((level+trend)+seasons[q])
        level+=trend+a*error;trend+=b*error;seasons[q]+=2.**-52*error
        trends.append(trend);errors.append(error)
    window=len(y)if len(y)<=10 else(len(y)+1)//2
    tail=trends[-(window-1):];mean=sequential_total(tail)/len(tail)
    sd=math.sqrt(sequential_total((v-mean)**2 for v in tail)/len(tail))
    reset=abs(mean)<sd*window/len(tail)
    return level,trend,seasons,errors,reset


def slope_analysis():
    from .ets import SLOPES, slope
    cases=read(OUT/'step3-slope-map-frozen.json')['cases']
    parameters=grouped(read(OUT/'step3-slope-map-parameters.json'))
    forecasts=grouped(read(OUT/'step3-slope-map-forecast.json'))
    linest=read(OUT/'step3-slope-map-linest.json')['answers']
    names=[*SLOPES,'endpoint','differences','LINEST'];results={}
    for name in names:
        failures=[];by_n=Counter()
        for i,(case,stat,fo)in enumerate(zip(cases,parameters,forecasts)):
            y=case['y'];a,b=[value(stat[str(k)]['bits'])for k in [1,2]]
            if name in SLOPES:t=slope(y,name)
            elif name=='endpoint':t=(y[-1]-y[0])/(len(y)-1)
            elif name=='differences':t=sequential_total(v-u for u,v in zip(y,y[1:]))/(len(y)-1)
            else:t=value(linest[i]['bits'])
            level,trend,seasons,errors,reset=fit_with_trend(y,a,b,t)
            if reset:level=y[-1]-seasons[(len(y)-1)%2];trend=t
            center=(seasons[0]+seasons[1])/2;level+=center;seasons=[v-center for v in seasons]
            for h in [0.,1.,2.,5.]:
                pred=bits((level+h*trend)+seasons[(len(y)-1+int(h))%2])
                if pred!=fo[str(h)].get('bits'):
                    failures.append(dict(case=i,horizon=h,prediction=pred,answer=fo[str(h)],initial_trend=bits(t),reset=reset))
                    by_n[len(y)]+=1
        results[name]=dict(failures=len(failures),by_n=dict(by_n),first=failures[:12])
    save('step3-slope-first-family',dict(results=results,limits='Reported beta can be post-fit reset. Constant and n<7 cases are not affected by a hidden nonzero beta reset in the ordinary retained scope. Failures refute full candidate programs; measured STAT leaves make them conditional.'))
    print('slope first family refutations', {k:v['failures']for k,v in results.items()},flush=True)


def slope_impulses():
    cases=[]
    for n in range(3,11):
        for at in [-1,*range(n)]:
            y=[128.*j for j in range(n)]
            if at>=0:y[at]+=1/16
            cases.append(dict(y=y,season=0,kind='initial_trend_impulse',at=at))
    frozen('step3-impulses',dict(cases=cases,
           predictions=[bits(initial_trend(z['y']))for z in cases],
           method='Positive exact dyadic trend plus one impulse. Read back initial trend from affine forecast response; the inverse is diagnostic, not a bit-exact hidden value.'))
    c,s=service()
    save('step3-impulses-parameters',ask_series(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]))
    save('step3-impulses-forecast',ask_series(c,s,cases,'FORECAST.ETS',[0.,1.,2.,5.]))


def impulse_analysis():
    cases=read(OUT/'step3-impulses-frozen.json')['cases']
    parameters=grouped(read(OUT/'step3-impulses-parameters.json'))
    forecasts=grouped(read(OUT/'step3-impulses-forecast.json'))
    from .ets import numeric_forecast
    rows=[]
    for i,(case,st,fo)in enumerate(zip(cases,parameters,forecasts)):
        y=case['y'];a,b=[value(st[str(k)]['bits'])for k in [1,2]]
        # The fixed recurrence is affine in the initial slope. Two distant
        # supplied slopes avoid dividing a cancellation-sized difference.
        lo=0.;hi=256.
        f0=numeric_forecast(y,a,b,5.,lo);f1=numeric_forecast(y,a,b,5.,hi)
        observed=value(fo['5.0']['bits'])
        recovered=(observed-f0)/(f1-f0)*(hi-lo)
        rows.append(dict(case=i,n=len(y),at=case['at'],parameters=[a,b],
                         approximate_initial_trend=recovered,raw_ols=initial_trend(y),
                         predicted_forecast=bits(numeric_forecast(y,a,b,5.,initial_trend(y))),observed=fo['5.0']))
    save('step3-impulse-readback',dict(rows=rows,limits='Approximate inverse response, conditioned on the recurrence and measured STAT. No exact hidden slope is inferred from rounded forecasts.'))
    for n in range(3,11):
        base=next(r for r in rows if r['n']==n and r['at']==-1)['approximate_initial_trend']
        print(n,'weights',[(r['at'],round((r['approximate_initial_trend']-base)*16,12))for r in rows if r['n']==n and r['at']>=0],flush=True)


def slope_stores():
    import itertools
    import numpy as np
    from .ets_loop import slope_graph,SLOPE_TAGS
    from .x87_programs import Native
    cases=read(OUT/'step3-slope-map-frozen.json')['cases']
    stat=grouped(read(OUT/'step3-slope-map-parameters.json'))
    forecasts=grouped(read(OUT/'step3-slope-map-forecast.json'))
    # A reported beta != .001 excludes the post-fit beta reset. This removes
    # a conditional-input ambiguity, not a failed observation from statistics.
    selected=[i for i,z in enumerate(cases)if value(stat[i]['2']['bits'])!=.001]
    native=Native(OUT/'x87-round40.so');groups={}
    for n in range(3,11):
        ids=[i for i in selected if len(cases[i]['y'])==n]
        if ids:groups[n]=(ids,np.array([[int(bits(v),16)for v in cases[i]['y']]for i in ids],dtype='uint64'))
    results=[];survivors=[]
    for form,cw,mask in itertools.product(range(3),[0x123f,0x133f],range(1<<len(SLOPE_TAGS))):
        failure=None
        for n,(ids,inputs)in groups.items():
            slopes=native.grid(slope_graph(n,form,True),inputs,cw,mask)[:,0]
            for i,raw in zip(ids,slopes):
                y=cases[i]['y'];a,b=[value(stat[i][str(k)]['bits'])for k in [1,2]]
                t=value(f'{int(raw):016x}');level,trend,seasons,errors,reset=fit_with_trend(y,a,b,t)
                center=(seasons[0]+seasons[1])/2;level+=center;seasons=[v-center for v in seasons]
                for h in [0.,1.,2.,5.]:
                    prediction=bits((level+h*trend)+seasons[(n-1+int(h))%2])
                    if prediction!=forecasts[i][str(h)]['bits']:
                        failure=dict(case=i,horizon=h,initial_slope=bits(t),prediction=prediction,answer=forecasts[i][str(h)])
                        break
                if failure:break
            if failure:break
        profile=dict(form=form,cw=cw,mask=mask)
        if failure:results.append(dict(**profile,first=failure))
        else:survivors.append(profile)
    save('step3-slope-store-refutations',dict(selected=selected,survivors=survivors,refutations=results,
         tags=SLOPE_TAGS,helper_source=(ROOT/'racer/ets_loop.py').read_text(),
         limits='Conditional on the measured recurrence; reported non-reset beta values only. Original excluded reset observations remain in the full mapping and every full-program loss remains counted.'))
    print('slope store programs',len(results)+len(survivors),'survivors',survivors[:8],flush=True)


def padded_history(y):
    t=initial_trend(y)
    return y+[y[j%2]+(j-j%2)*t for j in range(len(y),7)]


def exact_short(y):
    """Sufficient arithmetic proof for both moment regressions and extension."""
    from fractions import Fraction as Q
    def exact(x):return math.isfinite(float(x))and Q(float(x))==x
    def trend(v):
        n=len(v);sx=Q(n*(n+1),2);sxx=Q(n*(n+1)*(2*n+1),6)
        sy=Q(0);sxy=Q(0)
        for i,x in enumerate(v):
            product=Q(i+1)*Q(x)
            sy+=Q(x);sxy+=product
            if not all(map(exact,[product,sy,sxy])):return None
        left=n*sxy;right=sx*sy;num=left-right;den=n*sxx-sx*sx;t=num/den
        return float(t)if all(map(exact,[left,right,num,t]))else None
    t=trend(y)
    if t is None:return False
    v=list(y)
    for j in range(len(y),7):
        product=Q(j-j%2)*Q(t);total=Q(y[j%2])+product
        if not exact(product)or not exact(total):return False
        v.append(float(total))
    return trend(v)is not None


def short_predictions(cases):
    from .ets_optimizer import parameters
    from .ets import statistics_from_errors
    results=[]
    for z in cases:
        y=z['y'];extended=padded_history(y);a,b=parameters(extended)
        t=initial_trend(extended)
        level,trend,seasons,errors,reset=fit_with_trend(y,a,b,t)
        predictions=[y[0]]
        # Prediction j is y[j]-error[j], but that subtraction can lose a bit.
        l=y[0];tr=t;se=[0.,0.]
        for j,v in enumerate(y[1:],1):
            q=j%2;p=(l+tr)+se[q];e=v-p;predictions.append(p);l+=tr+a*e;tr+=b*e;se[q]+=2.**-52*e
        stats=statistics_from_errors(y,errors,predictions)
        if reset:level=y[-1]-seasons[(len(y)-1)%2];trend=t;b=.001
        center=(seasons[0]+seasons[1])/2;level+=center;seasons=[v-center for v in seasons]
        forecasts=[bits((level+h*trend)+seasons[(len(y)-1+int(h))%2])for h in [0.,1.,2.,5.]]
        results.append(dict(stat=[bits(v)for v in [a,b,2.**-52,*stats]],forecasts=forecasts,
                            alpha=a,beta=b,errors=errors,trend=t,reset=reset))
    return results


def short_proof():
    rng=random.Random(400351);cases=[];attempts=Counter()
    for n in range(3,8):
        while sum(len(z['y'])==n for z in cases)<64:
            attempts[n]+=1
            if attempts[n]>200000:raise ValueError('exact-history design exhausted')
            y=[rng.randrange(-2048,2049)/16 for _ in range(n)]
            if not exact_short(y):continue
            exponent=rng.randrange(-200,201);y=[math.ldexp(v,exponent)for v in y]
            if not exact_short(y):raise ValueError('dyadic scale changed exactness')
            cases.append(dict(y=y,season=0,kind='fresh_exact_short'))
    predictions=short_predictions(cases)
    frozen('step3-short-proof',dict(cases=cases,predictions=predictions,attempts=dict(attempts),
           method='Fixed first-two-phase extension to seven using raw OLS; both slopes and padding proved exact before asking. Existing fixed search family, shared innovations, and existing reset expression.',
           helper_sources={name:(ROOT/'racer'/name).read_text()for name in ['ets.py','ets_optimizer.py']}))
    c,s=service()
    for name,function,variants in [('parameters','FORECAST.ETS.STAT',list(range(1,8))),('forecast','FORECAST.ETS',[0.,1.,2.,5.]),('confidence','FORECAST.ETS.CONFINT',[(1.,.95),(2.,.95)])]:
        save('step3-short-proof-'+name,ask_series(c,s,cases,function,variants))
    failures=[]
    for name,key,variants in [('parameters','stat',list(range(1,8))),('forecast','forecasts',[0.,1.,2.,5.])]:
        document=read(OUT/('step3-short-proof-'+name+'.json'));expected=[v for p in predictions for v in p[key]]
        for i,(p,a)in enumerate(zip(expected,document['answers'])):
            if p!=a.get('bits'):failures.append(dict(kind=name,case=i//len(variants),variant=variants[i%len(variants)],prediction=p,answer=a))
    save('step3-short-proof-result',dict(cases=len(cases),selections=len(cases)*11,failures=failures))
    print('short fresh proof',len(cases)*11,'failures',len(failures),failures[:8],flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('step', choices=['reanalysis', 'normal-map', 'normal-refinement','confidence-reanalysis','slope-map','slope-analysis','slope-impulses','impulse-analysis','slope-stores','short-proof'])
    args=parser.parse_args()
    {'reanalysis':reanalysis, 'normal-map':quantile_mapping, 'normal-refinement':quantile_refinement,
     'confidence-reanalysis':confidence_reanalysis,'slope-map':slope_mapping,'slope-analysis':slope_analysis,
     'slope-impulses':slope_impulses,'impulse-analysis':impulse_analysis,'slope-stores':slope_stores,'short-proof':short_proof}[args.step]()
