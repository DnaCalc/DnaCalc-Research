"""Round40: executable retained-answer analysis and frozen, timed questions.

Original MIT research from Excel answers, IEEE754 and the existing mathematical
tables. No restricted implementation sources. R39 files remain immutable.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import random
import time

from . import ets
from .ets import bits, value, read, grouped, sequential_total, initial_trend
from .fingerprint import AuditedClientMixin, CachedService, load_service, sha
from .question_designer import Census

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('/home/ubuntu/codex-runs/round40')
RETAINED = ROOT / 'evidence/round39/retained'
_LOADED_SOURCE = Path(__file__).read_text()


def save(name, document):
    path = OUT / (name + '.json')
    with path.open('x') as stream:
        json.dump(document, stream, separators=(',', ':'), allow_nan=False)
    return path


def frozen(name, document):
    document['producer_source'] = _LOADED_SOURCE
    document['source_sha256'] = hashlib.sha256(_LOADED_SOURCE.encode()).hexdigest()
    return save(name + '-frozen', document)


def service(offline=False):
    client, _ = load_service(ets.CHECKOUT)

    class TimedClient(AuditedClientMixin,client.Client):
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
    per_output=set()
    for path in sorted(RETAINED.glob('step*.json.gz')):
        document = read(path)
        if not isinstance(document,dict)or document.get('function')!='FORECAST.ETS.STAT'or not all(k in document for k in ['cases', 'labels', 'answers']): continue
        for i, row in enumerate(grouped(document)):
            for k in [1, 2, 3]:
                a = row.get(str(k))
                if a and a.get('type') == 'number':
                    census.observe([path.name, i, k], a['bits']);per_output.add((k,a['bits']))
    report = census.report()
    decimals = []
    for raw in report['counts']:
        x = value(json.loads(raw))
        decimals.append(dict(bits=bits(x), value=x, six_significant=float(format(x, '.6g'))==x,
                             nonseasonal_gamma=x==2.**-52))
    save('step1-full-reanalysis', dict(rounded_grid=rounded, representation_changes=representation,
         original=retained('step7-gamma-generalization-refutations')['matched'],
         rounded_matched=sum(row['matched'] for row in results), total=len(results), rows=results,
         census=report, per_output_patterns=sorted(per_output), parameters=decimals,
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


def ask_series(c, s, cases, function, variants, chunk_size=400, reuse_batches=False):
    """Existing ETS question construction with a throughput-sized batch."""
    rows=[];labels=[]
    for i,z in enumerate(cases):
        timeline=z.get('timeline',[float(j+1)for j in range(len(z['y']))])
        y=c.number_array([[v]for v in z['y']]);t=c.number_array([[v]for v in timeline]);season=c.number(float(z.get('season',0)))
        for v in variants(z)if callable(variants)else variants:
            if function=='FORECAST.ETS.STAT':args=[y,t,c.number(float(v)),season]
            elif function=='FORECAST.ETS':args=[c.number(z.get('targets',{}).get(str(v),timeline[-1]+v)),y,t,season]
            elif function=='FORECAST.ETS.CONFINT':
                h,conf=v;args=[c.number(z.get('targets',{}).get(str(h),timeline[-1]+h)),y,t,c.number(conf),season]
            else:args=[y,t]
            rows.append(args);labels.append([i,v])
    if reuse_batches:
        import hashlib
        answers=[]
        for start in range(0,len(rows),chunk_size):
            batch=rows[start:start+chunk_size]
            # Same canonical request identity as CachedService. Keep every
            # successful original batch; only unreturned work is made smaller.
            request=dict(function=function,rows=batch)
            identity=hashlib.sha256(json.dumps(request,sort_keys=True,separators=(',',':')).encode()).hexdigest()
            cached=(s.directory/(identity+'.json.gz')).exists()
            answers.extend(s.ask_rows(function,batch,chunk_size=chunk_size if cached else 250))
    else:answers=s.ask_rows(function,rows,chunk_size=chunk_size)
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


def parameter_census(resume=False):
    from .ets_optimizer import parameters
    rng=random.Random(400401);cases=[];predictions=[];raw=[]
    for i in range(10000):
        family=i%5;n=rng.randrange(8,33);exponent=rng.randrange(-200,201)
        if family==0:y=[rng.randrange(-4096,4097)/32 for _ in range(n)]
        elif family==1:
            drift=rng.randrange(-256,257)/16;curve=rng.randrange(-16,17)/128
            y=[drift*j+curve*j*j+rng.randrange(-8,9)/64 for j in range(n)]
        elif family==2:
            at=rng.randrange(n);size=rng.randrange(-256,257)/8
            y=[128.*j+(size if j>=at else 0.)for j in range(n)]
        elif family==3:
            period=rng.randrange(2,9);pattern=[rng.randrange(-128,129)/16 for _ in range(period)]
            drift=rng.randrange(-16,17)/16
            y=[pattern[j%period]+drift*j+rng.randrange(-4,5)/64 for j in range(n)]
        else:
            center=int(bits(rng.choice([1., math.pi,1/3])),16)
            y=[value(f'{center+rng.randrange(-32,33):016x}')for _ in range(n)]
        y=[math.ldexp(v,exponent)for v in y]
        case=dict(y=y,season=0,family=family)
        a,b=parameters(y);_,_,_,_,reset=fit_with_trend(y,a,b,initial_trend(y))
        predictions.append([bits(a),bits(.001 if reset else b),bits(2.**-52)])
        raw.append(dict(alpha=a,beta=b,reset=reset))
        cases.append(case)
    document=dict(cases=cases,predictions=predictions,raw_parameters=raw,
           families=['noise','curvature','step','periodic','adjacent-lattice'],
           method='Seed and all existing finite-search predictions frozen before asking. Census maps observed triples; no count selects a new point.',
           helper_source=(ROOT/'racer/ets_optimizer.py').read_text())
    if resume:
        original=read(OUT/'step4-census-frozen.json')
        if original['cases']!=cases or original['predictions']!=predictions:raise ValueError('census continuation changed questions or predictions')
        save('step4-census-continuation',dict(frozen_sha256=sha(OUT/'step4-census-frozen.json'),
             reason='One original 1000-selection service job exceeded 75 seconds. Successful batches are reused; remaining frozen selections use batches of 250.',producer_source=_LOADED_SOURCE))
    else:frozen('step4-census',document)
    c,s=service();document=ask_series(c,s,cases,'FORECAST.ETS.STAT',[1,2,3],chunk_size=1000,reuse_batches=resume)
    save('step4-census-answers',document)
    census=Census();failures=[];by_family=Counter();blocks=[]
    for i,(z,stat,pred)in enumerate(zip(cases,grouped(document),predictions)):
        actual=[stat[str(k)].get('bits',stat[str(k)])for k in[1,2,3]]
        census.observe(i,actual)
        if pred!=actual:
            failures.append(dict(case=i,family=z['family'],prediction=pred,answer=actual,raw=raw[i]))
            by_family[z['family']]+=1
        if (i+1)%500==0:blocks.append(dict(observations=i+1,distinct=len(census.counts)))
    save('step4-census-result',dict(census=census.report(),blocks=blocks,failures=failures,
         failures_by_family=dict(by_family),histories=len(cases),selections=len(cases)*3,
         limits='Near-constant initial-state rounding remains unresolved. A finite output map is not proof of global saturation or of a search algorithm. All parameter losses remain failures.'))
    print('parameter census',len(cases),'distinct',len(census.counts),'failed histories',len(failures),'families',dict(by_family),'blocks',blocks,flush=True)


def final_add_proof():
    import numpy as np
    from fractions import Fraction as Q
    from checker.exact import EXTENDED64 as E,BINARY64 as B
    from .ets_optimizer import parameters
    from .question_designer import choose_split
    rng=random.Random(400361);cases=[];attempts=0
    while len(cases)<24:
        attempts+=1
        if attempts>200000:raise ValueError('final-add search exhausted')
        n=3 if attempts%2 else 7
        y=[rng.randrange(-65536,65537)/32 for _ in range(n)]
        if not exact_short(y):continue
        a,b=parameters(padded_history(y));t=initial_trend(padded_history(y))
        l,tr,se,e,reset=fit_with_trend(y,a,b,t)
        if reset:l=y[-1]-se[(n-1)%2];tr=t
        center=(se[0]+se[1])/2;l+=center;se=[v-center for v in se]
        for h in [0.,1.,2.,5.,17.]:
            subtotal=l+h*tr;season=se[(n-1+int(h))%2]
            strict=bits(subtotal+season);pc64=bits(float(np.longdouble(subtotal)+np.longdouble(season)))
            if strict==pc64:continue
            certified=bits(float(B.round(E.round(Q(subtotal)+Q(season)))))
            if certified!=pc64:raise ValueError('native and rational addition disagree')
            cases.append(dict(y=y,season=0,horizon=h,strict=strict,pc64=pc64,
                              subtotal=bits(subtotal),seasonal_value=bits(season),alpha=a,beta=b))
            break
    selection=choose_split({name:[z[name]for z in cases]for name in ['strict','pc64']})
    frozen('step3-final-add',dict(cases=cases,predictions={name:[z[name]for z in cases]for name in ['strict','pc64']},
           first_split=selection,attempts=attempts,
           method='Fixed final seasonal addition. Mine its binary64/PC64 rounding-midpoint disagreements; independently certify every proposed answer with exact rationals. Existing parameter search is frozen.'))
    c,s=service();rows=[]
    for z in cases:
        rows.append([c.number(len(z['y'])+z['horizon']),c.number_array([[v]for v in z['y']]),c.number_array([[float(j+1)]for j in range(len(z['y']))]),c.number(0.)])
    answers=s.ask_rows('FORECAST.ETS',rows,chunk_size=1000)
    save('step3-final-add-proof',dict(cases=cases,answers=answers,receipts=s.receipts,
         counts={name:sum(z[name]==a.get('bits')for z,a in zip(cases,answers))for name in ['strict','pc64']}))
    print('fresh final-add', {name:sum(z[name]==a.get('bits')for z,a in zip(cases,answers))for name in ['strict','pc64']},flush=True)


def final_product_proof():
    import numpy as np
    from fractions import Fraction as Q
    from checker.exact import EXTENDED64 as E,BINARY64 as B
    from .ets_optimizer import parameters
    from .ets_confidence import variance
    from .question_designer import product_residue
    z=value(retained('step6-quantile-coupling-normal')['answers'][5]['bits'])
    rng=random.Random(400371);cases=[];attempts=0
    while len(cases)<24:
        attempts+=1
        if attempts>200000:raise ValueError('confidence-product search exhausted')
        y=[rng.randrange(-65536,65537)/32 for _ in range(3)]
        a,b=parameters(padded_history(y));t=initial_trend(padded_history(y))
        _,_,_,errors,reset=fit_with_trend(y,a,b,t);pb=.001 if reset else b
        sigma=sequential_total(e*e for e in errors)/len(y)
        for h in [1.,2.,5.,17.]:
            root=math.sqrt(sigma*variance(a,pb,h,39))
            product,residue=product_residue(z,root)
            strict=bits(product);pc64=bits(float(np.longdouble(z)*np.longdouble(root)))
            if strict==pc64:continue
            if bits(float(B.round(E.round(Q(z)*Q(root)))))!=pc64:raise ValueError('product replay differs')
            cases.append(dict(y=y,season=0,horizon=h,strict=strict,pc64=pc64,z=bits(z),root=bits(root),residue=str(residue)))
            break
    frozen('step3-final-product',dict(cases=cases,predictions={name:[z[name]for z in cases]for name in['strict','pc64']},
           attempts=attempts,method='Fixed final confidence product; exact rational midpoint residues select fresh distinguishing histories. Independent retained normal quantile at p=.975.'))
    c,s=service();rows=[]
    for cs in cases:
        rows.append([c.number(3+cs['horizon']),c.number_array([[v]for v in cs['y']]),c.number_array([[1.],[2.],[3.]]),c.number(.95),c.number(0.)])
    answers=s.ask_rows('FORECAST.ETS.CONFINT',rows,chunk_size=1000)
    save('step3-final-product-proof',dict(cases=cases,answers=answers,receipts=s.receipts,
         counts={name:sum(cs[name]==a.get('bits')for cs,a in zip(cases,answers))for name in['strict','pc64']}))
    print('fresh confidence product',{name:sum(cs[name]==a.get('bits')for cs,a in zip(cases,answers))for name in['strict','pc64']},flush=True)


def search_structures():
    """Refute fixed familiar search programs; never append observed points."""
    from .ets_optimizer import objective
    from .question_designer import choose_split
    document=read(OUT/'step4-census-answers.json');stat=grouped(document)
    ids=[i for i,z in enumerate(document['cases'])if z['family']<4][:80]
    # Complete raw axes from the owner-supplied structural hypothesis. All
    # feasible Cartesian pairs are scored, rather than extending the sparse
    # observed list by whichever extra point happened to be missing.
    axes=[0.,1/12,.1,.125,1/6,.25,1/3,.5,2/3,.75,5/6,.9,1.]
    pool=sorted(set((float(format(min(.998,max(.002,a)),'.6g')),
                     float(format(min(max(.001,b),min(.998,max(.002,a))-.001),'.6g')))
                    for a in axes for b in axes))
    grid_failures=[]
    for i in ids:
        y=document['cases'][i]['y'];a,b=min(pool,key=lambda pair:objective(y,*pair))
        _,_,_,_,reset=fit_with_trend(y,a,b,initial_trend(y));pb=.001 if reset else b
        actual=[stat[i][str(k)]['bits']for k in[1,2]]
        if [bits(a),bits(pb)]!=actual:grid_failures.append(dict(case=i,predicted=[a,pb],answer=actual))

    def project(v):
        a=min(.998,max(.002,v[0]));b=min(a-.001,max(.001,v[1]));g=min(1-a-.001,max(.001,v[2]))
        return[a,b,g]
    def nm(y,step,iterations,stored):
        # Nonseasonal gamma is the independently measured epsilon. Keeping
        # three coordinates tests the centroid-divide-by-three hypothesis.
        def evaluate(v):return objective(y,*project(v)[:2])
        def move(v):return project(v)if stored else list(v)
        start=[.9,.001,.001]
        points=[start]+[[v+(step if j==k else 0.)for j,v in enumerate(start)]for k in range(3)]
        points=[move(v)for v in points]
        scored=[(evaluate(v),v)for v in points]
        for _ in range(iterations):
            scored.sort(key=lambda p:p[0]);worst=scored[-1][1]
            center=[sum(p[1][j]for p in scored[:-1])/3 for j in range(3)]
            reflected=move([2*c-w for c,w in zip(center,worst)]);fr=evaluate(reflected)
            if fr<scored[0][0]:
                expanded=move([c+2*(r-c)for c,r in zip(center,reflected)]);fe=evaluate(expanded)
                scored[-1]=(fe,expanded)if fe<fr else(fr,reflected)
            elif fr<scored[-2][0]:scored[-1]=(fr,reflected)
            else:
                target=reflected if fr<scored[-1][0]else worst
                contracted=move([(c+t)/2 for c,t in zip(center,target)]);fc=evaluate(contracted)
                if fc<min(fr,scored[-1][0]):scored[-1]=(fc,contracted)
                else:
                    best=scored[0][1]
                    scored=[scored[0]]+[(evaluate(v),v)for v in[move([(a+b)/2 for a,b in zip(best,p)])for _,p in scored[1:]]]
        a,b,_=project(min(scored,key=lambda p:p[0])[1]);return float(format(a,'.6g')),float(format(b,'.6g'))
    programs=[]
    for step in [.1,.5,1.]:
        for stored in[False,True]:
            for iterations in[24,64]:
                failures=[]
                for i in ids[:24]:
                    y=document['cases'][i]['y'];scale=math.ldexp(1.,math.frexp(max(abs(v)for v in y))[1]);normal=[v/scale for v in y]
                    a,b=nm(normal,step,iterations,stored)
                    _,_,_,_,reset=fit_with_trend(y,a,b,initial_trend(y));pb=.001 if reset else b
                    actual=[stat[i][str(k)]['bits']for k in[1,2]]
                    if [bits(a),bits(pb)]!=actual:failures.append(dict(case=i,predicted=[a,pb],answer=actual))
                programs.append(dict(step=step,stored_projection=stored,iterations=iterations,failures=failures))
    save('step4-search-structures',dict(cartesian_axes=axes,cartesian_pairs=pool,
         cartesian_cases=ids,cartesian_failures=grid_failures,simplex_programs=programs,
         limits='These are finite program refutations, not a proof against every Nelder-Mead seed or termination rule. The sparse original family is not expanded; unidentified search behavior remains open.'))
    print('Cartesian grid refutes',len(grid_failures),'/',len(ids),'simplex refutations',[len(p['failures'])for p in programs],flush=True)


def reset_boundaries():
    """Candidate-predicted crossings, then balanced frozen boundary questions."""
    import itertools
    import numpy as np
    from .ets_optimizer import parameters
    from .question_designer import Boundary,vector_path,choose_split
    design=read(OUT/'step4-census-frozen.json')['cases']
    bases=[]
    for n in[8,12,32]:
        for z in design:
            if len(z['y'])==n and z['family']==0 and parameters(z['y'])[0]==.002:
                exponent=math.frexp(max(abs(v)for v in z['y']))[1]
                bases.append([math.ldexp(v,-exponent)for v in z['y']]);break
    bases.append([-2.,8.,-13.,10.,-14.,12.,-17.,-11.,-29.,27.,30.,27.,-29.,-11.,-17.,12.,-14.,10.,-13.,8.,-2.])
    bases += [[10.+j*j/16 for j in range(n)]for n in[8,12]]
    profiles=list(itertools.product(['current','states','stat','full','transitions'],['population','sample'],['window','count']))
    def predicate(y,a,b,profile,center=False):
        mode,denom,multiplier=profile;n=len(y);initial=initial_trend(y);l=y[0];t=initial;se=[0.,0.];trends=[initial]
        for j,v in enumerate(y[1:],1):
            q=j%2;e=v-((l+t)+se[q]);l+=t+a*e;t+=b*e;se[q]+=2.**-52*e;trends.append(t)
        w=n if n<=10 else(n+1)//2
        if mode=='current':vs=trends[-(w-1):]
        elif mode=='states':vs=trends[-w:]
        elif mode=='stat':vs=trends[-(n if n<=10 else min(int(.4*n),10)):]
        elif mode=='full':vs=trends
        else:vs=trends[1:]
        count=len(vs);mean=sequential_total(vs)/count
        variance=sequential_total((v-mean)**2 for v in vs)/(count if denom=='population'else count-1)
        factor=w/(w-1)if multiplier=='window'else(count+1)/count
        return -mean if center else abs(mean)<math.sqrt(variance)*factor
    cases=[];decision_predictions={str(p):[]for p in profiles};traces=[]
    last=1<<52
    for base_id,base in enumerate(bases):
        a,b=parameters(base);n=len(base)
        span=max(1.,2*max(abs(v)for v in base)/n)
        left=[v-span*j for j,v in enumerate(base)];right=[v+span*j for j,v in enumerate(base)]
        def path(index):return vector_path(left,right,index,last)
        indices={0,last,last//2}
        for profile in profiles:
            # Locate both candidate crossings using a bounded offline scan;
            # then bisect each opposite-decision path to adjacent indices.
            inside=round((predicate(base,a,b,profile,True)/span+1)*last/2)
            grid=sorted(set([j*last//64 for j in range(65)]+[max(0,min(last,inside))]))
            flags=[predicate(path(j),a,b,profile)for j in grid]
            for lo,hi,lv,hv in zip(grid,grid[1:],flags,flags[1:]):
                if lv==hv:continue
                edge=Boundary(lo,hi,lv,hv)
                while (mid:=edge.question())is not None:edge.observe(mid,predicate(path(mid),a,b,profile))
                indices.update([edge.low,edge.high])
                traces.append(dict(base=base_id,profile=profile,low=edge.low,high=edge.high))
        ordered=sorted(indices);pool=sorted(set(ordered+[(a+b)//2 for a,b in zip(ordered,ordered[1:])]))
        predicted={str(p):[predicate(path(j),a,b,p)for j in pool]for p in profiles}
        groups=[list(predicted)];selected=set()
        for _ in range(16):
            options=[(len(group),group,choose_split(predicted,group))for group in groups]
            options=[x for x in options if x[2]is not None]
            if not options:break
            _,group,choice=max(options,key=lambda x:x[0]);index=choice['index'];selected.add(index)
            groups.remove(group)
            for flag in[False,True]:
                bucket=[name for name in group if predicted[name][index]==flag]
                if bucket:groups.append(bucket)
        for index in sorted(selected):
            j=pool[index];y=path(j)
            case=dict(y=y,season=0,base=base_id,path_index=j,raw_alpha=a,raw_beta=b)
            cases.append(case)
            for p in profiles:decision_predictions[str(p)].append(predicted[str(p)][index])
    frozen('step5-reset-boundaries',dict(cases=cases,predictions=decision_predictions,profiles=profiles,
           candidate_brackets=traces,method='Offline path bisection to adjacent dyadic indices; choose questions that split candidate decisions. No Excel answer chooses a window or coefficient.'))
    c,s=service();stats=ask_series(c,s,cases,'FORECAST.ETS.STAT',[1,2],chunk_size=250)
    forecasts=ask_series(c,s,cases,'FORECAST.ETS',[0.,1.],chunk_size=250)
    save('step5-reset-stat',stats);save('step5-reset-forecast',forecasts)
    decisions=[];unknown=[];changed=[]
    for i,(z,st,fo)in enumerate(zip(cases,grouped(stats),grouped(forecasts))):
        y=z['y'];a,b=z['raw_alpha'],z['raw_beta']
        if st['1']['bits']!=bits(a):changed.append(i);decisions.append(None);continue
        l,t,se,e,_=fit_with_trend(y,a,b,initial_trend(y));center=(se[0]+se[1])/2
        seasons=[v-center for v in se];states=[(l+center,t),(y[-1]-se[(len(y)-1)%2]+center,initial_trend(y))]
        choices=[]
        for decision,(level,trend)in enumerate(states):
            pred=[bits(float(np.longdouble(level+h*trend)+np.longdouble(seasons[(len(y)-1+int(h))%2])))for h in[0.,1.]]
            if pred==[fo[str(h)]['bits']for h in[0.,1.]]:choices.append(bool(decision))
        if len(choices)==1:decisions.append(choices[0])
        else:decisions.append(None);unknown.append(dict(case=i,choices=choices))
    survivors=[name for name,values in decision_predictions.items()if all(actual is None or pred==actual for pred,actual in zip(values,decisions))]
    save('step5-reset-result',dict(histories=len(cases),decisions=decisions,survivors=survivors,
         changed_search=changed,unidentified_state_arithmetic=unknown,
         limits='Candidate reset decisions are conditioned on the frozen raw search parameters and recurrence. Unknown state readbacks and changed search premises remain failures, not ignored passes.'))
    print('reset boundaries',len(cases),'survivors',survivors,'unknown',len(unknown),'changed search',len(changed),flush=True)


def seasonal_initialization(partial=False):
    from .ets_seasonal import moving_initial,fit
    import itertools
    import numpy as np
    rng=random.Random(400901 if partial else 400601);cases=[]
    name='step9-partial-cycles'if partial else'step6-seasonal-initialization'
    for p in range(2,9):
        lengths=sorted(set([2*p+1,2*p+2,3*p-1]if partial else[3*p,3*p+1,4*p]))
        pattern=[p*rng.randrange(-16,17)/4 for _ in range(p)]
        for n in lengths:
            for at in[-1,0,n//2,n-1]:
                y=[128.*j+pattern[j%p]for j in range(n)]
                if at>=0:y[at]=math.nextafter(y[at],math.inf)if at==n//2 else y[at]+1/16
                cases.append(dict(y=y,season=p,at=at,horizons=[float(h)for h in range(2*p+1)]))
    frozen(name+'-design',dict(cases=cases,predictions=None,
           method='Positive dyadic trend, distinct periodic phases, impulses and one-bit observation steps; census first, no optimizing parameter injected into the application.'))
    c,s=service();parameters=ask_series(c,s,cases,'FORECAST.ETS.STAT',[1,2,3],chunk_size=200)
    save(name+'-parameters',parameters);st=grouped(parameters)
    programs=list(itertools.product(['original','pad3p','pad3p1'],['direct','subtract_add','add_difference'],['centered','raw']))
    predicted={};initials={}
    for extension,order,adjustment in programs:
        key=str((extension,order,adjustment));out=[];starts=[]
        for cs,ss in zip(cases,st):
            y=cs['y'];p=cs['season'];init_y=y
            if extension!='original':
                minimum=3*p+int(extension=='pad3p1');trend=initial_trend(y)
                init_y=y+[y[j%p]+(j-j%p)*trend for j in range(len(y),minimum)]
            if len(init_y)<3*p:
                starts.append(None);out.extend([None]*len(cs['horizons']));continue
            init=moving_initial(init_y,p,order,'centered'if adjustment=='centered'else'raw')
            a,b,g=[value(ss[str(k)]['bits'])for k in[1,2,3]]
            l,t,se,errors,pr=fit(y,p,a,b,g,init)
            out.extend(bits(float(np.longdouble(l+h*t)+np.longdouble(se[(len(y)-1+int(h))%p])))for h in cs['horizons'])
            starts.append(init)
        predicted[key]=out;initials[key]=starts
    frozen(name,dict(cases=cases,predictions=predicted,parameters=parameters,initial_states=initials,
           helper_source=(ROOT/'racer/ets_seasonal.py').read_text(),
           limits='Conditional on reported STAT leaves and the existing recursion. Positive drift separates ordinary trend reset. Original and two structurally motivated phase extensions are fixed before forecast answers.'))
    forecasts=ask_series(c,s,cases,'FORECAST.ETS',lambda z:z['horizons'],chunk_size=200)
    save(name+'-answers',forecasts);actual=[a.get('bits')for a in forecasts['answers']];results={}
    for key,pred in predicted.items():
        failures=[dict(selection=i,prediction=p,answer=forecasts['answers'][i],case=forecasts['labels'][i][0])for i,p in enumerate(pred)if p!=actual[i]]
        results[key]=dict(matched=len(actual)-len(failures),failures=len(failures),first=failures[:12],by_period=dict(Counter(cases[z['case']]['season']for z in failures)))
    census=Census()
    for i,ss in enumerate(st):census.observe(i,[ss[str(k)].get('bits')for k in[1,2,3]])
    save(name+'-result',dict(selections=len(actual),programs=results,parameter_census=census.report(),
         limits='Every unavailable initializer and every mismatch fails. No seasonal search or initial-state body is installed from conditional leaves.'))
    print(name,len(actual),{k:v['failures']for k,v in results.items()},flush=True)


def seasonal_exact_proof(resume=False):
    from .ets_seasonal import initial,fit
    rng=random.Random(400641);cases=[];forecast=[];stat=[];ci=[]
    for i in range(168):
        p=2+i%7;cycles=rng.randrange(3,32//p+1);n=p*cycles
        exponent=rng.randrange(-300,301);offset=rng.randrange(-256,257)/8
        pattern=[p*rng.randrange(-32,33)/8 for _ in range(p)]
        trend=rng.randrange(-32,33)/16
        y=[math.ldexp(offset+pattern[j%p]+trend*j,exponent)for j in range(n)]
        l,t,se,errors,pr=fit(y,p,.751,.001,.001,initial(y,p,'firsttwo','first_detrend'))
        if any(e!=0. for e in errors):raise ValueError('exact seasonal design has an innovation')
        cases.append(dict(y=y,season=p,horizons=[0.,1.,float(p),float(p+1),float(2*p+5)]))
        forecast.extend(bits((l+h*t)+se[(n-1+int(h))%p])for h in cases[-1]['horizons'])
        stat.extend(bits(v)for v in[.751,.001,.001,0.,0.,0.,0.,1.])
        ci.extend([bits(0.)]*3)
    document=dict(cases=cases,predictions=dict(forecast=forecast,stat=stat,confidence=ci),
           method='Whole-cycle affine patterns, all periods 2..8, at least three cycles, varied signs/offsets/exponents. The shared recurrence calculates zero innovations; no observed parameter leaf is supplied.')
    if resume:
        original=read(OUT/'step6-exact-seasonal-frozen.json')
        if original['cases']!=cases or original['predictions']!=document['predictions']:raise ValueError('seasonal continuation changed predictions')
        save('step6-exact-seasonal-continuation',dict(original_frozen_sha256=sha(OUT/'step6-exact-seasonal-frozen.json'),producer_source=_LOADED_SOURCE,
             reason='The v2 relay reply lacked the legacy session receipt. Preserve it unqualified; repeat the same frozen questions through the audited touchstone socket.'))
    else:frozen('step6-exact-seasonal',document)
    c,s=service();results={}
    for name,function,variants,pred in [('forecast','FORECAST.ETS',lambda z:z['horizons'],forecast),('stat','FORECAST.ETS.STAT',list(range(1,9)),stat),('confidence','FORECAST.ETS.CONFINT',[(1.,.25),(2.,.95),(5.,math.nextafter(1.,0.))],ci)]:
        d=ask_series(c,s,cases,function,variants,chunk_size=200);save('step6-exact-seasonal-'+name,d)
        failures=[dict(selection=i,prediction=p,answer=a,case=d['labels'][i][0])for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
        results[name]=dict(selections=len(pred),failures=failures)
    save('step6-exact-seasonal-result',dict(results=results))
    print('fresh exact seasonal',{k:(v['selections'],len(v['failures']))for k,v in results.items()},flush=True)


def detection_proof():
    import itertools
    from .ets_seasonal import detection_candidate
    rng=random.Random(400801);cases=[]
    for p in range(2,13):
        for n in range(2*p-1,2*p+6):
            for kind in['spike','cosine']:
                y=[float(p-1 if j%p==0 else-1)for j in range(n)]if kind=='spike'else[math.cos(2*math.pi*j/p)for j in range(n)]
                cases.append(dict(y=y,kind=kind,period=p,n=n))
    noise=[rng.randrange(-128,129)/16 for _ in range(32)]
    for exponent in range(-16,17):
        amp=2.**exponent
        cases.append(dict(y=[v+amp*(1 if j%2==0 else-1)for j,v in enumerate(noise)],kind='noise_threshold',amplitude=amp))
    for exponent in range(-8,9):
        a=2.**(exponent/4)
        cases.append(dict(y=[a*[1.,0.,-1.,0.][j%4]+[1.,1.,1.,1.,-1.,-1.,-1.,-1.][j%8]for j in range(32)],kind='competing_peaks',amplitude=a))
    programs=list(itertools.product(['mean','linear','difference'],['biased','unbiased','pearson'],['first','global'],[1,2,3],['zero','two','normal95']))
    predictions={}
    for pre,norm,search,offset,threshold in programs:
        pred=[]
        for z in cases:
            cutoff=0. if threshold=='zero'else(2. if threshold=='two'else 1.96)/math.sqrt(len(z['y']))
            pred.append(detection_candidate(z['y'],pre,norm,search,offset,cutoff)[0])
        predictions[str((pre,norm,search,offset,threshold))]=pred
    frozen('step8-detection',dict(cases=cases,predictions=predictions,
           helper_source=(ROOT/'racer/ets_seasonal.py').read_text(),
           method='Full short-length census around 2p+1..2p+3, two distinct periodic shapes, noise exponent edges, and competing local peaks. Fixed correlation normalizations and cutoffs.'))
    c,s=service();d=ask_series(c,s,cases,'FORECAST.ETS.SEASONALITY',[0],chunk_size=200)
    save('step8-detection-answers',d);actual=[value(a['bits'])if a.get('type')=='number'else a for a in d['answers']]
    results={key:dict(failures=sum(p!=a for p,a in zip(pred,actual)),first=[dict(case=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,actual))if p!=a][:8])for key,pred in predictions.items()}
    bounds=[]
    for p in range(2,13):
        for kind in['spike','cosine']:
            rows=[dict(n=z['n'],answer=a)for z,a in zip(cases,actual)if z.get('period')==p and z['kind']==kind]
            bounds.append(dict(period=p,kind=kind,observations=rows))
    save('step8-detection-result',dict(selections=len(cases),programs=results,boundaries=bounds,
         survivors=[key for key,v in results.items()if not v['failures']],
         limits='A finite length/noise map does not prove a universal lag bound. Refuted correlation programs do not identify a detector or justify expanding automatic forecasts.'))
    print('detector fresh',len(cases),'programs',len(programs),'survivors',[key for key,v in results.items()if not v['failures']],flush=True)


def calendar_proof():
    import calendar
    import datetime as dt
    origin=dt.date(1899,12,30)
    serial=lambda date:float((date-origin).days)
    cases=[];predictions=[]
    for kind in['first','end','middle','mixed','quarter','year','day28','day29','day30','day31']:
        for n in[8,12,24]:
            start=dt.date(2019,9,1)
            def date_at(j):
                if kind.startswith('day'):return start+dt.timedelta(days=int(kind[3:])*j)
                month=start.year*12+start.month-1+j*(3 if kind=='quarter'else 12 if kind=='year'else 1)
                year,month0=divmod(month,12);month=month0+1
                day=calendar.monthrange(year,month)[1]if kind=='end'else 15 if kind=='middle'else(1 if j%2==0 else 15)if kind=='mixed'else 1
                return dt.date(year,month,day)
            timeline=[serial(date_at(j))for j in range(n)]
            targets={str(h):serial(date_at(n-1+int(h)))for h in[0.,1.,2.,5.]}
            y=[2.+128.*j for j in range(n)]
            cases.append(dict(y=y,season=0,kind=kind,timeline=timeline,targets=targets))
            predictions.extend(bits(2.+128.*(n-1+h))for h in[0.,1.,2.,5.])
    frozen('step10-calendar',dict(cases=cases,predictions=predictions,
           method='Calendar-index affine model: first/end/middle/mixed month days, quarters, years and 28..31-day numeric controls across leap dates. Original numeric serials, no workbook settings are changed.'))
    c,s=service();stats=ask_series(c,s,cases,'FORECAST.ETS.STAT',[8],chunk_size=200)
    forecasts=ask_series(c,s,cases,'FORECAST.ETS',[0.,1.,2.,5.],chunk_size=200)
    save('step10-calendar-stat',stats);save('step10-calendar-answers',forecasts)
    failures=[dict(selection=i,prediction=p,answer=a,case=forecasts['labels'][i][0])for i,(p,a)in enumerate(zip(predictions,forecasts['answers']))if p!=a.get('bits')]
    save('step10-calendar-result',dict(selections=len(predictions),failures=failures,
         steps=[dict(kind=z['kind'],n=len(z['y']),answer=a)for z,a in zip(cases,stats['answers'])],
         limits='Calendar index is a structural candidate. Every changed numeric prediction and every typed error remains a failure. Date1904 and mixed-day inference need independent evidence; no calendar preparation body is installed.'))
    print('calendar fresh',len(predictions),'failures',len(failures),flush=True)


def calendar_followup():
    """Separate calendar units and residual-day horizon trees on new dates."""
    import calendar
    import datetime as dt
    origin=dt.date(1899,12,30)
    serial=lambda date:float((date-origin).days)
    def month_at(start,months,day):
        year,month=divmod(start.year*12+start.month-1+months,12);month+=1
        return dt.date(year,month,min(day,calendar.monthrange(year,month)[1]))
    cases=[];predictions={name:[]for name in ['calendar_index','month_plus_fixed31','month_plus_last_length','numeric_serial']}
    for start in [dt.date(2000,1,1),dt.date(2020,2,1),dt.date(2024,11,1)]:
        for step in [1,2,3,4,5,6,7,8,9,10,11,12,18,24,36]:
            for day in [1,31]:
                n=8;dates=[month_at(start,j*step,day)for j in range(n)];last=dates[-1]
                target_dates={str(float(h)):month_at(start,(n-1+h)*step,day)for h in [1,2,5]}
                z=dict(y=[2.+128.*j for j in range(n)],season=0,month_step=step,start=str(start),day=day,
                       timeline=[serial(date)for date in dates],targets={key:serial(date)for key,date in target_dates.items()})
                cases.append(z)
                for h in [1,2,5]:
                    date=target_dates[str(float(h))];delta=date.day-last.day
                    predictions['calendar_index'].append(bits(z['y'][-1]+128.*h))
                    predictions['month_plus_fixed31'].append(bits(z['y'][-1]+128.*(h+delta/(31.*step))))
                    predictions['month_plus_last_length'].append(bits(z['y'][-1]+128.*(h+delta/(float(calendar.monthrange(last.year,last.month)[1])*step))))
                    predictions['numeric_serial'].append(bits(z['y'][-1]+128.*(serial(date)-serial(last))/(z['timeline'][-1]-z['timeline'][-2])))
    units={name:[]for name in ['jan2001','jan2000','dec1899','average_month','thirty_day']}
    for z in cases:
        step=z['month_step']
        for name,anchor in [('jan2001',dt.date(2001,1,1)),('jan2000',dt.date(2000,1,1)),('dec1899',dt.date(1899,12,30))]:
            units[name].append(bits(float((month_at(anchor,step,anchor.day)-anchor).days)))
        units['average_month'].append(bits(365.25*step/12.));units['thirty_day'].append(bits(30.*step))
    frozen('step10-calendar-followup',dict(cases=cases,predictions=predictions,units=units,
           method='Fresh start years, leap boundaries, two days and fifteen whole-month intervals. Fixed Gregorian anchor rules and day-residual horizon trees; no fitted date-dependent coefficients.'))
    c,s=service();stats=ask_series(c,s,cases,'FORECAST.ETS.STAT',[8],chunk_size=200)
    forecasts=ask_series(c,s,cases,'FORECAST.ETS',[1.,2.,5.],chunk_size=200)
    save('step10-calendar-followup-stat',stats);save('step10-calendar-followup-answers',forecasts)
    def assess(pred,actual):
        losses=[dict(selection=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,actual))if p!=a.get('bits')]
        return dict(selections=len(actual),failures=len(losses),first=losses[:12])
    save('step10-calendar-followup-result',dict(units={key:assess(pred,stats['answers'])for key,pred in units.items()},
         forecasts={key:assess(pred,forecasts['answers'])for key,pred in predictions.items()},
         limits='Numeric-unit and forecast hypotheses are independent. A surviving unit rule alone does not identify timeline preparation, aggregation or forecast arithmetic. Every mismatch and typed error fails.'))
    print('calendar followup',len(cases),{key:assess(pred,forecasts['answers'])['failures']for key,pred in predictions.items()},flush=True)


def detection_followup():
    """Fresh phase/length questions for the newly observed 2n/5 lag envelope."""
    import itertools
    from .ets_seasonal import detection_candidate
    cases=[]
    for p in range(2,13):
        for n in sorted(set([max(5,(5*p)//2-1),(5*p+1)//2,(5*p+1)//2+2,3*p+1])):
            for kind in ['shifted_spike','shifted_cosine']:
                y=[float(p-1 if (j+1)%p==0 else -1)for j in range(n)]if kind=='shifted_spike'else[math.cos(2*math.pi*(j+1)/p)for j in range(n)]
                cases.append(dict(y=y,period=p,n=n,kind=kind))
    predictions={}
    for pre,norm,search,threshold in itertools.product(['mean','linear','difference'],['biased','unbiased','pearson'],['first','global'],['zero','two','normal95']):
        key=str(('acf_2n5',pre,norm,search,threshold));predictions[key]=[]
        for z in cases:
            n=len(z['y']);cutoff=0. if threshold=='zero'else(2. if threshold=='two'else 1.96)/math.sqrt(n)
            predictions[key].append(detection_candidate(z['y'],pre,norm,search,n-2*((2*n)//5),cutoff)[0])
    for pre,padding,rounding in itertools.product(['mean','linear','difference'],[1,4],['nearest','floor','ceil']):
        values=[]
        for z in cases:
            y=z['y'];n=len(y);mean=sequential_total(y)/n;t=initial_trend(y)
            x=[v-mean for v in y]if pre=='mean'else[v-(mean+(j-(n-1)/2)*t)for j,v in enumerate(y)]if pre=='linear'else[v-u for u,v in zip(y,y[1:])]
            size=padding*n;low=max(1,math.ceil(size/max(2,(2*n)//5)));high=size//2
            spectrum=[(sequential_total(v*math.cos(2*math.pi*k*j/size)for j,v in enumerate(x))**2+
                       sequential_total(v*math.sin(2*math.pi*k*j/size)for j,v in enumerate(x))**2,k)for k in range(low,high+1)]
            if not spectrum:
                values.append(0);continue
            _,frequency=max(spectrum,key=lambda q:q[0]);period=size/frequency
            out=round(period)if rounding=='nearest'else math.floor(period)if rounding=='floor'else math.ceil(period)
            values.append(out)
        predictions[str(('periodogram',pre,padding,rounding))]=values
    frozen('step8-detection-followup',dict(cases=cases,predictions=predictions,
           helper_source=(ROOT/'racer/ets_seasonal.py').read_text(),
           method='Fresh phase-shifted patterns at the new length envelope. Fixed correlation and direct periodogram families; no fitted cutoff or extra peak correction.'))
    c,s=service();d=ask_series(c,s,cases,'FORECAST.ETS.SEASONALITY',[0],chunk_size=200)
    save('step8-detection-followup-answers',d)
    actual=[value(a['bits'])if a.get('type')=='number'else a for a in d['answers']]
    programs={}
    for key,pred in predictions.items():
        failures=[dict(case=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,actual))if p!=a]
        programs[key]=dict(failures=len(failures),first=failures[:8])
    save('step8-detection-followup-result',dict(selections=len(cases),programs=programs,
         survivors=[key for key,p in programs.items()if not p['failures']],
         limits='The observed 2n/5 envelope is not an implementation. These two complete families are stopped on refutation; no detector scope is expanded.'))
    print('detector followup',len(cases),'survivors',[key for key,p in programs.items()if not p['failures']],flush=True)


def audit_evidence():
    """Revalidate retained transport evidence and list every slow live request."""
    from .fingerprint import accept
    import csv
    counts=Counter();types=Counter();functions=Counter();builds=Counter();by_id={};raw=[];refusals=[];rejected=[]
    for path in sorted((OUT/'service-cache').glob('*.json.gz')):
        document=read(path);request,reply=document['request'],document['reply']
        if not reply.get('ok'):
            refusals.append(dict(file=path.name,id=reply.get('id'),error=reply.get('error'),calls=len(request['rows'])));continue
        accept(reply,request['function'],request['rows'],audit=document.get('audit'),result_index=request.get('result_index'))
        transport='touchstone/1'if document.get('audit')else'session-host-v1'
        counts[transport]+=1;functions[request['function']]+=len(reply['rows']);builds[str(reply['excel_build'])]+=1
        types.update(a['type']for a in reply['rows']);by_id[reply['id']]=path.name
        raw.append(dict(file=path.name,sha256=sha(path),id=reply['id'],transport=transport,calls=len(reply['rows'])))
    for path in sorted((OUT/'unqualified-service-cache').glob('*.json.gz')):
        document=read(path);request,reply=document['request'],document['reply']
        try:accept(reply,request['function'],request['rows'],audit=document.get('audit'),result_index=request.get('result_index'))
        except ValueError as error:rejected.append(dict(file=path.name,sha256=sha(path),id=reply['id'],calls=len(reply['rows']),reason=str(error)))
        else:raise ValueError('unqualified acquisition unexpectedly accepted')
        by_id[reply['id']]='unqualified-service-cache/'+path.name
    stages={}
    for path in sorted(OUT.glob('*.json')):
        document=read(path)
        if not isinstance(document,dict):continue
        for receipt in document.get('receipts',[]):
            if 'file'in receipt:stages.setdefault(receipt['file'],[]).append(path.name)
    waits=[json.loads(line)for line in(OUT/'service-times.jsonl').read_text().splitlines()]
    with (OUT/'service-times-final.jsonl').open('x')as stream:
        for record in waits:stream.write(json.dumps(record)+'\n')
    slow=[]
    for index,record in enumerate(waits,1):
        if record['seconds']<=3.:continue
        name=by_id.get(record.get('id'),'')
        row=dict(sequence=index,function=record['function'],result_index=record.get('result_index'),calls=record['calls'],seconds=record['seconds'],
                 request_id=record.get('id',''),question_cache=name,stages=';'.join(stages.get(name,[])),error=record.get('detail',''))
        if not name and not row['error']:raise ValueError('slow successful request has no retained reply')
        slow.append(row)
    with (OUT/'excel-waits.tsv').open('x',newline='')as stream:
        writer=csv.DictWriter(stream,fieldnames=['sequence','function','result_index','calls','seconds','request_id','question_cache','stages','error'],delimiter='\t')
        writer.writeheader();writer.writerows(slow)
    sources={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),ROOT/'racer/fingerprint.py',ROOT/'racer/question_designer.py',ROOT/'racer/ets_proof.py']}
    save('final-service-audit',dict(successful_replies=len(raw),selections=sum(functions.values()),transports=counts,functions=functions,
         answer_types=types,excel_builds=builds,raw=raw,service_refusals=refusals,unqualified_rejected=rejected,
         timed_live_requests=len(waits),slow_requests=len(slow),slow_seconds=sum(row['seconds']for row in slow),
         longest_seconds=max(row['seconds']for row in slow),waits_sha256=sha(OUT/'excel-waits.tsv'),sources=sources,
         limits='All qualified normalized replies are revalidated against unchanged typed input echoes and worker provenance. Client v2 frame/environment hashes are distinct from server receipts. Service errors and unqualified replies never enter numeric proof counts. Every timed live request over three seconds appears in excel-waits.tsv.'))
    print('service audit',len(raw),sum(functions.values()),'slow waits',len(slow),flush=True)


def reconcile_evidence():
    """Join immutable original questions and reconcile inclusive measurements."""
    baseline=RETAINED/'application-replay-pushed-complete';current=OUT/'final-round39-replay'
    transitions={};losses=[];new_discrepancies=[];bindings={}
    def lines(path):
        bindings[str(path)]=sha(path)
        with gzip.open(path,'rt')as stream:return [json.loads(line)for line in stream]
    for path in sorted(current.glob('*.inputs.json.gz')):
        name=path.name.removesuffix('.inputs.json.gz')
        before_inputs={r['id']:r for r in read(baseline/path.name)['witnesses']}
        before={r['id']:r for r in lines(baseline/(name+'.predictions.jsonl.gz'))}
        after={r['id']:r for r in lines(current/(name+'.predictions.jsonl.gz'))}
        tally=Counter()
        for row in read(path)['witnesses']:
            old=before_inputs[row['id']];a=before[row['id']];b=after[row['id']]
            if row['args']!=old['args']or row['expected']!=old['expected']or row['args']!=a['args']or row['args']!=b['args']:
                raise ValueError('changed original registered question')
            was=a['calc_func_prediction']==row['expected'];now=b['calc_func_prediction']==row['expected']
            tally[str((was,now))]+=1
            if was and not now:losses.append(dict(function=name,**row,before=a,after=b))
            if a['calc_func_prediction'].startswith('fault:')and not b['calc_func_prediction'].startswith('fault:')and not now:
                new_discrepancies.append(dict(function=name,**row,before=a['calc_func_prediction'],after=b['calc_func_prediction']))
        transitions[name]=tally
    before=read(OUT/'before-corpus/summary.json');after=read(OUT/'final-corpus/summary.json')
    old={w['name']:w for w in before['workbooks']};new={w['name']:w for w in after['workbooks']}
    if old.keys()!=new.keys()or any(old[k]['sha256']!=new[k]['sha256']for k in old):raise ValueError('changed corpus packages')
    changed=[dict(workbook=k,before=old[k]['counts'],after=new[k]['counts'])for k in old if old[k]['counts']!=new[k]['counts']]
    old_bank=read(OUT/'before-bank.json');new_bank=read(OUT/'final-bank.json')
    if old_bank!=new_bank:raise ValueError('bank observations changed; inspect differences before claiming equality')
    frozen('final-reconciliation',dict(registered_round39=dict(transitions=transitions,previous_matches_lost=losses,
           formerly_refused_discrepancies=new_discrepancies,inputs_and_traces=bindings),
           corpus=dict(before=before['totals'],after=after['totals'],changed_workbooks=changed,unchanged_packages=len(old),
                       before_run=before['run'],after_run=after['run']),
           bank=dict(all_fields_equal=True,totals=new_bank['totals'],functions=len(new_bank['functions'])),
           limits='Original IDs, operands and expected answers are joined before transition counts. Inclusive corpus counts and all bank fields reconcile; workbook totals alone do not prove per-cell transition counts.'))
    print('R39 gained',sum(v.get(str((False,True)),0)for v in transitions.values()),'lost',len(losses),
          'new numeric discrepancies',len(new_discrepancies),'corpus',after['totals']['counts']['pass'],flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('step', choices=['reanalysis', 'normal-map', 'normal-refinement','confidence-reanalysis','slope-map','slope-analysis','slope-impulses','impulse-analysis','slope-stores','short-proof','final-add-proof','final-product-proof','census','resume-census','search-structures','reset-boundaries','seasonal-initialization','seasonal-exact-proof','resume-seasonal-exact','partial-cycles','detection','detection-followup','calendar','calendar-followup','audit-evidence','reconcile-evidence'])
    args=parser.parse_args()
    {'reanalysis':reanalysis, 'normal-map':quantile_mapping, 'normal-refinement':quantile_refinement,
     'confidence-reanalysis':confidence_reanalysis,'slope-map':slope_mapping,'slope-analysis':slope_analysis,
     'slope-impulses':slope_impulses,'impulse-analysis':impulse_analysis,'slope-stores':slope_stores,'short-proof':short_proof,'final-add-proof':final_add_proof,'final-product-proof':final_product_proof,'census':parameter_census,'resume-census':lambda:parameter_census(True),'search-structures':search_structures,'reset-boundaries':reset_boundaries,'seasonal-initialization':seasonal_initialization,'seasonal-exact-proof':seasonal_exact_proof,'resume-seasonal-exact':lambda:seasonal_exact_proof(True),'partial-cycles':lambda:seasonal_initialization(True),'detection':detection_proof,'detection-followup':detection_followup,'calendar':calendar_proof,'calendar-followup':calendar_followup,'audit-evidence':audit_evidence,'reconcile-evidence':reconcile_evidence}[args.step]()
