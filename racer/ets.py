"""ETS identification from tier-1 Excel answers and IEEE754. Original MIT code.

STAT-supplied parameters prove recursion conditionally, never the optimiser.
Exact comparisons stay exact. No tier-3/4 implementation source is used.
"""
import argparse
from collections import Counter
import gzip
import itertools
import json
import math
from pathlib import Path
import random
from .fingerprint import CachedService, load_service, bits, value, sha, write

ROOT=Path(__file__).resolve().parents[1]
CHECKOUT=ROOT.parent/'DnaCalc'
OUT=Path('/home/ubuntu/codex-runs/round39')

def read(path):
    with (gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)) as stream:
        return json.load(stream)

def slope(y,kind='centered'):
    n=len(y); x=list(range(n)); mx=(n-1)/2; my=sum(y)/n
    if kind=='centered':return sum((i-mx)*(v-my) for i,v in enumerate(y))/sum((i-mx)**2 for i in x)
    if kind=='centered_y':return sum((i-mx)*v for i,v in enumerate(y))/sum((i-mx)**2 for i in x)
    if kind=='moments0':return (n*sum(i*v for i,v in enumerate(y))-sum(x)*sum(y))/(n*sum(i*i for i in x)-sum(x)**2)
    if kind=='moments1':return (n*sum((i+1)*v for i,v in enumerate(y))-sum(i+1 for i in x)*sum(y))/(n*sum((i+1)**2 for i in x)-sum(i+1 for i in x)**2)
    if kind=='mean0':return (sum(i*v for i,v in enumerate(y))-n*mx*my)/(sum(i*i for i in x)-n*mx*mx)
    if kind=='mean1':return (sum((i+1)*v for i,v in enumerate(y))-n*(mx+1)*my)/(sum((i+1)**2 for i in x)-n*(mx+1)**2)
    raise ValueError(kind)

SLOPES=('centered','centered_y','moments0','moments1','mean0','mean1')

def sequential_total(values):
    result=0.
    for v in values:result=result+v
    return result

def initial_trend(y):
    n=float(len(y));sx=sequential_total(float(i+1)for i in range(len(y)))
    sy=sequential_total(y);sxx=sequential_total(float(i+1)**2 for i in range(len(y)))
    sxy=sequential_total(float(i+1)*v for i,v in enumerate(y))
    return (n*sxy-sx*sy)/(n*sxx-sx*sx)

def numeric_fit(y,a,b,trend=None):
    """Fixed binary64 innovations tree, including two epsilon seasonal states."""
    level=y[0];trend=slope(y) if trend is None else trend
    seasons=[0.,0.];errors=[0.];predictions=[y[0]]
    for j,v in enumerate(y[1:],1):
        q=j%2;p=(level+trend)+seasons[q];e=v-p
        level=level+(trend+a*e)
        trend=trend+b*e
        seasons[q]=seasons[q]+2.**-52*e
        errors.append(e);predictions.append(p)
    center=(seasons[0]+seasons[1])/2
    level=level+center;seasons=[s-center for s in seasons]
    return level,trend,seasons,errors,predictions

def numeric_forecast(y,a,b,h,trend=None):
    l,t,s,e,p=numeric_fit(y,a,b,trend)
    return (l+h*t)+s[(len(y)-1+int(h))%2]

def numeric_statistics(y,a,b,trend=None):
    l,t,s,e,p=numeric_fit(y,a,b,trend)
    return statistics_from_errors(y,e,p)

def statistics_from_errors(y,e,p):
    n=len(y);tail=n if n<=10 else min(int(.4*n),10);start=n-tail
    mae=sequential_total(abs(v)for v in e[start:])/tail
    rmse=math.sqrt(sequential_total(v*v for v in e[start:])/tail)
    smape=sequential_total(abs(v)/((abs(y0)+abs(p0))/2)if y0 or p0 else 0. for v,y0,p0 in zip(e[start:],y[start:],p[start:]))/tail
    prefix=start or n
    scale=sequential_total(abs(v-u)for u,v in zip(y[:prefix],y[1:prefix]))/(prefix-1)
    return [mae/scale if scale else 0.,smape,mae,rmse]

def fit(y,a,b,sk='centered',lk='weighted',tk='innovations'):
    level=y[0]; trend=slope(y,sk); errors=[0.];predictions=[y[0]]
    for v in y[1:]:
        p=level+trend; e=v-p;old=level
        level=a*v+(1-a)*p if lk=='weighted' else p+a*e
        trend=trend+b*e if tk=='innovations' else b*(level-old)+(1-b)*trend
        errors.append(e);predictions.append(p)
    return level,trend,errors,predictions

def metrics(y,e,p,start=0):
    n=len(y)-start
    mae=sum(abs(v) for v in e[start:])/n
    rmse=math.sqrt(sum(v*v for v in e[start:])/n)
    smape=sum(abs(v)/((abs(y0)+abs(p0))/2) if y0 or p0 else 0. for v,y0,p0 in zip(e[start:],y[start:],p[start:]))/n
    scale=sum(abs(v-u) for u,v in zip(y,y[1:]))/(len(y)-1)
    return [mae/scale if scale else 0.,smape,mae,rmse]

def service(offline=False):
    c,_=load_service(CHECKOUT)
    return c,CachedService(None if offline else c.Client(timeout=300),OUT/'service-cache',offline=offline)

def ask_cases(c,s,cases,what,variants):
    rows=[];labels=[]
    for i,z in enumerate(cases):
        y=c.number_array([[v]for v in z['y']]);t=c.number_array([[float(j+1)]for j in range(len(z['y']))]);season=c.number(float(z.get('season',0)))
        for v in variants:
            if what=='FORECAST.ETS.STAT':args=[y,t,c.number(float(v)),season]
            elif what=='FORECAST.ETS':args=[c.number(len(z['y'])+v),y,t,season]
            elif what=='FORECAST.ETS.CONFINT':
                h,conf=v;args=[c.number(len(z['y'])+h),y,t,c.number(conf),season]
            else:args=[y,t]
            rows.append(args);labels.append([i,v])
    answers=s.ask_rows(what,rows,chunk_size=100)
    return dict(cases=cases,function=what,labels=labels,answers=answers,receipts=s.receipts,producer_sha256=sha(__file__))

def grouped(doc):
    result=[{}for _ in doc['cases']]
    for (i,v),a in zip(doc['labels'],doc['answers']):result[i][str(v)]=a
    return result

def reanalysis():
    retained=ROOT/'evidence/round38/retained'
    design=read(ROOT/'evidence/round37/retained/ets-design.json.gz')['cases']
    stats=read(retained/'ets/FORECAST.ETS.STAT.answers.json.gz')['records']
    conf=read(retained/'ets/FORECAST.ETS.CONFINT.answers.json.gz')['records']
    bystat={};byconf={}
    for r in stats:bystat.setdefault(r['case'],{})[r['variant']]=r['answer']
    for r in conf:byconf.setdefault(r['case'],{})[r['variant']]=r['answer']
    results=[]
    for i,z in enumerate(design):
        if z['seasonality']!=0 or z['data_completion']!=1 or z['aggregation']!=1 or len(z['y']) not in (8,12,24) or z['pattern'] not in ('noisy','seasonal'):continue
        y=z['y'];st=bystat[i];a=value(st[1]['bits']);b=value(st[2]['bits']);l,t,e,p=fit(y,a,b)
        target=[value(st[k]['bits'])for k in range(4,8)]
        windows=[(j,metrics(y,e,p,j))for j in range(len(y))]
        start,m=min(windows,key=lambda r:max(abs(u-v)for u,v in zip(r[1][1:],target[1:])))
        actual=value(byconf[i][.95]['bits']);pred=1.959963984540054*math.sqrt(sum(v*v for v in e)/len(y))
        results.append(dict(case=i,n=len(y),pattern=z['pattern'],alpha=a,beta=b,tail=len(y)-start,metrics_actual=target,metrics_predicted=m,confidence_actual=actual,confidence_predicted=pred,max_relative=max(abs(u-v)/max(1.,abs(u))for u,v in zip(target[1:]+[actual],m[1:]+[pred]))))
    if len(results)!=6 or any(r['max_relative']>1e-11 for r in results):raise ValueError(results)
    write(OUT/'step1-reanalysis.json',dict(results=results,producer_sha256=sha(__file__)))
    print(json.dumps(results,indent=2))

def arithmetic(offline=False):
    rng=random.Random(3902);cases=[]
    for i in range(64):
        n=rng.randrange(5,33)
        if i<20:
            base=rng.choice([1.,2.,math.pi,1e8]);y=[value(f'{int(bits(base),16)+rng.randrange(-10,11):016x}')for _ in range(n)]
        elif i<40:
            exp=rng.randrange(-400,400);y=[math.ldexp(rng.uniform(-8.,8.),exp) for _ in range(n)]
        else:y=[rng.uniform(-30,30)+j*rng.uniform(-2,2)for j in range(n)]
        cases.append(dict(y=y,season=0))
    c,s=service(offline);d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,6,7]);write(OUT/'step2-stat.json',d)
    f=ask_cases(c,s,cases,'FORECAST.ETS',[0.,1.,2.,5.]);write(OUT/'step2-forecasts.json',f)
    analyze_arithmetic(d,f)

def analyze_arithmetic(d,f):
    st=grouped(d);fo=grouped(f);results=[]
    for sk,lk,tk in itertools.product(SLOPES,('weighted','error'),('innovations','classical')):
        failures=[];exact=0
        for i,z in enumerate(d['cases']):
            if st[i]['1']['type']!='number' or st[i]['2']['type']!='number':continue
            l,t,e,p=fit(z['y'],value(st[i]['1']['bits']),value(st[i]['2']['bits']),sk,lk,tk)
            for h in (0.,1.,2.,5.):
                actual=fo[i][str(h)];pred=bits(l+h*t)
                if actual.get('bits')==pred:exact+=1
                else:failures.append(dict(case=i,h=h,predicted=pred,actual=actual))
        results.append(dict(slope=sk,level=lk,trend=tk,exact=exact,failures=failures))
    write(OUT/'step2-analysis.json',dict(results=results,producer_sha256=sha(__file__)))
    print([(r['slope'],r['level'],r['trend'],r['exact'],len(r['failures']))for r in results])

def parameter_candidate(y):
    """Frozen discovery candidate; the fresh receipt records every refutation."""
    n=len(y);start=0 if n<=10 else n-min(int(.4*n),10)
    scores=[]
    for a,b0 in itertools.product([.9,.1,.25,.5,.75,.002,.126],[.001,.1,.25,.5,.75,.9]):
        b=min(b0,a-.001);l=y[0];t=slope(y);errors=[0.]
        for j,v in enumerate(y[1:],1):
            p=l+t;e=v-p;errors.append(e)
            if not start or j<start:l=a*v+(1-a)*p;t+=b*e
            else:l=p
        scores.append((sum(v*v for v in errors),a,b))
    return min(scores,key=lambda r:r[0])

def fresh_parameters(offline=False):
    rng=random.Random(3945);cases=[];predictions=[]
    for i in range(100):
        n=rng.randrange(7,33);kind=i%5
        y=[rng.uniform(-1,1)+(j/5)**kind for j in range(n)] if kind<3 else [rng.uniform(-5,5)+5*math.sin(j/rng.choice([2.,3.]))for j in range(n)]if kind==3 else [rng.uniform(-5,5)for j in range(n)]
        cases.append(dict(y=y,season=0));best=parameter_candidate(y)
        predictions.append(dict(alpha=bits(best[1]),beta=bits(best[2]),objective=best[0]))
    path=OUT/'step45-frozen-100.json'
    frozen=dict(cases=cases,predictions=predictions,producer_sha256=sha(__file__),method='Hybrid SSE: update before tail, then project without updates; preselected coarse alpha and clamped beta lists. Predictions frozen before this request.')
    if path.exists():
        if read(path)['predictions']!=predictions:raise ValueError('frozen predictions changed')
    else:write(path,frozen)
    c,s=service(offline);d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3,6,7]);write(OUT/'step45-fresh-100.json',d)
    actual=grouped(d);failures=[]
    for i,p in enumerate(predictions):
        if actual[i]['1'].get('bits')!=p['alpha'] or actual[i]['2'].get('bits')!=p['beta']:
            failures.append(dict(case=i,prediction=p,actual=[actual[i]['1'],actual[i]['2']]))
    write(OUT/'step45-fresh-analysis.json',dict(cases=len(cases),failures=failures,producer_sha256=sha(__file__)))
    print('parameter failures',len(failures),failures[:8])

def parameter_mapping(offline=False):
    rng=random.Random(39551);cases=[]
    for n,a,b0 in itertools.product([8,12,24,32],[.05,.2,.5,.8],[.01,.05,.2,.5]):
        b=min(a,b0);l=10.;t=.25;y=[l]
        for j in range(1,n):
            e=rng.uniform(-1.,1.);p=l+t;y.append(p+e);l=p+a*e;t+=b*e
        cases.append(dict(y=y,season=0,design=[n,a,b]))
    for n,noise in itertools.product([8,12,24,32],[0.,.1,.25,.5,1.,2.,4.,8.]):
        y=[10.+.05*j*j+noise*rng.uniform(-1.,1.)for j in range(n)]
        cases.append(dict(y=y,season=0,design=['quadratic',n,noise]))
    write(OUT/'step5-mapping-design.json',dict(cases=cases,producer_sha256=sha(__file__)))
    c,s=service(offline);d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);write(OUT/'step5-mapping-stat.json',d)
    print(Counter((value(r['1']['bits']),value(r['2']['bits']))for r in grouped(d)))

def smooth_leaves(offline=False):
    rng=random.Random(39226);cases=[]
    for i in range(100):
        n=rng.choice([7,8,11,12,13,16,17,24,25,32]);trend=128.
        half=[rng.randrange(-4,5)for _ in range((n+1)//2)]
        noise=half+half[-2::-1]if n%2 else half+half[::-1];total=sum(noise)
        curvature=rng.choice([1/256,1/64,1/16]);noise_scale=rng.choice([0.,1/64,1/16,.25])
        y=[trend*j+curvature*(3*(2*j-n+1)**2-(n*n-1))+noise_scale*(n*noise[j]-total)for j in range(n)]
        if any(slope(y,k)!=trend for k in SLOPES):raise ValueError('inexact designed trend')
        cases.append(dict(y=y,season=0,trend=trend,design=[curvature,noise_scale]))
    write(OUT/'step2-smooth-design.json',dict(cases=cases,producer_sha256=sha(__file__)))
    c,s=service(offline);d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',range(1,9));write(OUT/'step2-smooth-stat.json',d)
    print('smooth parameter leaves',Counter((value(r['1']['bits']),value(r['2']['bits']))for r in grouped(d)))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['reanalysis','arithmetic','parameters','mapping','smooth']);parser.add_argument('--offline',action='store_true');args=parser.parse_args()
    OUT.mkdir(exist_ok=True,parents=True)
    if args.command=='reanalysis':reanalysis()
    elif args.command=='arithmetic':arithmetic(args.offline)
    elif args.command=='parameters':fresh_parameters(args.offline)
    elif args.command=='mapping':parameter_mapping(args.offline)
    else:smooth_leaves(args.offline)
