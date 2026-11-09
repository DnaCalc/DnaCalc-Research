"""R47 ETS shape research. Original MIT; tier 1 Excel answers and IEEE754.

Extends the existing audited service, ETS recurrence, plateau and replay tools.
No restricted implementation source. Predictions are written before questions.
"""
import argparse
from collections import Counter
import ctypes
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import time

import numpy as np
from . import ets, round41
from .fingerprint import bits, value
from .ets_optimizer import GRID
from .question_designer import vector_path

ROOT = Path(__file__).resolve().parents[1]
OUT = Path('/home/ubuntu/codex-runs/work47')
POOL = GRID + [(.833333, .0833333), (.998, .249)]
# The four strict precedence observations on R41 adjacent boundaries constrain
# only these ties. Preserve earlier relative order everywhere else.
RANKED_POOL = [p for p in POOL if p != (.126,.001)]
RANKED_POOL.insert(RANKED_POOL.index((.666667,.5)),(.126,.001))


def save(name, obj):
    with (OUT/(name+'.json')).open('x') as f:
        json.dump(obj, f, separators=(',', ':'), allow_nan=False)


def frozen(name, obj):
    save(name+'-frozen', dict(obj, source=Path(__file__).read_text(),
                             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()))


def service():
    round41.OUT = OUT
    return round41.service(priority='interactive')


def ask(c, s, cases, function, variants):
    rows=[]; labels=[]
    for i,z in enumerate(cases):
        y=c.number_array([[v] for v in z['y']])
        t=c.number_array([[float(j+1)] for j in range(len(z['y']))])
        season=c.number(float(z.get('season',0)))
        for k in variants:
            if function=='FORECAST.ETS.STAT': args=[y,t,c.number(float(k)),season]
            elif function=='FORECAST.ETS': args=[c.number(len(z['y'])+k),y,t,season]
            elif function=='FORECAST.ETS.CONFINT': args=[c.number(len(z['y'])+k[0]),y,t,c.number(k[1]),season]
            else: args=[y,t]
            rows.append(args); labels.append([i,k])
    answers=s.ask_rows(function,rows,chunk_size=96)
    return dict(cases=cases,function=function,labels=labels,answers=answers,receipts=s.receipts)


def retained(round_name, name):
    paths=list((ROOT/'evidence'/round_name/'retained').rglob(name+'.json.gz'))
    if len(paths)!=1: raise ValueError((name,paths))
    return ets.read(paths[0])


C_SOURCE=r'''
#include <math.h>
static double xa(double a,double b){volatile double z=(long double)a+(long double)b;return z;}
static double xs(double a,double b){volatile double z=(long double)a-(long double)b;return z;}
static double xm(double a,double b){volatile double z=(long double)a*(long double)b;return z;}
static double xd(double a,double b){volatile double z=(long double)a/(long double)b;return z;}
double reg(const double*y,int n){double sy=0,sxy=0,sx=0,sxx=0;for(int j=0;j<n;j++){double x=j+1;sy+=y[j];sxy+=x*y[j];sx+=x;sxx+=x*x;}return (n*sxy-sx*sy)/(n*sxx-sx*sx);}
static double slope(const double*y,int n){if(n>=7)return reg(y,n);double z[7],t=reg(y,n);for(int j=0;j<7;j++)z[j]=j<n?y[j]:y[j%2]+(j-j%2)*t;return reg(z,7);}
void run(const double*y,int n,double a,double b,int validation,double initial,double*out){
 double l=y[0],t=initial,se[2]={0,0},tr[128],err[128]={0},pred[128],sse=0;pred[0]=y[0];
 for(int j=1;j<n;j++){int q=j%2;double p=(l+t)+se[q],e=y[j]-p;err[j]=e;pred[j]=p;sse+=e*e;
 if(validation && j>=validation)l+=t;else{l+=t+a*e;t+=b*e;se[q]+=0x1p-52*e;}tr[j-1]=t;}
 int w=n<=10?n:(n+1)/2,dw=n<=10?n:n/2;double total=0,sq=0;
 for(int j=n-w;j<n-1;j++){total=xa(total,tr[j]);sq=xa(sq,xm(tr[j],tr[j]));}
 double m=xd(total,dw-1),var=xs(xd(sq,dw-1),xm(m,m));int reset=fabs(m)<xd(xm(sqrt(var),dw),dw-1);
 if(reset){l=y[n-1]-se[(n-1)%2];t=initial;b=.001;}
 double center=(se[0]+se[1])/2;l+=center;se[0]-=center;se[1]-=center;
 out[0]=sse;out[1]=a;out[2]=b;out[3]=l;out[4]=t;out[5]=se[0];out[6]=se[1];out[7]=reset;
 int count=n<=10?n:(int)(n*.4)<10?(int)(n*.4):10,start=n-count,pre=start?start:n;double mae=0,mse=0,smape=0,scale=0;
 for(int j=start;j<n;j++){mae+=fabs(err[j]);mse+=err[j]*err[j];double den=(fabs(y[j])+fabs(pred[j]))/2;smape+=den==0?0:fabs(err[j])/den;}
 for(int j=1;j<pre;j++)scale+=fabs(y[j]-y[j-1]);scale/=pre-1;
 out[8]=scale==0?0:mae/count/scale;out[9]=smape/count;out[10]=mae/count;out[11]=sqrt(mse/count);
}
void search(const double*y,int n,const double*pool,int np,double*out){
 double z[128];int m=n<7?7:n;double t0=reg(y,n);for(int j=0;j<m;j++)z[j]=j<n?y[j]:y[j%2]+(j-j%2)*t0;
 int tail=m<=10?m:(int)(m*.4)<10?(int)(m*.4):10,start=m-tail;double initial=slope(z,m),best=INFINITY,tmp[12];int chosen=0;
 for(int k=0;k<np;k++){run(z,m,pool[2*k],pool[2*k+1],start,initial,tmp);out[12+k]=tmp[0];if(tmp[0]<best){best=tmp[0];chosen=k;}}
 run(y,n,pool[2*chosen],pool[2*chosen+1],0,slope(y,n),out);out[12+np]=chosen;
}
double finish(double l,double t,double s,double h){return xa(l+h*t,s);}
'''


class Native:
    def __init__(self):
        p=OUT/'ets47.c'
        if not p.exists():p.write_text(C_SOURCE)
        elif p.read_text()!=C_SOURCE:raise ValueError('changed native source')
        so=OUT/'ets47.so'
        if not so.exists():subprocess.run(['cc','-O0','-fno-fast-math','-ffp-contract=off','-shared','-fPIC',str(p),'-o',str(so),'-lm'],check=True,timeout=120)
        self.lib=ctypes.CDLL(str(so));ptr=ctypes.POINTER(ctypes.c_double)
        self.lib.search.argtypes=[ptr,ctypes.c_int,ptr,ctypes.c_int,ptr]
        self.lib.run.argtypes=[ptr,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_double,ptr]
        self.lib.search.restype=None;self.lib.run.restype=None
        self.lib.finish.argtypes=[ctypes.c_double]*4;self.lib.finish.restype=ctypes.c_double

    def fit(self,y,pool=POOL):
        if not 2<=len(y)<=128 or not pool or any(len(p)!=2 for p in pool):raise ValueError('bounded history and parameter pairs required')
        ys=np.array(y,dtype='float64');ps=np.array(pool,dtype='float64');out=np.zeros(13+len(pool))
        ptr=ctypes.POINTER(ctypes.c_double)
        self.lib.search(ys.ctypes.data_as(ptr),len(y),ps.ctypes.data_as(ptr),len(pool),out.ctypes.data_as(ptr))
        return out

    def forecast(self,y,fit,h):
        return self.lib.finish(fit[3],fit[4],fit[5+(len(y)-1+int(h))%2],h)

    def conditional(self,y,a,b,t):
        if not 2<=len(y)<=128:raise ValueError('bounded history required')
        ys=np.array(y,dtype='float64');out=np.zeros(12);ptr=ctypes.POINTER(ctypes.c_double)
        self.lib.run(ys.ctypes.data_as(ptr),len(y),a,b,0,t,out.ctypes.data_as(ptr))
        return out


def search_boundaries():
    native=Native();d=retained('round41','ets-prior-boundaries-frozen');r=retained('round41','ets-prior-boundaries-result');rows=[]
    for path in r['paths']:
        if path['focused_edge']is None:continue
        z=d['paths'][path['path']];right=z['right']if not path.get('anchor')else z['alternatives'][path['anchor']-1]
        obs=dict(path['observations'])
        for index in path['focused_edge']:
            y=vector_path(z['left'],right,index);f=native.fit(y);scores=f[12:12+len(POOL)];minimum=np.flatnonzero(scores==min(scores)).tolist()
            compatible=[]
            for k in minimum:
                ff=native.fit(y,[POOL[k]])
                if [bits(ff[1]),bits(ff[2]),bits(2.**-52)]==obs[index]:compatible.append(k)
            rows.append(dict(path=path['path'],index=index,y=y,answer=obs[index],minimum=minimum,compatible=compatible,prediction=[bits(f[1]),bits(f[2]),bits(2.**-52)]))
    save('search-boundaries',dict(pool=POOL,rows=rows,contradictions=sum(not r['compatible']for r in rows)))
    print('retained boundaries',len(rows),'incompatible',sum(not r['compatible']for r in rows),flush=True)


def search_fresh():
    native=Native();cases=list(round41.histories(512,470101));pred=[]
    for z in cases:
        f=native.fit(z['y']);pred.extend([bits(f[1]),bits(f[2]),bits(2.**-52)])
    frozen('search-fresh',dict(cases=cases,pool=POOL,predictions=pred,method='Fixed25 table, first strict minimum, R46 raw moment reset; eight fixed history families including exponent and one-bit controls.'))
    c,s=service();d=ask(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);save('search-fresh-answers',d)
    failures=[dict(index=i,case=i//3,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
    save('search-fresh-result',dict(selections=len(pred),failures=failures))
    print('search fresh',len(pred),'failures',len(failures),failures[:4],flush=True)


def short_map():
    native=Native();rng=random.Random(470201);cases=[];forecasts=[];stats=[]
    for n in range(3,8):
        for i in range(64):
            ex=rng.randrange(-350,351)
            y=[math.ldexp(rng.randrange(-512,513)/16,ex)for _ in range(n)]
            cases.append(dict(y=y,season=0,kind='coarse_short',n=n))
            f=native.fit(y);stats.extend([bits(f[1]),bits(f[2]),bits(2.**-52)])
            forecasts.extend(bits(native.forecast(y,f,h))for h in [0.,1.,2.,5.])
    frozen('short-map',dict(cases=cases,stat_predictions=stats,forecast_predictions=forecasts,pool=POOL,
         method='Raw moment OLS, first-two-phase extension, second raw OLS. Coarse dyadic inputs keep moment sums exact; quotient and extension may round. Strict binary64 fixed candidate.'))
    c,s=service()
    for name,fun,variants,pred in [('stat','FORECAST.ETS.STAT',[1,2,3],stats),('forecast','FORECAST.ETS',[0.,1.,2.,5.],forecasts)]:
        d=ask(c,s,cases,fun,variants);save('short-map-'+name,d)
        fail=[dict(index=i,case=i//len(variants),prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
        save('short-map-'+name+'-result',dict(selections=len(pred),failures=fail));print('short',name,len(pred),'failed',len(fail),fail[:3],flush=True)


def detector_map():
    from .ets_round40 import padded_history
    from .ets_seasonal import detection_candidate
    rng=random.Random(470102);cases=[]
    for n in [5,6,7,8,10,11,12,15,16,20,24,31,32]:
        for k in range(24):
            ex=[-350,-100,0,100,350][k%5];p=2+k%7;phase=k%p;trend=[0.,.125,16.][k%3]
            pattern=[float(p-1 if j==0 else -1)for j in range(p)]
            y=[pattern[(j+phase)%p]+trend*j for j in range(n)]if k<12 else[rng.randrange(-32,33)/8 for _ in range(n)]
            cases.append(dict(y=[math.ldexp(v,ex)for v in y],n=n,p=p,kind='periodic'if k<12 else'random'))
    programs={}
    # Distinct fixed programs extend R40 with the measured short extension,
    # first/last/global ACF peaks, and a full-lag scan before the cycle limit.
    for pre,norm,scope,order in itertools.product(['mean','linear','difference'],['biased','unbiased','pearson'],['all','padded','prefix60'],['first','last','global']):
        key=str((pre,norm,scope,order));pred=[]
        for z in cases:
            y=z['y'];y=padded_history(y)if len(y)<7 else y
            if scope=='prefix60':y=y[:max(4,len(y)-(len(y)if len(y)<=10 else min(int(len(y)*.4),10)))]
            try:
                _,corr=detection_candidate(y,pre,norm,'first',0,0.)
                peaks=[j for j in range(2,len(corr)-1)if corr[j]>0 and corr[j]>corr[j-1]and corr[j]>=corr[j+1]]
                if scope=='padded':peaks=[j for j in peaks if j<=2*len(y)//5]
                chosen=(peaks[0]if order=='first'else peaks[-1]if order=='last'else max(peaks,key=lambda j:corr[j]))if peaks else 0
                if chosen>2*len(y)//5:chosen=0
            except (ValueError,OverflowError,ZeroDivisionError):chosen=None
            pred.append(chosen)
        programs[key]=pred
    frozen('detector-map',dict(cases=cases,predictions=programs,method='Measured short extension plus fixed ACF programs; no fitted cutoff or polynomial. Exponents and trends expose preprocessing; full-lag selection separates search from cycle-limit rejection.'))
    c,s=service();d=ask(c,s,cases,'FORECAST.ETS.SEASONALITY',[0]);save('detector-map-answers',d)
    actual=[value(a['bits'])if 'bits'in a else a for a in d['answers']]
    result={k:dict(contradictions=sum(p!=a for p,a in zip(ps,actual)),first=next((dict(index=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(ps,actual))if p!=a),None))for k,ps in programs.items()}
    save('detector-map-result',dict(programs=result,survivors=[k for k,v in result.items()if not v['contradictions']]))
    print('detector programs',len(result),'survivors',sum(not v['contradictions']for v in result.values()),flush=True)


def short_stores():
    from .ets_loop import slope_graph,SLOPE_TAGS
    from .x87_programs import Native as GraphNative,build_library
    native=Native();library=OUT/'x87-ets47.so';build_library(library);graph_native=GraphNative(library)
    d=ets.read(OUT/'short-map-frozen.json');cases=d['cases'];actual=ets.read(OUT/'short-map-forecast.json')['answers']
    outputs=[a.get('bits')for a in actual];fits=[native.fit(z['y'])for z in cases];survivors=[];programs=[]
    for form in range(3):
        graphs={n:slope_graph(n,form,extend=True)for n in range(3,8)}
        arrays={n:np.array([[int(bits(v),16)for v in z['y']]for z in cases if len(z['y'])==n],dtype='uint64')for n in range(3,8)}
        for mask in range(256):
            ts=[]
            for n in range(3,8):
                values=graph_native.grid(graphs[n],arrays[n],0x133f,mask)
                ts.extend(value(f'{int(r[0]):016x}')for r in values)
            failures=[]
            for i,(z,f,t)in enumerate(zip(cases,fits,ts)):
                y=z['y'];raw=int(f[12+len(POOL)]);a,b=POOL[raw];ff=native.conditional(y,a,b,t)
                predicted=[bits(native.forecast(y,ff,h))for h in [0.,1.,2.,5.]]
                if predicted!=outputs[i*4:i*4+4]:failures.append(dict(case=i,prediction=predicted,answer=outputs[i*4:i*4+4],initial=bits(t)))
            programs.append(dict(form=form,mask=mask,contradictions=len(failures),first=failures[:1]))
            if not failures:survivors.append(dict(form=form,mask=mask,stores=[t for i,t in enumerate(SLOPE_TAGS)if mask&(1<<i)]))
    save('short-stores',dict(programs=programs,survivors=survivors,limits='Conditional on the previously matching raw search pair; every initial slope program is replayed with the shared recurrence and R46 reset. No rounded-constant extension is inferred.'))
    print('short slope complete survivors',len(survivors),survivors[:8],flush=True)


def seasonal_map():
    from .ets_seasonal import moving_initial,fit
    rng=random.Random(470301);cases=[]
    grid=[[a,b,.001]for a in [.751,.501,.251,.126]for b in [.001,a-.001]]
    grid += [[.002,.001,g]for g in [.125,.25,.5,.75]]
    grid += [[a,b,.999-a]for a in [.1,.5,.75,.9]for b in [.001,a-.001]]
    grid += [[.998,b,.001]for b in [.099,.249,.499,.749,.899]]
    grid=[[float(format(v,'.6g'))for v in row]for row in grid]
    params=[];inits=[]
    for i in range(128):
        p=2 if i<96 else 4;n=rng.choice([8,12,24])if p==2 else rng.choice([12,16,24]);ex=rng.randrange(-250,251)
        pattern=[4*rng.randrange(-32,33)for _ in range(p)]
        y=[128.*j+pattern[j%p]for j in range(n)]
        if i%2:y=[v+rng.randrange(-16,17)/16 for v in y]
        else:y[rng.randrange(n)]+=1/16
        y=[math.ldexp(v,ex)for v in y];init=moving_initial(y,p);inits.append(init)
        start=n-(n if n<=10 else min(int(.4*n),10))
        def score(par):
            a,b,g=par;l,t,se=init[0],init[1],init[2].copy();sse=0.
            for j,v in enumerate(y[1:],1):
                q=j%p;e=v-((l+t)+se[q]);sse+=e*e
                if start and j>=start:l+=t
                else:l+=t+a*e;t+=b*e;se[q]+=g*e
            return sse
        chosen=min(grid,key=score);params.extend(bits(v)for v in chosen)
        cases.append(dict(y=y,season=p,kind='seasonal_innovations'))
    frozen('seasonal-map',dict(cases=cases,grid=grid,initials=inits,parameter_predictions=params,
         method='R39 forward2p moving averages and per-phase counts, rounded fixed25 seasonal points and hybrid SSE. Positive trends separate reset. All programs fixed before answers.'))
    c,s=service();d=ask(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);save('seasonal-map-stat',d)
    fail=[dict(index=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(params,d['answers']))if p!=a.get('bits')]
    forecasts=[]
    for i,z in enumerate(cases):
        a,b,g=[value(d['answers'][3*i+j]['bits'])for j in range(3)]
        l,t,se,e,pr=fit(z['y'],z['season'],a,b,g,inits[i])
        forecasts.extend(bits(float(np.longdouble(l+h*t)+np.longdouble(se[(len(z['y'])-1+int(h))%z['season']])))for h in [0.,1.,2.,5.])
    frozen('seasonal-conditional',dict(cases=cases,forecast_predictions=forecasts,parameters=d,
          method='Reported parameters condition the same moving initializer and innovations recurrence. Forecast arithmetic uses the already measured final PC64 store.'))
    fo=ask(c,s,cases,'FORECAST.ETS',[0.,1.,2.,5.]);save('seasonal-map-forecast',fo)
    ff=[dict(index=i,case=i//4,prediction=p,answer=a)for i,(p,a)in enumerate(zip(forecasts,fo['answers']))if p!=a.get('bits')]
    save('seasonal-map-result',dict(parameter_failures=fail,forecast_failures=ff,selections=len(params)+len(forecasts),
         stop='Fixed search and initializer programs have independent counterexamples. No extra point or arithmetic piece follows from their match count.'))
    print('seasonal parameter failures',len(fail),'conditional forecast failures',len(ff),flush=True)


def length_map():
    rng=random.Random(470401);cases=[];forecast=[];stats=[];confidence=[]
    # This fixed affine subset is exactly the first112-cell cycle-limit block.
    # Short complete cycles still enter seasonal initialization in Excel.
    for i in range(192):
        ex=rng.randrange(-350,351);n=4;p=2;offset=rng.randrange(-1024,1025)/8;trend=rng.randrange(-64,65)/16
        if i<32:trend=0.
        y=[math.ldexp(offset+trend*j,ex)for j in range(n)]
        cases.append(dict(y=y,season=p,trend=math.ldexp(trend,ex),kind='short_affine'))
        forecast.extend(bits(math.ldexp(offset+trend*(n-1+h),ex))for h in [0.,1.,2.,5.])
        stats.extend(bits(v)for v in [.751,.001,.001,0.,0.,0.,0.,1.]);confidence.extend([bits(0.)]*3)
    frozen('length-map',dict(cases=cases,forecasts=forecast,statistics=stats,confidence=confidence,
         method='Exact four-observation affine subset with explicit period2; no seasonal component. Fixed zero-innovation candidate, exponents, signs and zero trends.'))
    c,s=service();result={}
    for name,fun,variants,pred in [('forecast','FORECAST.ETS',[0.,1.,2.,5.],forecast),('stat','FORECAST.ETS.STAT',list(range(1,9)),stats),('confidence','FORECAST.ETS.CONFINT',[(1.,.95),(5.,.95),(10.,.5)],confidence)]:
        d=ask(c,s,cases,fun,variants);save('length-map-'+name,d)
        failures=[dict(index=i,case=i//len(variants),prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
        result[name]=dict(selections=len(pred),failures=failures)
        print('length4 affine',name,'failed',len(failures),'of',len(pred),flush=True)
    save('length-map-result',result)


def short_period_boundary():
    rng=random.Random(470402);cases=[];pairs=[]
    for n in range(3,13):
        for p in [2,3,4]:
            for kind in ['affine','seasonal','noise']:
                pattern=[rng.randrange(-16,17)/4 for _ in range(p)]
                ex=rng.choice([-100,0,100]);y=[math.ldexp(32*j+(0 if kind=='affine'else pattern[j%p])+(rng.randrange(-16,17)/16 if kind=='noise'else 0),ex)for j in range(n)]
                pairs.append(dict(n=n,period=p,kind=kind))
                cases.extend([dict(y=y,season=0),dict(y=y,season=p)])
    frozen('short-period-boundary',dict(cases=cases,pairs=pairs,hypothesis='Periods at length4 use the nonseasonal fit; extend mapping to identify the length boundary. No outputs are used as predictions for a numerical fit.'))
    c,s=service();results={}
    for name,fun,variants in [('stat','FORECAST.ETS.STAT',[1,2,3]),('forecast','FORECAST.ETS',[0.,1.,2.,5.])]:
        d=ask(c,s,cases,fun,variants);save('short-period-boundary-'+name,d);width=len(variants)
        results[name]=[dict(pair=z,equal=d['answers'][2*i*width:(2*i+1)*width]==d['answers'][(2*i+1)*width:(2*i+2)*width],left=d['answers'][2*i*width:(2*i+1)*width],right=d['answers'][(2*i+1)*width:(2*i+2)*width])for i,z in enumerate(pairs)]
    save('short-period-boundary-result',results)
    for n in range(3,13):print('period fit equality n',n,{k:sum(v['equal']for v in rows if v['pair']['n']==n)for k,rows in results.items()},flush=True)


def length_fresh():
    # Independent answers after the parameter mapping. This exact affine region
    # is structural, not an answer-selected list of histories.
    rng=random.Random(470403);cases=[];forecasts=[];stats=[];confidence=[]
    for i in range(256):
        ex=[-400,400][i%2]if i<16 else rng.randrange(-390,391)
        offset=rng.randrange(-512,513)/16;trend=rng.randrange(-64,65)/32
        y=[math.ldexp(offset+trend*j,ex)for j in range(4)]
        cases.append(dict(y=y,season=2,kind='four_point_affine_holdout',trend=math.ldexp(trend,ex)))
        forecasts.extend(bits(math.ldexp(offset+trend*(3+h),ex))for h in [0.,1.,2.,5.,12.])
        stats.extend(bits(v)for v in [.9,.001,2.**-52,0.,0.,0.,0.,1.]);confidence.extend([bits(0.)]*3)
    frozen('length-fresh',dict(cases=cases,forecasts=forecasts,statistics=stats,confidence=confidence,method='Independent exact four-point affine histories, explicit period2, nonseasonal reported parameters. Exponent endpoints plus fresh dyadic slopes and offsets.'))
    c,s=service();results={}
    for name,fun,variants,pred in [('forecast','FORECAST.ETS',[0.,1.,2.,5.,12.],forecasts),('stat','FORECAST.ETS.STAT',list(range(1,9)),stats),('confidence','FORECAST.ETS.CONFINT',[(1.,.95),(5.,.95),(10.,.5)],confidence)]:
        d=ask(c,s,cases,fun,variants);save('length-fresh-'+name,d)
        fail=[dict(index=i,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
        results[name]=dict(selections=len(pred),failures=fail);print('length fresh',name,'fail',len(fail),'/',len(pred),flush=True)
    save('length-fresh-result',results)


def search_order_proof():
    native=Native();d=retained('round41','ets-prior-boundaries-frozen');r=retained('round41','ets-prior-boundaries-result');cases=[];seen=set()
    for path in r['paths']:
        if path['focused_edge']is None:continue
        z=d['paths'][path['path']];right=z['right']if not path.get('anchor')else z['alternatives'][path['anchor']-1]
        for center in path['focused_edge']:
            for offset in [-8,-4,-2,-1,0,1,2,4,8]:
                index=max(0,min(1<<52,center+offset));y=vector_path(z['left'],right,index);key=tuple(bits(v)for v in y)
                if key not in seen:cases.append(dict(y=y,season=0,path=path['path'],index=index));seen.add(key)
    cases.extend(round41.histories(256,470104));pred=[];detail=[]
    for z in cases:
        f=native.fit(z['y'],RANKED_POOL);pred.extend([bits(f[1]),bits(f[2]),bits(2.**-52)])
        scores=f[12:12+len(RANKED_POOL)];detail.append(dict(selected=int(f[-1]),ties=np.flatnonzero(scores==min(scores)).tolist()))
    frozen('search-order-proof',dict(cases=cases,pool=RANKED_POOL,predictions=pred,scores=detail,
         method='Fixed25 first strict minimum; move126/.001 before666667/.5 per retained exact-score tie. No other order inferred. Adjacent input mapping plus independent eight-family holdout.'))
    c,s=service();d=ask(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);save('search-order-proof-answers',d)
    fail=[dict(index=i,case=i//3,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
    save('search-order-proof-result',dict(selections=len(pred),failures=fail,unresolved='Unobserved table precedence and score-load neighbors remain undetermined; no continuous optimizer family is claimed excluded.'))
    print('search order',len(pred),'failed',len(fail),fail[:3],flush=True)


def short_period_fresh():
    rng=random.Random(470404);cases=[]
    for i in range(256):
        ex=rng.randrange(-350,351)
        if i<192:y=[math.ldexp(rng.randrange(-1024,1025)/32,ex)for _ in range(4)]
        else:
            raw=int(bits(math.ldexp(1.5,ex)),16);y=[value(f'{raw+rng.randrange(-16,17):016x}')for _ in range(4)]
        cases.extend([dict(y=y,season=0),dict(y=y,season=2)])
    frozen('short-period-fresh',dict(cases=cases,prediction='Every consecutive pair has identical answers. Four prepared observations and explicit period2 choose the nonseasonal fit even when the values are not affine.',method='Independent coarse dyadic and single-bit histories with exponent variation. Equality is a routing observation; existing nonseasonal numerical refusals are unchanged.'))
    c,s=service();result={}
    for name,fun,variants in [('stat','FORECAST.ETS.STAT',list(range(1,9))),('forecast','FORECAST.ETS',[0.,1.,2.,5.]),('confidence','FORECAST.ETS.CONFINT',[(1.,.5),(5.,.95)])]:
        d=ask(c,s,cases,fun,variants);save('short-period-fresh-'+name,d);w=len(variants)
        fail=[dict(pair=i,left=d['answers'][2*i*w:(2*i+1)*w],right=d['answers'][(2*i+1)*w:(2*i+2)*w])for i in range(256)if d['answers'][2*i*w:(2*i+1)*w]!=d['answers'][(2*i+1)*w:(2*i+2)*w]]
        result[name]=dict(pairs=256,selections=len(d['answers']),failures=fail);print('short period paired',name,'failed',len(fail),flush=True)
    save('short-period-fresh-result',result)


def short_automatic():
    rng=random.Random(470405);cases=[]
    for n in [2,3,4]:
        for i in range(64):
            ex=rng.randrange(-350,351);y=[math.ldexp(rng.randrange(-512,513)/16,ex)for _ in range(n)]
            cases.extend([dict(y=y,season=0),dict(y=y,season=1)])
    frozen('short-automatic',dict(cases=cases,prediction='Each automatic/nonseasonal pair agrees for2..4 prepared points, consistent with the already measured zero detected seasonality. This is routing evidence, independent of numerical initialization.'))
    c,s=service();result={}
    for name,fun,variants in [('stat','FORECAST.ETS.STAT',list(range(1,9))),('forecast','FORECAST.ETS',[0.,1.,2.,5.]),('confidence','FORECAST.ETS.CONFINT',[(1.,.5),(5.,.95)])]:
        d=ask(c,s,cases,fun,variants);save('short-automatic-'+name,d);w=len(variants)
        fail=[dict(pair=i,left=d['answers'][2*i*w:(2*i+1)*w],right=d['answers'][(2*i+1)*w:(2*i+2)*w])for i in range(len(cases)//2)if d['answers'][2*i*w:(2*i+1)*w]!=d['answers'][(2*i+1)*w:(2*i+2)*w]]
        result[name]=dict(pairs=len(cases)//2,selections=len(d['answers']),failures=fail);print('short automatic',name,'failed',len(fail),flush=True)
    save('short-automatic-result',result)


def fit_load_plateau():
    from .ets_loop import graph
    from .plateau_probe import probe
    from .x87_programs import Native as GraphNative,build_library
    native=Native();path=OUT/'x87-ets47.so'
    if not path.exists():build_library(path)
    gn=GraphNative(path)
    class PhaseNative:
        def grid(self,g,inputs,cw,mask):return gn.grid(g,inputs,cw,mask,inner_cw=0x133f)
    rng=random.Random(470105);groups={};attempts=12000
    for i in range(attempts):
        n=rng.choice([8,9,10,11,12,13,16,24,32]);ex=rng.randrange(-300,301)
        half=[rng.randrange(-256,257)/16 for _ in range((n+1)//2)]
        residual=half+half[-2::-1]if n%2 else half+half[::-1]
        y=[math.ldexp(128*j+v,ex)for j,v in enumerate(residual)];t=math.ldexp(128.,ex)
        if bits(ets.initial_trend(y))!=bits(t):raise ValueError('constructed OLS not exact')
        f=native.fit(y,RANKED_POOL)
        if f[7]:continue
        pair=int(f[-1]);key=(n,pair);g=groups.setdefault(key,[])
        if len(g)<64:g.append((y,t,f))
    found={};selected={};details=[]
    for (n,pair),rows in sorted(groups.items(),key=lambda kv:-len(kv[1])):
        if all((pair,k,d)in found for k in ['alpha','beta']for d in [-1,1]):continue
        a,b=RANKED_POOL[pair];gr=graph(n,form=3,phases=('publication',))
        for node,v in [(n+2,a),(n+3,b)]:gr['nodes'][node]={'op':'const','bits':bits(v)}
        # Use a one-step forecast; input layout remains y,h,initial trend.
        output=gr['outputs'][n%2];gr['outputs']=[output,output]
        ins=np.array([[int(bits(v),16)for v in y+[1.,t]]for y,t,_ in rows],dtype='uint64')
        shape=dict(graph=gr,constants={'alpha':n+2,'beta':n+3},mask=(1<<gr['sites'])-1,cw=0x123f)
        baseline=PhaseNative().grid(gr,ins,shape['cw'],shape['mask'])
        for index,(y,t,f)in enumerate(rows):
            if int(baseline[index,0])!=int(bits(native.forecast(y,f,1.)),16):raise ValueError('plateau graph disagrees with fixed fit')
        report=probe(shape,ins,PhaseNative(),keep=1);details.append(dict(n=n,pair=pair,search_rows=len(rows),constants=report))
        for constant in report:
            for neighbor in constant['neighbors']:
                key=(pair,constant['name'],neighbor['direction'])
                if key in found or not neighbor['examples']:continue
                example=neighbor['examples'][0];idx=next(i for i,row in enumerate(ins)if [f'{int(v):016x}'for v in row]==example['inputs'])
                y,t,f=rows[idx];identity=tuple(bits(v)for v in y);selected.setdefault(identity,dict(y=y,season=0,pair=pair,prediction=example['baseline']))
                found[key]=dict(pair=pair,parameter=constant['name'],direction=neighbor['direction'],example=example)
    chosen=list(selected.values());status=[dict(pair=i,parameters=RANKED_POOL[i],neighbors=[found.get((i,k,d),dict(parameter=k,direction=d,status='undetermined_within_search'))for k in ['alpha','beta']for d in [-1,1]])for i in range(len(RANKED_POOL))]
    frozen('fit-load-plateau',dict(cases=chosen,neighbors=status,groups=details,attempts=attempts,
         method='Existing plateau tool separates immediate binary64 neighbors in the fixed recurrence. Constructed exact OLS, selected raw table pair and no reset. This tests fitted-state loads; validation-score-only loads remain distinct and undetermined.'))
    if chosen:
        c,s=service();d=ask(c,s,chosen,'FORECAST.ETS',[1.]);save('fit-load-plateau-answers',d)
        st=ask(c,s,chosen,'FORECAST.ETS.STAT',[1,2]);save('fit-load-plateau-stat',st)
        fail=[dict(case=i,prediction=z['prediction'],answer=a)for i,(z,a)in enumerate(zip(chosen,d['answers']))if z['prediction']!=a.get('bits')]
        param_fail=[dict(case=i)for i,z in enumerate(chosen)if [a.get('bits')for a in st['answers'][2*i:2*i+2]]!=[bits(v)for v in RANKED_POOL[z['pair']]]]
        save('fit-load-plateau-result',dict(selections=len(chosen),forecast_failures=fail,parameter_failures=param_fail,separated_neighbors=len(found),total_neighbors=4*len(RANKED_POOL)))
        print('plateau selected',len(chosen),'separated',len(found),'forecast failures',len(fail),'parameter failures',len(param_fail),flush=True)


def search_finish():
    native=Native();d=ets.read(OUT/'search-order-proof-frozen.json');cases=[];pred=[];stats=[]
    for z in d['cases']:
        old=native.fit(z['y'],GRID);new=native.fit(z['y'],RANKED_POOL)
        if [bits(v)for v in old[1:3]]==[bits(v)for v in new[1:3]]:continue
        cases.append(z);pred.extend(bits(native.forecast(z['y'],new,h))for h in [0.,1.,2.,5.]);stats.extend(bits(v)for v in new[8:12])
    frozen('search-finish',dict(cases=cases,forecasts=pred,statistics=stats,pool=RANKED_POOL,
         method='All changed reported parameter pairs in the frozen order proof, final-state forward predictions before asking. Matching parameters alone never certifies their initializer or consumers.'))
    c,s=service();result={}
    for name,fun,variants,predictions in [('forecast','FORECAST.ETS',[0.,1.,2.,5.],pred),('stat','FORECAST.ETS.STAT',[4,5,6,7],stats)]:
        answers=ask(c,s,cases,fun,variants);save('search-finish-'+name,answers)
        fail=[dict(index=i,case=i//4,prediction=p,answer=a)for i,(p,a)in enumerate(zip(predictions,answers['answers']))if p!=a.get('bits')]
        result[name]=dict(selections=len(predictions),failures=fail);print('search finish',name,'fail',len(fail),'/',len(predictions),flush=True)
    save('search-finish-result',result)


def search_ci():
    cases=ets.read(OUT/'search-finish-frozen.json')['cases'];c,s=service()
    variants=[(h,conf)for h in [1.,2.,5.,10.]for conf in [.5,.95,math.nextafter(1.,0.)]]
    rows=[]
    for i,z in enumerate(cases):
        for h,confidence in variants:
            args=[len(z['y'])+h,[[v]for v in z['y']],[[float(j+1)]for j in range(len(z['y']))],confidence,0.]
            def encode(v):return [encode(x)for x in v]if isinstance(v,list)else'0x'+bits(v)
            rows.append(dict(id=f'r47-ci-{i}-{h}-{bits(confidence)}',args=encode(args)))
    inp=OUT/'search-ci.inputs.json';out=OUT/'search-ci.predictions.jsonl'
    with inp.open('x')as f:json.dump(dict(function='FORECAST.ETS.CONFINT',witnesses=rows),f,separators=(',',':'))
    exe=ROOT/'target/debug/symbolic-racer'
    provenance=json.loads(subprocess.check_output([str(exe),'provenance'],text=True,timeout=60))
    subprocess.run([str(exe),'invoke-inputs',str(inp),str(out)],check=True,timeout=300)
    with out.open()as stream:predictions=[json.loads(line)['calc_func_prediction']for line in stream]
    frozen('search-ci',dict(cases=cases,variants=variants,predictions=predictions,compiled_provenance=provenance,
          method='Registered application predictions frozen before fresh confidence answers. All changed reported-pair cases from the order proof.'))
    d=ask(c,s,cases,'FORECAST.ETS.CONFINT',variants);save('search-ci-answers',d)
    fail=[dict(index=i,case=i//len(variants),prediction=p,answer=a)for i,(p,a)in enumerate(zip(predictions,d['answers']))if p!='0x'+a.get('bits','')]
    save('search-ci-result',dict(selections=len(predictions),failures=fail));print('search CI',len(predictions),'failed',len(fail),fail[:3],flush=True)


def tie_constraints():
    native=Native();p=OUT/'route-before/FORECAST.ETS.STAT.inputs.json.gz';rows=ets.read(p)['witnesses'];observations={}
    for r in rows:
        a=r['args']
        if len(a)<4 or a[3]!='0x'+bits(0.)or len(a[0])<8 or len(a[0])>32:continue
        k=value(a[2][2:])
        if k not in [1.,2.]or not r['expected'].startswith('0x'):continue
        key=json.dumps(a[:2],separators=(',',':'));entry=observations.setdefault(key,dict(y=[value(v[0][2:])for v in a[0]],observed={}))
        entry['observed'][int(k)]=value(r['expected'][2:])
    constraints=set();ties=[];incompatible=[]
    for i,z in enumerate(observations.values()):
        if 1 not in z['observed']:continue
        f=native.fit(z['y'],POOL);scores=f[12:12+len(POOL)];ms=np.flatnonzero(scores==min(scores)).tolist()
        if len(ms)<2:continue
        a=z['observed'][1];b=z['observed'].get(2);cs=[j for j in ms if POOL[j][0]==a and (b is None or b==.001 or b==POOL[j][1])]
        row=dict(y=z['y'],observed=z['observed'],minimum=ms,compatible=cs)
        if not cs:incompatible.append(row)
        elif len(cs)==1:
            constraints.update((cs[0],j)for j in ms if j!=cs[0]);ties.append(row)
    remaining=set(range(len(POOL)));order=[]
    while remaining:
        ready=next((i for i in sorted(remaining)if not any(b==i and a in remaining for a,b in constraints)),None)
        if ready is None:break
        remaining.remove(ready);order.append(ready)
    save('tie-constraints',dict(histories=len(observations),pool=POOL,constraints=sorted(constraints),order=order,
         cycle=bool(remaining),ties=ties,incompatible=incompatible,
         limits='Necessary stable-first ordering constraints from exact-score ties only. Missing compatible minima keep the complete search/initialization program refuted; no numeric coefficient is chosen by counts.'))
    print('tie constraints',len(constraints),'cycle',bool(remaining),'incompatible ties',len(incompatible),'order',[POOL[i]for i in order],flush=True)


def tie_fresh():
    native=Native();rank=ets.read(OUT/'tie-constraints.json')
    if rank['cycle']:raise ValueError('no stable-first order')
    pool=[POOL[i]for i in rank['order']];cases=list(round41.histories(512,470106));rng=random.Random(470107)
    for i in range(1024):
        n=8+i%25;ex=[-100,0,100][i%3];base=math.ldexp([.5,2/3,math.pi][i%3],ex);raw=int(bits(base),16);width=[2,4,16,64][i%4]
        y=[value(f'{raw+rng.randrange(-width,width+1):016x}')for _ in range(n)]
        if i%2:y=[v+math.ldexp(j,ex-48)for j,v in enumerate(y)]
        cases.append(dict(y=y,season=0,kind='fresh_tie_cloud',n=n,width=width))
    pred=[]
    for z in cases:
        f=native.fit(z['y'],pool);pred.extend([bits(f[1]),bits(f[2]),bits(2.**-52)])
    frozen('tie-fresh',dict(cases=cases,pool=pool,predictions=pred,constraints=rank['constraints'],
         method='Stable-first table order satisfying all34 retained necessary constraints. Independent eight-family controls and single-bit clouds at25 lengths, exponent shifts and binade/interior constants. No new pair or rank is added after these answers.'))
    c,s=service();d=ask(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);save('tie-fresh-answers',d)
    fail=[dict(index=i,case=i//3,prediction=p,answer=a)for i,(p,a)in enumerate(zip(pred,d['answers']))if p!=a.get('bits')]
    save('tie-fresh-result',dict(selections=len(pred),failures=fail));print('tie fresh',len(pred),'failed',len(fail),fail[:4],flush=True)


def affine_edges():
    cases=[];predictions=[];horizons=[-3.,-2.,-1.,0.,1.,2.,5.,float(2**31-1),float(2**31)]
    for ex in [-901,-900,-899,899,900,901]:
        for offset,trend in [(1.,0.),(1.,1.),(7.,-1.),(-7.,1.),(-3.,-3.),(1.,3.)]:
            y=[math.ldexp(offset+trend*j,ex)for j in range(4)]
            for season in [0,2]:
                cases.append(dict(y=y,season=season,quantum_exponent=ex,offset=offset,trend=trend))
                predictions.extend(dict(type='error',code='#NUM!')if h==2**31 else dict(type='number',bits=bits(math.ldexp(offset+trend*(3+h),ex)))for h in horizons)
    frozen('affine-edges',dict(cases=cases,horizons=horizons,predictions=predictions,
         method='Independent signed affine controls at both actual quantum limits, their immediate exponent neighbors, every historical target, and both sides of the signed32-bit horizon limit. Period0/2 pairs are frozen together. Outside the product magnitude/exactness region still counts as refused.'))
    c,s=service();d=ask(c,s,cases,'FORECAST.ETS',horizons);save('affine-edges-forecast',d)
    fail=[dict(index=i,case=i//len(horizons),prediction=p,answer=a)for i,(p,a)in enumerate(zip(predictions,d['answers']))if p!=a]
    results=dict(forecast=dict(selections=len(predictions),failures=fail))
    for name,fun,variants in [('stat','FORECAST.ETS.STAT',list(range(1,9))),('confidence','FORECAST.ETS.CONFINT',[(1.,.95),(float(2**31-1),.95),(float(2**31),.95)])]:
        ans=ask(c,s,cases,fun,variants);save('affine-edges-'+name,ans);w=len(variants)
        failed=[dict(pair=i,left=ans['answers'][2*i*w:(2*i+1)*w],right=ans['answers'][(2*i+1)*w:(2*i+2)*w])for i in range(len(cases)//2)if ans['answers'][2*i*w:(2*i+1)*w]!=ans['answers'][(2*i+1)*w:(2*i+2)*w]]
        results[name]=dict(selections=len(ans['answers']),paired_failures=failed)
    save('affine-edges-result',results);print('affine edges',len(predictions),'formula failures',len(fail),'paired',[(k,len(v.get('paired_failures',[])))for k,v in results.items()],flush=True)


if __name__=='__main__':
    steps={'search-boundaries':search_boundaries,'search-fresh':search_fresh,'short-map':short_map,'detector-map':detector_map,'short-stores':short_stores,'seasonal-map':seasonal_map,'length-map':length_map,'short-period-boundary':short_period_boundary,'length-fresh':length_fresh,'search-order-proof':search_order_proof,'short-period-fresh':short_period_fresh,'short-automatic':short_automatic,'fit-load-plateau':fit_load_plateau,'search-finish':search_finish,'search-ci':search_ci,'tie-constraints':tie_constraints,'tie-fresh':tie_fresh,'affine-edges':affine_edges}
    p=argparse.ArgumentParser();p.add_argument('step',choices=steps);a=p.parse_args();steps[a.step]()
