"""Separate the ETS state correction by adding an exact linear trend."""
from .ets import *
from .ets_proof import freeze

def correction(y,a,b,offset=0,factor='window'):
    l=y[0];t=initial_trend(y);se=[0.,0.];trends=[];n=len(y)
    for j,v in enumerate(y[1:],1):
        q=j%2;e=v-((l+t)+se[q]);l=l+(t+a*e);t=t+b*e;se[q]+=2.**-52*e;trends.append(t)
    window=n if n<=10 else min(int(.4*n),10);vs=trends[max(0,n-window-1+offset):]
    mean=sequential_total(vs)/len(vs);sd=math.sqrt(sequential_total((v-mean)**2 for v in vs)/len(vs))
    mul={'window':window/(window-1),'n':n/(n-1),'count':len(vs)/(len(vs)-1)}[factor]
    return abs(mean)<sd*mul

def separate():
    base=[10.+.05*j*j for j in range(12)];cases=[];profiles=list(itertools.product([-1,0,1],['window','n','count']));pred={str(p):[]for p in profiles}
    for k in range(-64,-15):
        shift=k/64;y=[v+shift*j for j,v in enumerate(base)];cases.append(dict(y=y,season=0,added_trend=shift))
        for p in profiles:pred[str(p)].append(bits(.001 if correction(y,.9,.899,*p)else .899))
    freeze('step5-correction-window',dict(cases=cases,predictions=pred,profiles=profiles,extra_source=Path(__file__).read_text()))
    c,s=service();d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2]);write(OUT/'step5-correction-window-answers.json',d);actual=[v['2']['bits']for v in grouped(d)];survivors=[p for p,vs in pred.items()if vs==actual]
    write(OUT/'step5-correction-window-proof.json',dict(survivors=survivors,alpha_changes=[i for i,v in enumerate(grouped(d))if value(v['1']['bits'])!=.9]));print('correction window survivors',survivors)

def edges(multiple=False):
    c,s=service();base=[0.,0.,0.,0.,0.,0.,1.,2.]
    intervals=[[-.875,-.75,False,True],[.125,.25,True,False]];trace=[];bases=[base,base];filename='step5-correction-edges.json'
    if multiple:
        d=read(OUT/'step5-correction-multiple-map.json');base_by={z['kind']:z['y']for z in d['cases']if z['added_trend']==0.}
        bases=[base_by[k]for k in ['n8','n8','n12','n12','n32','n32']]
        intervals=[[-.625,-.5625,False,True],[.125,.1875,True,False],[-1.0625,-1.,False,True],[-.6875,-.625,True,False],[-.5625,-.5,False,True],[0.,.0625,True,False]];filename='step5-correction-multiple-edges.json'
    for step in range(48):
        mids=[(lo+hi)/2 for lo,hi,_,_ in intervals]
        cases=[dict(y=[v+x*j for j,v in enumerate(base)],season=0,added_trend=x)for base,x in zip(bases,mids)]
        d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[2]);got=[value(a['bits'])==.001 for a in d['answers']]
        trace.append(dict(step=step,intervals=[list(v)for v in intervals],mids=mids,answers=d))
        for k,(mid,answer)in enumerate(zip(mids,got)):
            if answer==intervals[k][2]:intervals[k][0]=mid
            else:intervals[k][1]=mid
        write(OUT/filename,dict(intervals=intervals,trace=trace,producer_sha256=sha(__file__)))
        if step%8==0:print(step,intervals,flush=True)
    print('edges',intervals,flush=True)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='separate':separate()
    else:edges(len(sys.argv)>1 and sys.argv[1]=='multiple')
