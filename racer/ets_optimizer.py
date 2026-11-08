"""Finite objective and optimiser candidates, with retained counterexamples."""
import itertools
from .ets import *

# The measured point family is frozen before the next random histories. Pairs
# introduced solely by the later trend correction are not optimisation points.
GRID=[(.9,.001),(.9,.899),(.833333,.75),(.75,.001),(.75,.749),
      (.666667,.166667),(.666667,.5),(.5,.001),(.5,.25),(.5,.499),
      (.333333,.166667),(.25,.001),(.25,.249),(.166667,.0833333),
      (.126,.001),(.126,.125),(.1,.001),(.1,.099),(.002,.001),
      (.998,.099),(.998,.499),(.998,.749),(.998,.899)]

def parameters(y):
    return min(((objective(y,a,b),a,b)for a,b in GRID),key=lambda z:z[0])[1:]

def corrected(y,a,b):
    l=y[0];t=initial_trend(y);se=[0.,0.];trends=[];errors=[0.];predictions=[y[0]];n=len(y)
    for j,v in enumerate(y[1:],1):
        q=j%2;p=(l+t)+se[q];e=v-p;l=l+(t+a*e);t=t+b*e;se[q]+=2.**-52*e;trends.append(t);errors.append(e);predictions.append(p)
    window=n if n<=10 else (n+1)//2;vs=trends[-(window-1):]
    mean=sequential_total(vs)/len(vs);sd=math.sqrt(sequential_total((v-mean)**2 for v in vs)/len(vs))
    correction=abs(mean)<sd*(len(vs)+1)/len(vs)
    if correction:l=y[-1]-se[(n-1)%2];t=initial_trend(y);b=.001
    center=(se[0]+se[1])/2;l+=center;se=[v-center for v in se]
    return l,t,se,b,errors,predictions,correction

def fresh():
    from .ets_proof import freeze
    name='step45-final-100';rng=random.Random(395919);cases=[];pred=[];raw=[]
    for i in range(100):
        n=rng.randrange(7,33);ex=rng.randrange(-120,121);trend=rng.randrange(-32,33)/8
        if i<75:y=[math.ldexp(rng.randrange(-256,257)/16+trend*j,ex)for j in range(n)]
        else:y=[math.ldexp((j*j)/32+trend*j+rng.randrange(-8,9)/64,ex)for j in range(n)]
        cases.append(dict(y=y,season=0));a,b=parameters(y);l,t,se,pb,e,p,corr=corrected(y,a,b)
        pred.extend([bits(a),bits(pb)]);raw.append(dict(alpha=a,beta=b,correction=corr))
    freeze(name,dict(cases=cases,predictions=pred,raw_parameters=raw,grid=GRID,extra_source=Path(__file__).read_text(),method='One fixed finite pair family, hybrid SSE, fixed trend correction expression; all random histories frozen before the questions.'))
    c,s=service();d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',range(1,9));write(OUT/(name+'-answers.json'),d)
    st=grouped(d);bad=[]
    for i in range(100):
        actual=[st[i][str(j)]['bits']for j in [1,2]]
        if actual!=pred[i*2:i*2+2]:bad.append(dict(case=i,predicted=pred[i*2:i*2+2],actual=actual,raw=raw[i]))
    write(OUT/(name+'-proof.json'),dict(failures=bad,cases=100,producer_sha256=sha(__file__)));print('fresh parameters failures',len(bad),bad[:10],flush=True)

def objective(y,a,b,kind='hybrid',omit=0):
    n=len(y);tail=n if n<=10 else min(int(.4*n),10);start=n-tail
    l=y[0];t=initial_trend(y);se=[0.,0.];errors=[0.]
    for j,v in enumerate(y[1:],1):
        q=j%2;p=(l+t)+se[q];e=v-p;errors.append(e)
        if start and j>=start and kind.startswith('hybrid'):l=l+t
        else:l=l+(t+a*e);t=t+b*e;se[q]=se[q]+2.**-52*e
    begin=start if kind.endswith('tail')else 0
    return sequential_total(v*v for v in errors[begin:n-omit or n])

def documents():
    for name in ['step2-stat.json','step3-prefix-stat.json','step5-mapping-stat.json','step2-smooth-stat.json','step45-fresh-100.json']:
        d=read(OUT/name);st=grouped(d)
        for i,(z,s)in enumerate(zip(d['cases'],st)):
            if len(z['y'])<7:continue
            yield name,i,z,value(s['1']['bits']),value(s['2']['bits'])

def analyze():
    ds=list(documents());pool=sorted(set((a,b)for _,_,_,a,b in ds));result=[]
    for kind,omit in itertools.product(['hybrid','hybrid_tail','full','full_tail'],[0,1]):
        bad=[]
        for name,i,z,a,b in ds:
            actual=objective(z['y'],a,b,kind,omit);score,x,y=min((objective(z['y'],x,y,kind,omit),x,y)for x,y in pool)
            if score<actual*(1-1e-12):bad.append(dict(source=name,case=i,actual=[a,b,actual],separating=[x,y,score]))
        result.append(dict(kind=kind,omit_last=omit,refutations=bad));print(kind,omit,'refutations',len(bad),'examples',bad[:3])
    write(OUT/'step4-objective-exact-tree.json',dict(pool=pool,results=result,producer_sha256=sha(__file__)))

if __name__=='__main__':
    import sys
    fresh()if len(sys.argv)>1 and sys.argv[1]=='fresh'else analyze()
