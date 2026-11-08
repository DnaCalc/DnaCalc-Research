"""ETS confidence variance: fixed impulse shape, explicit expression rounding."""
from itertools import product
import numpy as np
from .ets import *
from .x87_programs import Builder,Native
from .ets_proof import freeze,finish

def variance(a,b,h,form=1):
    if form==39:
        f=1+(h-1)*((h*a*b+h*(2*h-1)*b*b/6)+a*a)
        k=math.floor((h-1)/2);g=2.**-52
        return f+g*k*((2*a+g)+2*b*(k+1))
    if form==1:f=1+(h-1)*(a*a+a*b*h+b*b*h*(2*h-1)/6)
    elif form==4:f=1+(h-1)*(a*a+b*(a*h+b*h*(2*h-1)/6))
    else:f=1+(h-1)*(a*a+b*h*(a+b*(2*h-1)/6))
    k=math.floor((h-1)/2);g=2.**-52
    return f+g*k*(2*a+g+2*b*(k+1))

POLY_TAGS=['aa','ab','bb','sum','main','season_inner','season_product','factor','scale','root','publication']
def polynomial_graph(h,form):
    g=Builder(4);sigma,a,b,z=range(4);zero=g.constant(bits(0.));eps=g.constant(bits(2.**-52))
    H=g.constant(bits(float(h)));HM=g.constant(bits(float(h-1)));TWOH=g.constant(bits(float(2*h-1)));six=g.constant(bits(6.));two=g.constant(bits(2.));k=g.constant(bits(float((h-1)//2)));kp=g.constant(bits(float((h-1)//2+1)))
    def op(kind,x,y,tag):
        n=g.op(kind,x,y,store=False);g.nodes[n].update(site=POLY_TAGS.index(tag),tag=tag);return n
    aa=op('mul',a,a,'aa')
    if form==1:
        ab=op('mul',op('mul',a,b,'ab'),H,'ab');bb=op('div',op('mul',op('mul',op('mul',b,b,'bb'),H,'bb'),TWOH,'bb'),six,'bb');inner=op('add',op('add',aa,ab,'sum'),bb,'sum')
    elif form==4:
        bh=op('mul',a,H,'ab');bb=op('div',op('mul',op('mul',b,H,'bb'),TWOH,'bb'),six,'bb');inner=op('add',aa,op('mul',b,op('add',bh,bb,'sum'),'ab'),'sum')
    else:
        bb=op('div',op('mul',b,TWOH,'bb'),six,'bb');inner=op('add',aa,op('mul',op('mul',b,H,'ab'),op('add',a,bb,'sum'),'ab'),'sum')
    f=op('add',g.one,op('mul',HM,inner,'main'),'main')
    si=op('add',op('add',op('mul',two,a,'season_inner'),eps,'season_inner'),op('mul',op('mul',two,b,'season_inner'),kp,'season_inner'),'season_inner');sg=op('mul',op('mul',eps,k,'season_product'),si,'season_product');f=op('add',f,sg,'factor')
    root=op('sqrt',op('mul',sigma,f,'scale'),zero,'root');out=op('mul',z,root,'publication');g.sites=len(POLY_TAGS);return g.graph(out,out)

def arithmetic_search():
    native=Native(OUT/'x87-ets-loop.so');d=read(OUT/'step2-final-recurrence-parameters.json');st=grouped(d);fo=grouped(read(OUT/'step6-confidence-polynomial-answers.json'));z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits']);rows=[]
    for cs,s in zip(d['cases'][:16],st):
        a=value(s['1']['bits']);b=value(s['2']['bits']);_,_,_,e,_=numeric_fit(cs['y'],a,b,cs['trend']);sigma=sequential_total(v*v for v in e)/len(e);rows.append([int(bits(v),16)for v in [sigma,a,b,z]])
    rows=np.array(rows,dtype='uint64');expected={h:np.array([int(v[str([float(h),.95])]['bits'],16)for v in fo],dtype='uint64')for h in range(1,11)};survivors=[]
    for form in [1,4,5]:
        graphs={h:polynomial_graph(h,form)for h in range(1,11)}
        for mask in range(1<<len(POLY_TAGS)):
            if all(np.array_equal(native.grid(graphs[h],rows,0x133f,mask)[:,0],expected[h])for h in [3,7,10,1,2,4,5,6,8,9]):survivors.append(dict(form=form,mask=mask,stores=[v for i,v in enumerate(POLY_TAGS)if mask&(1<<i)]))
    write(OUT/'step6-polynomial-arithmetic.json',dict(survivors=survivors,tags=POLY_TAGS,producer_sha256=sha(__file__)));print('polynomial arithmetic survivors',len(survivors),survivors[:4])

def product_expressions(factors):
    """All multiplication associations of one fixed factor bag, deduplicated.

    The inverse-six leaf permits either a literal reciprocal or division.
    This changes expression order, never polynomial degree or sample partitions.
    """
    cache={}
    def visit(mask):
        if mask in cache:return cache[mask]
        found={}
        if mask&(mask-1)==0:
            i=mask.bit_length()-1;a,e=factors[i];found[a.tobytes()]=(a,e)
        else:
            sub=(mask-1)&mask
            while sub:
                other=mask^sub
                if sub<other:
                    for av,ae in visit(sub).values():
                        for bv,be in visit(other).values():
                            choices=[(av*bv,'('+ae+'*'+be+')')]
                            if ae=='(1/6)':choices.append((bv/6,'('+be+'/6)'))
                            if be=='(1/6)':choices.append((av/6,'('+ae+'/6)'))
                            for v,e in choices:
                                key=v.tobytes()
                                if key not in found or len(e)<len(found[key][1]):found[key]=(v,e)
                sub=(sub-1)&mask
        cache[mask]=found;return found
    return list(visit((1<<len(factors))-1).values())

def finite_trees():
    d=read(OUT/'step2-final-recurrence-parameters.json');st=grouped(d);f=read(OUT/'step6-confidence-polynomial-answers.json');z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits']);rows=[]
    for cs,s in zip(d['cases'][:16],st):
        a=value(s['1']['bits']);b=value(s['2']['bits']);_,_,_,e,_=numeric_fit(cs['y'],a,b,cs['trend']);sigma=sequential_total(v*v for v in e)/len(e)
        for h in range(1,11):rows.append([a,b,float(h),sigma])
    a,b,h,sigma=np.array(rows).T;k=np.floor((h-1)/2);g=2.**-52;actual=np.array([int(v['bits'],16)for v in f['answers']],dtype='uint64')
    aa=a*a;ab=product_expressions([(a,'a'),(b,'b'),(h,'h')]);bb=product_expressions([(b,'b'),(b,'b'),(h,'h'),(2*h-1,'(2*h-1)'),(np.full_like(h,1/6),'(1/6)')]);sg=[]
    for order in itertools.permutations(range(3)):
        terms=[(2*a,'2*a'),(np.full_like(a,g),'g'),(2*b*(k+1),'2*b*(k+1)')];inner=(terms[order[0]][0]+terms[order[1]][0])+terms[order[2]][0];expr='(('+terms[order[0]][1]+'+'+terms[order[1]][1]+')+'+terms[order[2]][1]+')'
        sg.extend([(g*k*inner,'g*k*'+expr),(g*(k*inner),'g*(k*'+expr+')')])
    best=0;survivors=[]
    print('finite terms',len(ab),len(bb),len(sg),flush=True)
    for (abv,abe),(bbv,bbe),order,(sgv,sge),add,scale in itertools.product(ab,bb,list(itertools.permutations(range(3))),sg,range(2),range(2)):
        terms=[(aa,'a*a'),(abv,abe),(bbv,bbe)];inner=(terms[order[0]][0]+terms[order[1]][0])+terms[order[2]][0];outer=(h-1)*inner;factor=(1+outer)+sgv if not add else 1+(outer+sgv);var=sigma*factor if not scale else sigma*(1+outer)+sigma*sgv
        pred=z*np.sqrt(var);same=int(np.count_nonzero(pred.view('uint64')==actual));best=max(best,same)
        if same==len(actual):survivors.append(dict(ab=abe,bb=bbe,order=order,season=sge,add=add,scale=scale))
    write(OUT/'step6-finite-trees.json',dict(survivors=survivors,best=best,producer_sha256=sha(__file__)));print('finite tree survivors',len(survivors),'best diagnostic',best,survivors[:3])

def final_proof():
    name='step6-final-confidence';d=read(OUT/'step2-final-recurrence-parameters.json');st=grouped(d);cases=d['cases'][16:48];z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits']);pred=[]
    for cs,s in zip(cases,st[16:48]):
        a=value(s['1']['bits']);b=value(s['2']['bits']);_,_,_,e,_=numeric_fit(cs['y'],a,b,cs['trend']);sigma=sequential_total(v*v for v in e)/len(e)
        for h in range(1,11):pred.append(bits(z*math.sqrt(sigma*variance(a,b,float(h),39))))
    freeze(name,dict(cases=cases,predictions=pred,conditional_inputs='Measured STAT smoothing leaves; exact constructed slope; identified closed variance expression order.',extra_source=Path(__file__).read_text()))
    finish(name,cases,pred,'FORECAST.ETS.CONFINT',[(float(h),.95)for h in range(1,11)])

TAGS=['multiple','impulse','season','square','accumulate','factor','scale','root','publication']
def graph(h,direct=True,seed=True,season_order=0):
    g=Builder(4);sigma,a,b,z=range(4);zero=g.constant(bits(0.));eps=g.constant(bits(2.**-52))
    def op(kind,x,y,tag):
        k=g.op(kind,x,y,store=False);g.nodes[k].update(site=TAGS.index(tag),tag=tag);return k
    acc=g.one if seed else zero;imp=a
    for j in range(1,h):
        bj=op('mul',b,g.constant(bits(float(j))),'multiple')
        imp=op('add',a,bj,'impulse')if direct else op('add',imp,b,'impulse')
        w=imp
        if j%2==0:
            if season_order==0:w=op('add',imp,eps,'season')
            else:w=op('add',op('add',a,eps,'season'),bj,'impulse')
        acc=op('add',acc,op('mul',w,w,'square'),'accumulate')
    if not seed:acc=op('add',g.one,acc,'factor')
    v=op('mul',sigma,acc,'scale');r=op('sqrt',v,zero,'root');out=op('mul',z,r,'publication');g.sites=len(TAGS);return g.graph(out,out)

def discovery():
    native=Native(OUT/'x87-ets-loop.so');d=read(OUT/'step6-initial-confidence.json');st=grouped(read(OUT/'step5-staircase-fullstat.json'));fo=grouped(d)
    rows={h:[]for h in [1,2,3]};ex={h:[]for h in [1,2,3]};z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits'])
    for i,cs in enumerate(d['cases'][:8]):
        a=value(st[i]['1']['bits']);b=value(st[i]['2']['bits']);_,_,_,e,_=numeric_fit(cs['y'],a,b,initial_trend(cs['y']));sigma=sequential_total(v*v for v in e)/len(e)
        for h in rows:
            rows[h].append([int(bits(v),16)for v in [sigma,a,b,z]]);ex[h].append(int(fo[i][str([float(h),.95])]['bits'],16))
    survivors=[]
    for direct,seed,order in product([True,False],[True,False],[0,1]):
        for mask in range(1<<len(TAGS)):
            if all(np.array_equal(native.grid(graph(h,direct,seed,order),np.array(rows[h],dtype='uint64'),0x133f,mask)[:,0],ex[h])for h in rows):survivors.append(dict(direct=direct,seed=seed,season_order=order,mask=mask))
    write(OUT/'step6-variance-programs.json',dict(programs=survivors,tags=TAGS,producer_sha256=sha(__file__)))
    print('confidence programs',len(survivors),survivors[:4])

def proof():
    name='step6-confidence-polynomial';d=read(OUT/'step2-final-recurrence-parameters.json');st=grouped(d);cases=d['cases'][:16]
    programs=[1,4,5];allpred={};z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits']);rows=[]
    for cs,s in zip(cases,st):
        a=value(s['1']['bits']);b=value(s['2']['bits']);_,_,_,e,_=numeric_fit(cs['y'],a,b,cs['trend']);sigma=sequential_total(v*v for v in e)/len(e);rows.append([sigma,a,b,z])
    for form in programs:
        allpred[str(form)]=[bits(z*math.sqrt(sigma*variance(a,b,float(h),form)))for sigma,a,b,z in rows for h in range(1,11)]
    freeze(name,dict(cases=cases,predictions=allpred,programs=programs,conditional_inputs='STAT smoothing leaves and exact constructed initial slopes; full error SSE/n; live NORM.S.INV(.975).',extra_source=Path(__file__).read_text()))
    c,s=service();d=ask_cases(c,s,cases,'FORECAST.ETS.CONFINT',[(float(h),.95)for h in range(1,11)]);write(OUT/(name+'-answers.json'),d);actual=[a.get('bits')for a in d['answers']]
    surviving=[k for k,p in allpred.items()if p==actual];counts={k:sum(a==b for a,b in zip(p,actual))for k,p in allpred.items()};write(OUT/(name+'-proof.json'),dict(survivors=surviving,counts=counts,selections=len(actual)));print('fresh variance survivors',len(surviving),'best diagnostic',max(counts.values())if counts else 0)

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='proof':proof()
    elif len(sys.argv)>1 and sys.argv[1]=='arithmetic':arithmetic_search()
    elif len(sys.argv)>1 and sys.argv[1]=='trees':finite_trees()
    elif len(sys.argv)>1 and sys.argv[1]=='final':final_proof()
    else:discovery()
