"""ETS recurrence identification using the existing ISA graph interpreter.

One store/precision decision per expression repeats across every observation.
Original MIT code; inputs and parameters come from retained tier-1 measurements.
"""
import ctypes as C
from itertools import product
import json
import numpy as np
from .ets import OUT, read, grouped, bits, write, sha
from .ets import service, ask_cases
from .x87_programs import Builder, Native, build_library

TAGS=['complement','prediction','error','difference','alpha_product','trend_product',
      'delta','level','beta_product','trend','gamma_product','season','mean',
      'level_center','season_center','horizon','forecast','publication']

SLOPE_TAGS=['sum_y','product_xy','sum_xy','numerator_product','numerator','quotient','pad_product','pad_sum']


def slope_graph(n, form=0, extend=False):
    """Original moment-OLS/padded-OLS candidates, with coherent store sites.

    An insufficient history can be extended to seven by its first two phases
    and the original OLS trend, then regressed again. The shape is inferred
    from impulse weights; no impulse coefficient is installed as a table.
    """
    g=Builder(n);zero=g.constant(bits(0.))
    def c(v):return g.constant(bits(float(v)))
    def op(kind,a,b,tag):
        node=g.op(kind,a,b,store=False);g.nodes[node].update(site=SLOPE_TAGS.index(tag),tag=tag);return node
    def regression(y):
        count=len(y);sx=count*(count+1)/2;sxx=count*(count+1)*(2*count+1)/6
        sy=zero;sxy=zero
        for i,v in enumerate(y):
            sy=op('add',sy,v,'sum_y');sxy=op('add',sxy,op('mul',c(i+1),v,'product_xy'),'sum_xy')
        if form==0:
            numerator=op('sub',op('mul',c(count),sxy,'numerator_product'),op('mul',c(sx),sy,'numerator_product'),'numerator');denominator=c(count*sxx-sx*sx)
        elif form==1:
            numerator=op('sub',sxy,op('div',op('mul',c(sx),sy,'numerator_product'),c(count),'numerator_product'),'numerator');denominator=c(sxx-sx*sx/count)
        elif form==2:
            numerator=op('sub',sxy,op('mul',c(sx/count),sy,'numerator_product'),'numerator');denominator=c(sxx-sx*sx/count)
        else:
            # R48: a textbook centered covariance can round before the
            # division even when raw products and totals are exact.
            mean=op('div',sy,c(count),'numerator_product')
            numerator=zero
            for i,v in enumerate(y):
                centered=op('sub',v,mean,'numerator') if form==3 else v
                numerator=op('add',numerator,op('mul',c(i+1-sx/count),centered,'product_xy'),'sum_xy')
            denominator=c(sxx-sx*sx/count)
        return op('div',numerator,denominator,'quotient')
    y=list(range(n));t=regression(y)
    if extend and n<7:
        for j in range(n,7):
            y.append(op('add',j%2,op('mul',c(j-j%2),t,'pad_product'),'pad_sum'))
        t=regression(y)
    g.sites=len(SLOPE_TAGS)
    return g.graph(t,t)

def graph(n,form=0,ef=0,phases=(),extended_initial=False,mse=False):
    g=Builder(n+6 if extended_initial else n+4);zero=g.constant(bits(0.));two=g.constant(bits(2.));gamma=g.constant(bits(2.**-52));errors=[zero]
    def op(kind,a,b,tag):
        g.phase=int(tag in phases);node=g.op(kind,a,b,store=False)
        g.nodes[node].update(site=TAGS.index(tag),tag=tag)
        return node
    h,t,a,b=n,n+1,n+2,n+3;l=0;s=[zero,zero]
    if extended_initial:
        g.phase=1;t=g.op('div',n+4,n+5,store=False);g.phase=0
    comp=op('sub',g.one,a,'complement')
    for j in range(1,n):
        q=j%2;p=op('add',l,t,'prediction')
        if ef==0:e=op('sub',j,op('add',p,s[q],'error'),'error')
        elif ef==1:e=op('sub',op('sub',j,p,'error'),s[q],'error')
        elif ef==2:e=op('sub',op('sub',j,s[q],'error'),p,'error')
        elif ef==3:e=op('sub',op('sub',op('sub',j,l,'error'),t,'error'),s[q],'error')
        else:e=op('sub',op('sub',j,l,'error'),op('add',t,s[q],'error'),'error')
        errors.append(e)
        if form==0:
            d=op('sub',j,s[q],'difference');ae=op('mul',a,d,'alpha_product');tp=op('mul',comp,p,'trend_product');l=op('add',ae,tp,'level')
        elif form in (1,2):
            d=op('sub',op('sub',j,l,'difference'),s[q],'difference');ae=op('mul',a,d,'alpha_product');tp=op('mul',comp,t,'trend_product')
            if form==1:l=op('add',l,op('add',ae,tp,'delta'),'level')
            else:l=op('add',op('add',l,ae,'delta'),tp,'level')
        else:
            ae=op('mul',a,e,'alpha_product')
            if form==3:l=op('add',l,op('add',t,ae,'delta'),'level')
            elif form==4:l=op('add',p,ae,'level')
            elif form==5:l=op('add',op('add',l,ae,'delta'),t,'level')
        t=op('add',t,op('mul',b,e,'beta_product'),'trend')
        s[q]=op('add',s[q],op('mul',gamma,e,'gamma_product'),'season')
    mean=op('div',op('add',s[0],s[1],'mean'),two,'mean')
    l=op('add',l,mean,'level_center');s=[op('sub',v,mean,'season_center')for v in s]
    prediction=op('add',l,op('mul',h,t,'horizon'),'forecast')
    outputs=[op('add',prediction,v,'publication')for v in s]
    if mse:
        acc=zero
        for e in errors:
            g.phase=0;sq=g.op('mul',e,e,store=False);g.nodes[sq]['site']=18
            acc=g.op('add',acc,sq,store=False);g.nodes[acc]['site']=19
        acc=g.op('div',acc,g.constant(bits(float(n))),store=False);g.nodes[acc]['site']=20
        outputs[1]=acc
        g.sites=21;return g.graph(*outputs)
    g.sites=len(TAGS);return g.graph(*outputs)

def run():
    path=OUT/'x87-ets-loop.so';build_library(path);native=Native(path)
    d=read(OUT/'step2-neutral-stat.json');st=grouped(d);f=grouped(read(OUT/'step2-neutral-forecast.json'))
    chosen=[i for i,z in enumerate(d['cases'])if len(z['y'])==13];inputs=[];expected=[];selector=[]
    for i in chosen:
        z=d['cases'][i]
        for h in [0.,1.,2.,5.]:
            inputs.append([int(bits(v),16)for v in z['y']+[h,z['design'][2],.002,.001]])
            expected.append(int(f[i][str(h)]['bits'],16));selector.append((12+int(h))%2)
    inputs=np.array(inputs,dtype='uint64');expected=np.array(expected,dtype='uint64');selector=np.array(selector)
    indices=np.array([12,13,14,15,28,29,30,31]);probe=np.ascontiguousarray(inputs[indices]);probe_out=np.empty((8,2),dtype='uint64');out=np.empty((len(inputs),2),dtype='uint64');ptr=C.POINTER(C.c_uint64)
    results=[];counts=[]
    for form,ef in product([0,1,2],[0,1]):
        g=graph(13,form,ef);nodes=native.nodes(json.dumps(g,sort_keys=True));present={n.get('tag')for n in g['nodes']};tags=[t for t in TAGS if t in present and t!='publication'];mask=0;previous=0;passed=0
        for k in range(1<<len(tags)):
            code=k^(k>>1);changed=code^previous
            if changed:mask^=1<<TAGS.index(tags[changed.bit_length()-1])
            previous=code
            rc=native.lib.x87_grid(nodes,len(nodes),probe.ctypes.data_as(ptr),17,8,0x133f133f,mask,*g['outputs'],probe_out.ctypes.data_as(ptr))
            if rc:raise ValueError('invalid loop graph')
            if not np.array_equal(probe_out[np.arange(8),selector[indices]],expected[indices]):continue
            passed+=1
            rc=native.lib.x87_grid(nodes,len(nodes),inputs.ctypes.data_as(ptr),17,len(inputs),0x133f133f,mask,*g['outputs'],out.ctypes.data_as(ptr))
            if rc:raise ValueError('invalid loop graph')
            if np.array_equal(out[np.arange(len(inputs)),selector],expected):results.append(dict(form=form,error_form=ef,stores=[t for t in tags if mask&(1<<TAGS.index(t))],mask=mask,graph=g))
        counts.append(dict(form=form,error_form=ef,prefilter_survivors=passed));print(counts[-1],'complete',len(results),flush=True)
        if results:break
    write(OUT/'step2-loop-search.json',dict(survivors=results,counts=counts,producer_sha256=sha(__file__)))

def proof(smooth=False):
    prefix='step2-smooth'if smooth else'step2-affine-mirror'
    native=Native(OUT/'x87-ets-loop.so');d=read(OUT/(prefix+'-stat.json'));st=grouped(d)
    profiles=read(OUT/'step2-neutral-all-lengths.json')['profiles']
    if smooth:
        allowed=set(read(OUT/'step2-affine-mirror-proof.json')['survivors']);profiles=[p for p in profiles if p['mask']in allowed]
    predictions={}
    for p in profiles:
        result=[None]*(len(d['cases'])*4)
        for n in sorted({len(z['y'])for z in d['cases']}):
            inds=[i for i,z in enumerate(d['cases'])if len(z['y'])==n];inputs=[]
            for i in inds:
                z=d['cases'][i];a=int(st[i]['1']['bits'],16);b=int(st[i]['2']['bits'],16)
                for h in [0.,1.,2.,5.]:inputs.append([int(bits(v),16)for v in z['y']+[h,z['trend']]]+[a,b])
            out=native.grid(graph(n,1),np.array(inputs,dtype='uint64'),0x133f,p['mask'])
            for j,i in enumerate(inds):
                for k,h in enumerate([0.,1.,2.,5.]):result[i*4+k]=f'{int(out[j*4+k,(n-1+int(h))%2]):016x}'
        predictions[str(p['mask'])]=result
    frozen=OUT/(prefix+'-frozen.json')
    document=dict(predictions=predictions,profiles=profiles,cases=d['cases'],producer_sha256=sha(__file__),method='All surviving shared-expression programs frozen before the new forecast questions; STAT parameters are conditional measured inputs; constructed histories have exactly represented OLS trend16.')
    if frozen.exists():
        if read(frozen)['predictions']!=predictions:raise ValueError('frozen predictions changed')
    else:write(frozen,document)
    c,s=service();answers=ask_cases(c,s,d['cases'],'FORECAST.ETS',[0.,1.,2.,5.]);write(OUT/(prefix+'-forecast.json'),answers)
    expected=[a.get('bits')for a in answers['answers']];survivors=[int(mask)for mask,pred in predictions.items()if pred==expected]
    counts={mask:sum(a==b for a,b in zip(pred,expected))for mask,pred in predictions.items()}
    write(OUT/(prefix+'-proof.json'),dict(survivors=survivors,counts=counts,cases=len(d['cases']),forecasts=len(expected),producer_sha256=sha(__file__)))
    print('fresh complete survivors',len(survivors),'forecasts',len(expected),'best diagnostic',max(counts.values()))

if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='proof':proof()
    elif len(sys.argv)>1 and sys.argv[1]=='smooth':proof(True)
    else:run()
