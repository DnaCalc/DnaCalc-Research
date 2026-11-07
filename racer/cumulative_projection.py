"""R11 necessary singleton constraints with one arbitrary Ext80 principal.

Extends the existing exact preimage tools; does not fit or install a principal.
Only end-timing singleton CUMPRINC observations are used in this projection.
"""
from collections import defaultdict, Counter
from fractions import Fraction as Q
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time
import datetime
import struct

from intervals import B64, E64, E53, Interval, inverse_op, inverse_add, rounding_inverse, log2floor, number, pow2

MODELS=('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53')


def direct_graph(function, count, beginning, first, logarithm, association, discount):
    """Ordinary closed-form loan algebra; coefficients are only0,1 and input counts.

    End principal at k is -r*(PV+FV)*q^(k-1)/(q^N-1).
    Summing its geometric progression gives CUMPRINC without calling PMT/PPMT.
    Beginning payments move the power by one, with the first payment separate.
    """
    from x87_programs import Builder
    g=Builder(count);one=g.one
    def op(kind,a,b=0,site=None):
        i=g.op(kind,a,b,store=False)
        if site is not None:g.nodes[i]['site']=site
        return i
    def add(a,b,site=6):return op('add',a,b,site)
    def sub(a,b,site=6):return op('sub',a,b,site)
    def mul(a,b,site=8):return op('mul',a,b,site)
    def div(a,b,site=8):return op('div',a,b,site)
    def neg(a):return op('neg',a)
    base=add(one,0,0)
    log=op('fyl2xp1',one,0,1)if logarithm==0 else op('fyl2x',one,base,1)
    def power(exponent):
        t=mul(exponent,log,2)
        integer=op('frndint',t)
        fraction=sub(t,integer,3)
        mantissa=add(one,op('f2xm1',fraction,site=3),3)
        return op('fscale',mantissa,integer,4)
    def combine(a,b,denominator):
        return div(mul(a,b,7),denominator,8)if association==0 else mul(a,div(b,denominator,8),7)
    if function=='PPMT':
        period,nper,pv,fv=1,2,3,4
    else:
        nper,pv=1,2
        fv=g.constant('0000000000000000')
    total_power=power(neg(nper)if discount else nper)
    denominator=sub(one,total_power,5)if discount else sub(total_power,one,5)
    def payment():
        numerator=add(pv,mul(fv,total_power,7),7)if discount else add(mul(pv,total_power,7),fv,7)
        value=neg(combine(0,numerator,denominator))
        return div(value,base,9)if beginning else value
    if function=='PPMT':
        if beginning and first:
            result=payment()
        else:
            degree=add(sub(nper,period),one)if discount else sub(period,one)
            if beginning:degree=add(degree,one)if discount else sub(degree,one)
            factor=power(neg(degree)if discount else degree)
            result=neg(combine(mul(0,add(pv,fv),7),factor,denominator))
    else:
        start,end=3,4
        if beginning and first:
            first_payment=payment()
            length=sub(end,one)
            degree=add(sub(nper,end),one)if discount else g.constant('0000000000000000')
        else:
            first_payment=None
            length=add(sub(end,start),one)
            degree=sub(nper,end)if discount else sub(start,one)
            if beginning:degree=add(degree,one)if discount else sub(degree,one)
        span=power(neg(length)if discount else length)
        difference=sub(one,span,5)if discount else sub(span,one,5)
        factor=power(neg(degree)if discount else degree)
        result=neg(combine(pv,mul(factor,difference,7),denominator))
        if first_payment is not None:result=add(first_payment,result,9)
    g.sites=10
    graph=g.graph(result,result)
    if len(graph['nodes'])>128:raise ValueError('closed-form graph exceeds native interpreter bound')
    return graph


def direct_run(input_path,library,destination):
    """Reject complete finite programs on exact observed outputs, never pick per row."""
    from x87_programs import Native,WORDS,digest
    dest=Path(destination);dest.mkdir()
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    frozen={str(p):digest(p)for p in (Path(input_path),Path(library),Path(__file__),Path(__file__).with_name('x87_programs.py'),Path(__file__).with_name('x87_native.c'))}
    native=Native(library);counts={};graphs={};scope=Counter();excluded=[]
    profiles=[dict(id=i,logarithm=l,association=a,discount=d,cw=cw,mask=mask)
        for i,(l,a,d,cw,mask)in enumerate(itertools.product(range(2),range(2),(False,True),WORDS,range(1024)))]
    with(dest/'profiles.json').open('x')as out:json.dump(profiles,out,separators=(',',':'));out.write('\n')
    for function in ('CUMPRINC','PPMT'):
        selected=[]
        for index,row in enumerate(rows):
            if row['function']!=function:continue
            args=row['args']
            try:
                if row['result_index']is not None or not row['expected'].startswith('0x')or not all(isinstance(v,str)and v.startswith('0x')for v in args):raise ValueError('surface_or_result')
                if function=='PPMT':
                    if not 4<=len(args)<=6:raise ValueError('arity')
                    args=args+['0x0000000000000000']*(6-len(args))
                    rate,period,nper,pv,fv,t=map(number,args)
                    if period!=int(period)or not 1<=period<=nper:raise ValueError('period')
                    first=period==1
                else:
                    if len(args)!=6:raise ValueError('arity')
                    rate,nper,pv,start,end,t=map(number,args)
                    if not(pv>0 and start==int(start)and end==int(end)and 1<=start<=end<=nper):raise ValueError('interval')
                    first=start==1
                if not(0<rate<=1 and nper==int(nper)and 1<=nper<=1000 and t in (0,1)):raise ValueError('rate_count_timing')
            except (ValueError,TypeError,OverflowError)as error:
                scope[function+':'+str(error)]+=1;excluded.append(dict(trace_index=index,reason=str(error)));continue
            selected.append(dict(index=index,args=args,first=first,beginning=t==1))
        if not selected:raise ValueError('empty closed-form scope '+function)
        alive={p['id']for p in profiles};executions=0
        with(dest/(function+'-first-kills.jsonl')).open('x')as out:
            for order,row in enumerate(selected):
                if not alive:break
                raw=rows[row['index']];inputs=list(map(lambda v:int(v,16),row['args']))
                for l,a,d,cw in itertools.product(range(2),range(2),(False,True),WORDS):
                    wanted=[p for p in profiles if p['id']in alive and(p['logarithm'],p['association'],p['discount'],p['cw'])==(l,a,d,cw)]
                    if not wanted:continue
                    key=(function,row['beginning'],row['first'],l,a,d)
                    if key not in graphs:graphs[key]=direct_graph(function,6,row['beginning'],row['first'],l,a,d)
                    values=native.masks(graphs[key],inputs,cw)
                    for p in wanted:
                        offset=20*p['mask'];value=values[offset:offset+10]
                        got=f'0x{native.publish(value,cw):016x}';executions+=1
                        if got!=raw['expected']:
                            alive.remove(p['id'])
                            out.write(json.dumps(dict(candidate_id=p['id'],trace_index=row['index'],order=order,
                                expected=raw['expected'],got=got,root_ext80_le=value.hex()))+'\n')
        with(dest/(function+'-survivors.json')).open('x')as out:json.dump(sorted(alive),out);out.write('\n')
        counts[function]=dict(programs=len(profiles),available_scope=len(selected),executions=executions,survivors=len(alive))
    with(dest/'graphs.json').open('x')as out:json.dump([dict(key=k,graph=g)for k,g in graphs.items()],out,separators=(',',':'));out.write('\n')
    if frozen!={p:digest(p)for p in frozen}:raise ValueError('closed-form source/input changed')
    result=dict(complete=True,counts=counts,scope_exclusions=scope,excluded=excluded,input_sha256=frozen,
        grammar='16384 labelled complete programs per function: two logarithm instructions, two numerator/division associations, growth/discount algebra, PC53/PC64 and1024 coherent store masks. Labels can repeat physical programs.',
        limits='All observations are discovery. Native x87 leaves are measured on this CPU, not a Windows claim. '
               'A first contradiction excludes its fixed program before other rows are needed; no accumulation order, internal payment, full kernel or product gain is inferred.')
    with(dest/'summary.json').open('x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def summation_questions(checkout, destination, cache, offline=False, cohort='round34'):
    """Read every subrange of controlled loans, independently of a growth model."""
    import math
    try:
        from .fingerprint import load_service, CachedService
    except ImportError:
        from fingerprint import load_service, CachedService
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    destination = Path(destination)
    destination.mkdir(exist_ok=False)
    loans, requests = [], []
    if cohort not in ('round34','validation35','dyadic35'):raise ValueError('unknown summation cohort')
    if cohort=='dyadic35':
        # In exact arithmetic these loans have PMT=-base**periods (end),
        # and per-period principal -base**(period-1). Public answers remain
        # the observations; this identity is only an operand-reading control.
        for rate,periods,timing in itertools.product((1.,3.,7.),(8,16,24,40,52,53,54,60,64,80),(0,1)):
            base=1+int(rate);pv=float(Q(base**periods-1,int(rate)))
            loans.append(dict(id=len(loans),rate=rate,periods=periods,pv=pv,timing=timing))
    else:
        rates=(2.**-10,.05,.5,1.)if cohort=='round34'else(.007,.03125,.1875,.8)
        scales=(-40,0,40)if cohort=='round34'else(-32,0,32)
        for rate,scale,offset,timing in itertools.product(rates,scales,(-1,0,1),(0,1)):
            pv=math.ldexp(1000.,scale)
            if offset:pv=math.nextafter(pv,math.inf if offset>0 else 0.)
            loans.append(dict(id=len(loans),rate=rate,periods=8,pv=pv,timing=timing,scale=scale,offset=offset))
    for loan in loans:
        periods=loan['periods']
        if cohort=='dyadic35':
            ranges={(p,p)for p in range(1,periods+1)}
            for start in (1,2,max(1,periods-20),periods-2,periods):
                for end in (start,min(periods,start+1),periods):ranges.add((start,end))
        else:ranges={(start,end)for start in range(1,periods+1)for end in range(start,periods+1)}
        for start,end in sorted(ranges):
            requests.append(dict(loan=loan['id'],start=start,end=end,args=[client.number(v)for v in
                (loan['rate'],float(periods),loan['pv'],float(start),float(end),float(loan['timing']))]))
    inputs = dict(loans=loans, ranges=requests)
    (destination/'questions.json').write_text(json.dumps(inputs, separators=(',', ':'))+'\n')
    results = {}
    for function in ('CUMPRINC', 'CUMIPMT'):
        results[function] = service.ask_rows(function, [r['args'] for r in requests], chunk_size=2000)
    results['PMT'] = service.ask_rows('PMT', [[client.number(v) for v in
        (l['rate'], float(l['periods']), l['pv'], 0., float(l['timing']))] for l in loans])
    (destination/'answers.json').write_text(json.dumps(results, separators=(',', ':'))+'\n')
    result = dict(complete=True, cohort=cohort, loans=len(loans), ranges=len(requests),
        calls=sum(map(len, results.values())), receipts=service.receipts,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limits='Every observed subrange has its own original Excel answer. No growth, payment or numerator program is assumed; subrange observations can refute shared-term summation premises.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))



def validate_sum_observations(questions,answers):
    """A missing range or operand must not turn an empty test into a survivor."""
    loans=questions['loans'];ranges=questions['ranges']
    if not loans or not ranges:raise ValueError('nonempty cumulative loans and ranges required')
    if [l['id']for l in loans]!=list(range(len(loans))):raise ValueError('ordered cumulative loan IDs required')
    if len(answers['PMT'])!=len(loans)or any(len(answers[f])!=len(ranges)for f in ('CUMPRINC','CUMIPMT')):
        raise ValueError('unaligned cumulative observations')
    keys=[(r['loan'],r['start'],r['end'])for r in ranges]
    if len(set(keys))!=len(keys)or any(not 0<=r['loan']<len(loans)or not 1<=r['start']<=r['end']<=loans[r['loan']]['periods']for r in ranges):
        raise ValueError('invalid or duplicate cumulative range')

def paired_finish(source,destination):
    """Read CUMIPMT's finish from independently observed cumulative principal.

    The beginning first period is absent from both the payment count and the
    principal operand. This does not assume any principal-producing program.
    """
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checker'))
    from exact import operation as arithmetic
    source,destination=map(Path,(source,destination));destination.mkdir(exist_ok=False)
    q=json.loads((source/'questions.json').read_text());a=json.loads((source/'answers.json').read_text())
    validate_sum_observations(q,a)
    def val(answer):
        if answer['type']!='number':raise ValueError('paired finish requires numeric controls')
        return Q(number('0x'+answer['bits']))
    principal={(r['loan'],r['start'],r['end']):p for r,p in zip(q['ranges'],a['CUMPRINC'])}
    programs=[];survivors=[]
    for product,finish in itertools.product(MODELS,MODELS):
        first=None;checked=0;matches=0
        for index,(r,answer)in enumerate(zip(q['ranges'],a['CUMIPMT'])):
            loan=q['loans'][r['loan']]
            start=max(2,r['start'])if loan['timing']else r['start'];end=r['end']
            count=max(0,end-start+1)
            p=val(principal[r['loan'],start,end])if count else Q(0)
            total=arithmetic(Q(count),val(a['PMT'][r['loan']]),'mul',product)
            got=B64.rn(arithmetic(total,-p,'add',finish));wanted=val(answer);checked+=1;matches+=got==wanted
            if got!=wanted and first is None:
                first=dict(index=index,original=r,principal_range=[start,end]if count else None,
                    payment=a['PMT'][r['loan']],principal='0x'+struct.pack('>d',float(p)).hex(),
                    expected=answer,prediction='0x'+struct.pack('>d',float(got)).hex())
        identity=len(programs)
        if first is None:survivors.append(identity)
        programs.append(dict(id=identity,product_model=product,finish_model=finish,rows=checked,matched=matches,first_contradiction=first))
    result=dict(complete=True,programs=programs,survivors=survivors,
        input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in(source/'questions.json',source/'answers.json',Path(__file__),Path(__file__).parents[1]/'checker/exact.py',Path(__file__).with_name('intervals.py'))},
        limits='The principal operand is the independent public CUMPRINC answer for the adjusted interval, not a sum of PPMT and not a fitted private operand. Every declared arithmetic word is compared exactly. Multiple survivors stay indistinguishable on this cohort; the cumulative principal producer remains open.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(rows=len(q['ranges']),survivors=survivors)))


def long_finish_questions(checkout,output,cache,cohort='validation35'):
    """Use integer-count product ties to discriminate the paired interest finish."""
    import math
    from fingerprint import CachedService,load_service,bits,write,sha
    from exact import operation as arithmetic,binary64
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=180),cache)
    if cohort not in ('discovery35','validation35'):raise ValueError('unknown finish cohort')
    rates=(2.**-14,.0001,.0003)if cohort=='discovery35'else(2.**-13,.00013,.00027)
    balances=(1.,math.nextafter(1000.,math.inf),1234.56789,math.nextafter(1.,0.))if cohort=='discovery35'else(1.3,math.nextafter(8000.,math.inf),500.25,.75)
    n=16384. if cohort=='discovery35'else 12001.
    loans=[dict(rate=r,periods=n,pv=pv,timing=t)for r,pv,t in itertools.product(rates,balances,(0.,1.))]
    payments=service.ask_rows('PMT',[[client.number(l[k])for k in('rate','periods','pv')]+[client.number(0.),client.number(l['timing'])]for l in loans])
    def round_integer(p,precision):
        shift=p.bit_length()-precision
        if shift<=0:return p
        q,r=divmod(p,1<<shift);half=1<<(shift-1)
        return(q+(r>half or r==half and q&1))<<shift
    rows=[]
    for loan,payment in zip(loans,payments):
        if payment['type']!='number':raise ValueError('finite PMT required')
        raw=int(payment['bits'],16);mantissa=(raw&((1<<52)-1))+(1<<52);selected=[]
        for count in range(2049,8192):
            product=mantissa*count
            if round_integer(product,53)!=round_integer(round_integer(product,64),53):selected.append(count)
            if len(selected)==3:break
        for count in selected:
            start=2 if loan['timing']else 1;end=start+count-1
            rows.append(dict(loan=loan,payment=payment,count=count,start=start,end=end,
                args=[client.number(v)for v in(loan['rate'],n,loan['pv'],float(start),float(end),loan['timing'])]))
    if len(payments)!=len(loans):raise ValueError('unaligned payment observations')
    if not rows:raise ValueError('no separating cumulative count-product ties')
    source=sha(__file__)
    write(str(output)+'.questions.json',dict(format='dnacalc-cumulative-count-ties/1',cohort=cohort,rows=rows,
        receipts=service.receipts.copy(),source_sha256=source,
        selection='Only the integer product of count and independently observed PMT selects the questions. No cumulative answer inspected.'))
    principals=service.ask_rows('CUMPRINC',[r['args']for r in rows],chunk_size=100)
    profiles=list(itertools.product(MODELS,MODELS));predictions=[]
    for row,p in zip(rows,principals):
        if p['type']!='number':raise ValueError('finite cumulative principal required')
        outcomes=[]
        for product,finish in profiles:
            total=arithmetic(row['count'],binary64(row['payment']['bits']),'mul',product)
            result=arithmetic(total,-binary64(p['bits']),'add',finish)
            outcomes.append(bits(float(B64.rn(result))))
        predictions.append(dict(row,principal=p,predictions=outcomes))
    write(str(output)+'.predictions.json',dict(format='dnacalc-cumulative-count-tie-predictions/1',
        records=predictions,profiles=profiles,receipts=service.receipts.copy(),source_sha256=source))
    answers=service.ask_rows('CUMIPMT',[r['args']for r in rows],chunk_size=100)
    counts=[0]*len(profiles)
    for row,answer in zip(predictions,answers):
        if answer['type']!='number':raise ValueError('finite cumulative interest required')
        row['interest']=answer
        for i,prediction in enumerate(row['predictions']):counts[i]+=prediction==answer['bits']
    if sha(__file__)!=source:raise ValueError('finish reader changed during observation')
    result=dict(format='dnacalc-cumulative-count-tie-finish/1',complete=True,cohort=cohort,
        records=predictions,profiles=profiles,matches=counts,survivors=[i for i,c in enumerate(counts)if c==len(rows)],
        receipts=service.receipts,source_sha256=source,
        limits='Conditional public PMT and adjusted public CUMPRINC operands, with25 fixed product/subtraction words. All predictions frozen before the new CUMIPMT answers. Remaining equal words are undetermined on this cohort. No principal-producing or summation program is selected.')
    write(output,result);print(json.dumps(dict(rows=len(rows),survivors=result['survivors'],matches=counts)))


def balance_sum_questions(checkout,source,output,cache,include_fv=False):
    """Condition backward-balance sum programs on independent public PV/IPMT."""
    from fingerprint import CachedService,load_service,bits,write,sha
    from exact import operation as arithmetic,binary64
    source,output=map(Path,(source,output));output.mkdir(exist_ok=False)
    q=json.loads((source/'questions.json').read_text());service_client,_=load_service(checkout)
    service=CachedService(service_client.Client(timeout=180),cache)
    answers=json.loads((source/'answers.json').read_text());validate_sum_observations(q,answers);loans=q['loans'];rows=[]
    for loan in loans:
        payment=answers['PMT'][loan['id']]
        if payment['type']!='number':raise ValueError('finite payment required')
        for period in range(1,loan['periods']+1):
            for timing in (0.,1.):
                rows.append(dict(loan=loan['id'],period=period,timing=timing,
                    args=[service_client.number(loan['rate']),service_client.number(float(loan['periods']-period+1)),
                          service_client.number(payment['bits']),service_client.number(0.),service_client.number(timing)]))
    write(output/'questions.json',dict(rows=rows,source_sha256=sha(__file__),
        hypothesis='Each remaining-period outstanding balance comes from public PV of the observed PMT; test summed balances, summed interest corrections and per-period principal separately.'))
    pv=service.ask_rows('PV',[r['args']for r in rows])
    part_rows=[dict(loan=l['id'],period=p,args=[service_client.number(v)for v in
        (l['rate'],float(p),float(l['periods']),l['pv'],0.,float(l['timing']))])for l in loans for p in range(1,l['periods']+1)]
    ipmt=service.ask_rows('IPMT',[r['args']for r in part_rows])
    if any(a['type']!='number'for a in pv+ipmt):raise ValueError('finite balance controls required')
    fv_rows=[];fv=[]
    if include_fv:
        for loan in loans:
            for period in range(1,loan['periods']+1):
                for timing in (0.,1.):
                    fv_rows.append(dict(loan=loan['id'],period=period,timing=timing,
                        args=[service_client.number(loan['rate']),service_client.number(float(period-1)),
                              service_client.number(answers['PMT'][loan['id']]['bits']),service_client.number(loan['pv']),service_client.number(timing)]))
        fv=service.ask_rows('FV',[r['args']for r in fv_rows])
        if len(fv)!=len(fv_rows)or any(a['type']!='number'for a in fv):raise ValueError('aligned finite FV controls required')
    write(output/'leaves.json',dict(pv=list(zip(rows,pv)),ipmt=list(zip(part_rows,ipmt)),fv=list(zip(fv_rows,fv)),receipts=service.receipts))
    balances={(r['loan'],r['period'],r['timing']):binary64(a['bits'])for r,a in zip(rows,pv)}
    interests={(r['loan'],r['period']):-binary64(a['bits'])for r,a in zip(part_rows,ipmt)}
    forward_balances={(r['loan'],r['period'],r['timing']):-binary64(a['bits'])for r,a in zip(fv_rows,fv)}
    def op(a,b,k,m):return arithmetic(Q(a),Q(b),k,m)
    def reduce(terms,model,order):
        if order=='balanced'and len(terms)>1:
            k=len(terms)//2;return op(reduce(terms[:k],model,order),reduce(terms[k:],model,order),'add',model)
        value=Q(0)
        for t in terms[::-1]if order=='reverse'else terms:value=op(value,t,'add',model)
        return value
    programs=[];survivors=[]
    for leaf,family,product,addition,summation,order in itertools.product(
        (('PV0','PV1','IPMT','FV0','FV1')if include_fv else('PV0','PV1','IPMT')),('sum_balances','sum_interest','sum_parts'),MODELS,MODELS,MODELS,('forward','reverse','balanced')):
        first=None;checked=0
        for index,(r,wanted)in enumerate(zip(q['ranges'],answers['CUMPRINC'])):
            loan=loans[r['loan']];rate=Q(loan['rate']);payment=binary64(answers['PMT'][r['loan']]['bits'])
            factor=rate if leaf=='PV0'else op(1,op(1,op(1,rate,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64')if leaf=='PV1'else Q(1)
            if leaf.startswith('FV'):
                factor=op(rate,op(1,rate,'add','stored_pc64'),'div','stored_pc64')if loan['timing']else rate
            terms=[]
            for period in range(r['start'],r['end']+1):
                no_interest=loan['timing']and period==1
                term=Q(0)if no_interest else interests[r['loan'],period]if leaf=='IPMT'else forward_balances[r['loan'],period,float(leaf=='FV1')]if leaf.startswith('FV')else balances[r['loan'],period,float(leaf=='PV1')]
                if family!='sum_balances':term=op(term,factor,'mul',product)
                if family=='sum_parts':term=op(payment,term,'add',addition)
                terms.append(term)
            total=reduce(terms,summation,order)
            if family=='sum_balances':total=op(total,factor,'mul',product)
            if family!='sum_parts':total=op(op(r['end']-r['start']+1,payment,'mul','stored_pc64'),total,'add',addition)
            got=bits(float(B64.rn(total)));checked+=1
            if got!=wanted['bits']:
                first=dict(index=index,original=r,expected=wanted,prediction=got);break
        identity=len(programs)
        if first is None:survivors.append(identity)
        programs.append(dict(id=identity,leaf=leaf,family=family,product=product,addition=addition,
            summation=summation,order=order,checked=checked,first_contradiction=first))
    result=dict(format='dnacalc-cumulative-backward-balance-sums/1',complete=True,programs=programs,
        survivors=survivors,receipts=service.receipts,input_sha256={str(p):sha(p)for p in(source/'questions.json',source/'answers.json',output/'leaves.json',Path(__file__))},
        limits='Public PV/FV of the observed PMT and public IPMT are explicit conditional leaves, recomputed independently at each remaining period. This differs from the previously refuted affine recurrences and principal-component sums. Every program refutes at its first exact contradiction; no private PV, payment or sum is identified by these tests.')
    write(output/'result.json',result);print(json.dumps(dict(programs=len(programs),survivors=survivors)))


def summation_analyze(source, destination):
    """Necessary summation constraints, with each observed singleton retained.

    Latent operands range over their complete binary64 publication cells. The
    closed interval enclosure deliberately admits both tie endpoints; failure
    therefore refutes even that enlarged fixed-term premise. Containment is
    necessary only, and never identifies an accumulator or private operand.
    """
    from bessel import operation
    source, destination = Path(source), Path(destination)
    destination.mkdir(exist_ok=False)
    questions = json.loads((source/'questions.json').read_text())
    answers = json.loads((source/'answers.json').read_text())
    requests = questions['ranges']
    def value(answer):
        if answer['type'] != 'number': raise ValueError('non-numeric sum observation')
        return Q(struct.unpack('>d', bytes.fromhex(answer['bits']))[0])
    principal, interest, payment = [[value(v) for v in answers[f]]
                                  for f in ('CUMPRINC', 'CUMIPMT', 'PMT')]
    singletons = {(r['loan'], r['start']): p for r, p in zip(requests, principal)
                  if r['start'] == r['end']}
    def publication(q): return '0x'+struct.pack('>d', float(q)).hex()
    def rounded(q, model):
        return operation('add', model, Q(0), q)
    def reduce(values, model, order):
        if order == 'balanced' and len(values)>1:
            k = len(values)//2
            return operation('add', model, reduce(values[:k], model, order),
                             reduce(values[k:], model, order))
        result = Q(0)
        for v in values[::-1] if order == 'reverse' else values:
            result = operation('add', model, result, v)
        return result
    profiles = []
    for model, order in itertools.product(MODELS, ('forward', 'reverse', 'balanced')):
        matched, contained, first, first_outside = 0, 0, None, None
        for index, (r, expected) in enumerate(zip(requests, principal)):
            terms = [singletons[r['loan'], i] for i in range(r['start'], r['end']+1)]
            got = B64.rn(reduce(terms, model, order))
            matched += got == expected
            if got != expected and first is None:
                first = dict(index=index, request=r, expected=publication(expected),
                             prediction=publication(got), singleton_terms=list(map(publication, terms)))
            cells = [B64.preimage(Interval(t, t)) for t in terms]
            if any(c is None for c in cells): raise ValueError('missing publication cell')
            low = B64.rn(reduce([c.lo for c in cells], model, order))
            high = B64.rn(reduce([c.hi for c in cells], model, order))
            inside = low <= expected <= high
            contained += inside
            if not inside and first_outside is None:
                first_outside = dict(index=index, request=r, expected=publication(expected),
                    low=publication(low), high=publication(high), singleton_terms=list(map(publication, terms)))
        profiles.append(dict(model=model, order=order, rows=len(requests), matched=matched,
            first_published_term_contradiction=first, latent_operand_enclosure_contains=contained,
            first_latent_operand_contradiction=first_outside))
    finish = []
    for model in MODELS:
        matched, first = 0, None
        for index, (r, p, expected) in enumerate(zip(requests, principal, interest)):
            count = Q(r['end']-r['start']+1)
            total = operation('mul', model, count, payment[r['loan']])
            got = B64.rn(operation('add', model, total, -p))
            matched += got == expected
            if got != expected and first is None:
                first = dict(index=index, request=r, expected=publication(expected), prediction=publication(got))
        finish.append(dict(model=model, rows=len(requests), matched=matched, first_contradiction=first))
    result = dict(complete=True, loans=len(questions['loans']), ranges=len(requests),
        sum_programs=profiles, paired_interest_finish=finish,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        inputs_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (source/'questions.json', source/'answers.json')},
        limits='Declared common per-period operand and deterministic reduction programs only. '
               'Closed latent publication cells enlarge the feasible domain. Containment is not proof; '
               'exclusions do not identify term production, loop precision or the complete cumulative kernel.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(sum_programs=[(p['model'], p['order'], p['matched'],
        p['latent_operand_enclosure_contains']) for p in profiles], paired_interest_finish=finish)))


def paired_sum_inverse(source, destination):
    """Read the principal operand through the paired interest subtraction.

    The public PMT value is an explicit premise. An interval intersection is a
    necessary constraint, never a supplied private result or an implementation.
    """
    from bessel import operation
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'checker'))
    from exact import EXTENDED64, BINARY64, Span
    source, destination = Path(source), Path(destination)
    destination.mkdir(exist_ok=False)
    questions = json.loads((source/'questions.json').read_text())
    answers = json.loads((source/'answers.json').read_text())
    def value(answer):
        if answer['type'] != 'number': raise ValueError('non-numeric paired observation')
        return Q(struct.unpack('>d', bytes.fromhex(answer['bits']))[0])
    principal, interest, payment = [[value(v) for v in answers[f]]
                                  for f in ('CUMPRINC', 'CUMIPMT', 'PMT')]
    profiles, recoveries = [], []
    for (product_id, product_model), (final_id, final_model) in itertools.product(enumerate(MODELS), enumerate(MODELS)):
        counts, first = Counter(), None
        payment_domains = {l['id']:B64.preimage(Interval(payment[l['id']], payment[l['id']]))
                           for l in questions['loans']}
        for index, (r, p, i) in enumerate(zip(questions['ranges'], principal, interest)):
            count = Q(r['end']-r['start']+1)
            total = operation('mul', product_model, count, payment[r['loan']])
            target = Interval(i, i)
            if final_model.startswith('resident'): target = B64.preimage(target)
            minus_principal = inverse_add(target, final_id, total)
            candidate = minus_principal.affine(Q(-1)) if minus_principal is not None else None
            published = B64.preimage(Interval(p, p))
            residual_cell = rounding_inverse(target, final_id)
            if residual_cell is None:
                payment_domains[r['loan']] = None
            elif payment_domains[r['loan']] is not None:
                total_cell = Interval(residual_cell.lo+published.lo, residual_cell.hi+published.hi,
                    residual_cell.lc and published.lc, residual_cell.hc and published.hc)
                possible_payment = inverse_op(total_cell, product_id, 'mul', count, True)
                payment_domains[r['loan']] = payment_domains[r['loan']].intersect(possible_payment)
            domain = published.intersect(candidate) if candidate is not None else None
            span = Span(domain.lo, domain.hi, domain.lc, domain.hc) if domain is not None else None
            extended = EXTENDED64.coordinates(span)
            binary = BINARY64.coordinates(span)
            counts['rows'] += 1
            counts['resident_principal_possible'] += extended is not None
            counts['stored_principal_possible'] += binary is not None
            counts['unique_resident_principal'] += extended is not None and extended[0] == extended[1]
            if extended is None and first is None:
                first = dict(index=index, request=r, principal=answers['CUMPRINC'][index],
                    interest=answers['CUMIPMT'][index], payment=answers['PMT'][r['loan']])
            if product_model == 'stored_pc64' and final_model == 'stored_pc64':
                recoveries.append(dict(index=index, loan=r['loan'], start=r['start'], end=r['end'],
                    principal=answers['CUMPRINC'][index], interest=answers['CUMIPMT'][index],
                    payment=answers['PMT'][r['loan']], extended_coordinates=extended,
                    binary64_coordinates=binary,
                    domain=None if domain is None else domain.json()))
        private_payments = []
        for loan, domain in payment_domains.items():
            span = Span(domain.lo, domain.hi, domain.lc, domain.hc) if domain is not None else None
            coordinates = EXTENDED64.coordinates(span)
            private_payments.append(dict(loan=loan, extended_coordinates=coordinates,
                binary64_coordinates=BINARY64.coordinates(span),
                domain=None if domain is None else domain.json()))
        profiles.append(dict(product_model=product_model, final_model=final_model,
            counts=counts, first_contradiction=first,
            common_private_payment_necessary_domains=private_payments,
            loans_with_common_private_payment=sum(p['extended_coordinates'] is not None for p in private_payments)))
    (destination/'recoveries.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in recoveries))
    result = dict(complete=True, loans=len(questions['loans']), profiles=profiles,
        inputs_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in (source/'questions.json', source/'answers.json')},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limits='Conditional on public PMT as the payment operand and the declared multiply/subtract words. '
               'All open/tie endpoints are retained. A nonempty private-principal cell does not identify its producer, '
               'the summation loop or the complete cumulative function, and no recovered operand is installed.')
    (destination/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps([(p['product_model'], p['final_model'], p['counts']) for p in profiles]))


def sum_factor_constraints(source, powers, discount_parts, destination):
    """Constrain one common factor around explicit sums of independently read powers.

    Closed geometric differences are deliberately absent. Each complete program
    uses the same sum word/order and final multiplication word for every range.
    Private factors are constrained only, never installed or selected by scores.
    """
    from bessel import operation
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'checker'))
    from exact import EXTENDED64, BINARY64, Span
    source, powers, discount_parts, destination = map(Path, (source, powers, discount_parts, destination))
    destination.mkdir(exist_ok=False)
    questions=json.loads((source/'questions.json').read_text())
    observed=json.loads((source/'answers.json').read_text())['CUMPRINC']
    public=json.loads(powers.read_text())
    core=json.loads(discount_parts.read_text())
    pub={(r['rate'],r['exponent']):number('0x'+a['bits'])
         for r,a in zip(public['rows'],public['reply']['rows'])}
    core={r['loan']['id']:r for r in core}
    def reduce(terms,model,order):
        if order=='balanced'and len(terms)>1:
            half=len(terms)//2
            return operation('add',model,reduce(terms[:half],model,order),reduce(terms[half:],model,order))
        result=Q(0)
        for term in terms[::-1]if order=='reverse'else terms:result=operation('add',model,result,term)
        return result
    programs=[];survivors=[]
    for family, sum_model, order, (final_id,final_model) in itertools.product(
            ('growth_public','discount_public','reciprocal_public','discount_core'), MODELS,
            ('forward','reverse','balanced'), enumerate(MODELS)):
        states={};first=None;excluded=Counter();scoped=0
        for index,(r,answer)in enumerate(zip(questions['ranges'],observed)):
            loan=questions['loans'][r['loan']]
            if loan['timing']==1 and r['start']==1:
                excluded['separate_first_beginning_payment']+=1;continue
            if family=='discount_core'and core[r['loan']].get('outside_observation_scope'):
                excluded['outside_existing_discount_core_branch']+=1;continue
            degrees=[i-1-loan['timing']for i in range(r['start'],r['end']+1)]
            if family=='growth_public':terms=[Q(pub[loan['rate'],i])for i in degrees]
            elif family=='discount_public':terms=[Q(pub[loan['rate'],i-loan['periods']])for i in degrees]
            elif family=='reciprocal_public':
                terms=[B64.rn(1/Q(pub[loan['rate'],loan['periods']-i]))for i in degrees]
            else:terms=[Q(number(core[r['loan']]['discounts'][i]))for i in degrees]
            total=reduce(terms,sum_model,order)
            if total<=0:raise ValueError('positive sum premise failed')
            p=Q(number('0x'+answer['bits']));target=Interval(p,p)
            if final_model.startswith('resident'):target=B64.preimage(target)
            required=inverse_op(target,final_id,'mul',total,True)
            before=states.get(r['loan'])
            if r['loan']not in states:states[r['loan']]=required
            elif before is not None:states[r['loan']]=before.intersect(required)
            scoped+=1
            if states[r['loan']]is None and first is None:
                first=dict(index=index,request=r,expected=answer,
                    prior=None if before is None else before.json(),
                    required=None if required is None else required.json())
        feasible=[]
        for loan,domain in states.items():
            span=Span(domain.lo,domain.hi,domain.lc,domain.hc)if domain is not None else None
            feasible.append(dict(loan=loan,extended=EXTENDED64.coordinates(span),binary64=BINARY64.coordinates(span),
                domain=None if domain is None else domain.json()))
        good=sum(r['extended']is not None for r in feasible)
        identity=len(programs)
        if good==len(states):survivors.append(identity)
        programs.append(dict(id=identity,family=family,sum_model=sum_model,order=order,final_model=final_model,
            rows_scoped=scoped,excluded=excluded,loans=len(states),loans_with_common_factor=good,
            first_contradiction=first,necessary_factor_domains=feasible if good==len(states)else None))
    result=dict(complete=True,original_rows=len(questions['ranges']),programs=programs,survivors=survivors,
        inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (source/'questions.json',source/'answers.json',powers,discount_parts)},
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limits='One arbitrary common PC64 factor per loan, explicit observed POWER/EXP terms and fixed summation programs. '
               'Necessary factor cells neither identify a numerical factor producer nor prove a complete cumulative implementation. '
               'First beginning-period payments and unsupported existing discount-core branches remain explicitly outside each scope, never passed.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(programs=len(programs),survivors=survivors)))


def stored_component_sums(source, components, destination):
    """Fixed sums of stored principal components, with timing applied after it.

    This differs from a common factor around a sum of raw powers: h*u is
    stored for each period before accumulation. Existing product h/u operands
    are explicit conditional leaves, never recovered from cumulative answers.
    """
    from bessel import operation,pair
    source,components,destination=map(Path,(source,components,destination));destination.mkdir(exist_ok=False)
    q=json.loads((source/'questions.json').read_text());a=json.loads((source/'answers.json').read_text())['CUMPRINC']
    core={r['loan']['id']:r for r in json.loads(components.read_text())}
    programs=[];survivors=[]
    def op(x,y,k,m):return operation(k,m,Q(x),Q(y))
    def reduce(terms,m,order):
        if order=='balanced'and len(terms)>1:
            half=len(terms)//2;return op(reduce(terms[:half],m,order),reduce(terms[half:],m,order),'add',m)
        value=Q(0)
        for term in terms[::-1]if order=='reverse'else terms:value=op(value,term,'add',m)
        return value
    for family,sum_model,final_model,order in itertools.product(
            ('components_then_timing','shifted_components_then_rate','published_parts'),
            MODELS,MODELS,('forward','reverse','balanced')):
        first=None;checked=0;excluded=Counter()
        for index,(row,expected)in enumerate(zip(q['ranges'],a)):
            loan=q['loans'][row['loan']];parts=core[row['loan']]
            if parts.get('outside_observation_scope'):
                excluded['outside_current_component_producer']+=1;continue
            h=Q(number(parts['quotient']));r=Q(loan['rate']);begin=bool(loan['timing'])
            due=op(1,op(1,op(1,r,'div','stored_pc64'),'add','stored_pc64'),'div','stored_pc64')
            factor=due if begin else r;initial=Q(0);terms=[]
            for period in range(row['start'],row['end']+1):
                if begin and period==1:
                    if family=='shifted_components_then_rate':initial=op(h,due,'mul','stored_pc64');continue
                    term=h
                else:
                    offset=period-2 if family=='shifted_components_then_rate'and begin else period-1
                    term=op(h,Q(number(parts['discounts'][offset])),'mul','stored_pc64')
                if family=='published_parts':term=op(term,factor,'mul',final_model)
                terms.append(term)
            total=reduce(terms,sum_model,order)
            if family=='shifted_components_then_rate':got=op(initial,op(total,r,'mul',final_model),'add',final_model)
            elif family=='published_parts':got=total
            else:got=op(total,factor,'mul',final_model)
            prediction='0x'+struct.pack('>d',float(B64.rn(got))).hex();wanted='0x'+expected['bits'];checked+=1
            if prediction!=wanted:
                first=dict(index=index,original=row,expected=expected,prediction=prediction,
                    loan=loan,terms=[pair(t)for t in terms],initial=pair(initial));break
        identity=len(programs)
        if first is None:survivors.append(identity)
        programs.append(dict(id=identity,family=family,sum_model=sum_model,final_model=final_model,
            order=order,checked_until_first_contradiction=checked,excluded=excluded,first_contradiction=first))
    result=dict(complete=True,original_rows=len(a),programs=programs,survivors=survivors,
        source_sha256={str(source/'questions.json'):hashlib.sha256((source/'questions.json').read_bytes()).hexdigest(),
            str(source/'answers.json'):hashlib.sha256((source/'answers.json').read_bytes()).hexdigest(),
            str(components):hashlib.sha256(components.read_bytes()).hexdigest(),str(Path(__file__).resolve()):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        limits='Finite stored-component sum hypotheses only. Existing h/u producers are conditional product leaves; no observed principal or selected factor is installed. Unsupported rate1 rows are outside every scope, never scored as passes.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(programs=len(programs),survivors=survivors)))


def recurrence_sum_domains(source,destination,balance=False,backward=False,interior=False):
    """Necessary initial-value cells for range-local geometric accumulation.

    Forward ranges share an arbitrary extended first term per (loan,start);
    backward ranges share an arbitrary last term per (loan,end). Each term is
    constrained by its exact singleton. No different range origin is required
    to use the same rounded producer.
    """
    from bessel import operation
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checker'))
    from exact import EXTENDED64,Span
    if balance and backward:raise ValueError('backward balance family is not specified')
    source,destination=map(Path,(source,destination));destination.mkdir(exist_ok=False)
    questions=json.loads((source/'questions.json').read_text());observations=json.loads((source/'answers.json').read_text());answers=observations['CUMPRINC']
    singleton={(r['loan'],r['start']):Q(number('0x'+v['bits']))for r,v in zip(questions['ranges'],answers)if r['start']==r['end']}
    profiles=[(m,m)for m in MODELS]+[('stored_pc64','strict')]
    programs=[];survivors=[]
    def op(a,b,k,m):return operation(k,m,Q(a),Q(b))
    descriptions=list(itertools.product(('principal','component'),('multiply_growth','add_interest'),profiles,MODELS))
    if backward:
        descriptions=list(itertools.product(('principal','component'),('divide_growth','multiply_reciprocal'),profiles,MODELS))
    if balance:
        descriptions=[(family,advance,profile,finish)for family in('principal_from_balance','interest_balances')
            for advance in('add_principal','multiply_growth_add_payment')for profile in profiles
            for finish in(('strict',)if family=='principal_from_balance'else MODELS)]
    for family,advance,(update_word,sum_word),finish_word in descriptions:
        domains={};first=None;checked=0;excluded=Counter();searches=0
        indexed=list(enumerate(zip(questions['ranges'],answers)))
        if backward:indexed.sort(key=lambda item:(item[1][0]['loan'],item[1][0]['end'],-item[1][0]['start']))
        for index,(row,answer)in indexed:
            loan=questions['loans'][row['loan']];start,end=row['start'],row['end']
            if interior and start==1:excluded['first_period_branch_outside_scope']+=1;continue
            if loan['timing']and start==1:excluded['separate_first_beginning_payment']+=1;continue
            r=Q(loan['rate']);growth=op(1,r,'add','stored_pc64');factor=r if family=='component'else Q(1);payment=Q(number('0x'+observations['PMT'][row['loan']]['bits']))
            key=(row['loan'],end if backward else start)
            if key not in domains:
                target=Interval(singleton[key],singleton[key])
                if balance and family=='principal_from_balance':
                    if sum_word.startswith('resident'):target=B64.preimage(target)
                    initial=inverse_add(target,MODELS.index(sum_word),Q(0))
                    initial=inverse_add(initial,MODELS.index(update_word),payment)if initial else None
                    initial=inverse_op(initial,MODELS.index(update_word),'mul',r,True)if initial else None
                elif balance:
                    if finish_word.startswith('resident'):target=B64.preimage(target)
                    first_payment=op(payment,1,'mul',finish_word)
                    initial=inverse_add(target,MODELS.index(finish_word),first_payment)
                    initial=inverse_op(initial,MODELS.index(finish_word),'mul',r,True)if initial else None
                    initial=inverse_add(initial,MODELS.index(sum_word),Q(0))if initial else None
                else:
                    if finish_word.startswith('resident'):target=B64.preimage(target)
                    initial=inverse_op(target,MODELS.index(finish_word),'mul',factor,True)
                    initial=inverse_add(initial,MODELS.index(sum_word),Q(0))if initial else None
                span=Span(initial.lo,initial.hi,initial.lc,initial.hc)if initial else None
                domains[key]=EXTENDED64.coordinates(span)
            prior=domains[key]
            expected=Q(number('0x'+answer['bits']))
            def evaluate(coordinate):
                term=EXTENDED64.value(coordinate)
                if balance:
                    total=Q(0)
                    for step in range(end-start+1):
                        principal=op(payment,op(term,r,'mul',update_word),'add',update_word)
                        total=op(total,principal if family=='principal_from_balance'else term,'add',sum_word)
                        if step<end-start:
                            term=op(term,principal,'add',update_word)if advance=='add_principal'else op(op(term,growth,'mul',update_word),payment,'add',update_word)
                    if family=='principal_from_balance':return B64.rn(total)
                    return B64.rn(op(op(payment,end-start+1,'mul',finish_word),op(total,r,'mul',finish_word),'add',finish_word))
                total=op(0,term,'add',sum_word)
                for _ in range(end-start):
                    if backward:
                        term=op(term,growth,'div',update_word)if advance=='divide_growth'else op(term,op(1,growth,'div','stored_pc64'),'mul',update_word)
                    else:
                        term=op(term,growth,'mul',update_word)if advance=='multiply_growth'else op(term,op(term,r,'mul',update_word),'add',update_word)
                    total=op(total,term,'add',sum_word)
                return B64.rn(op(total,factor,'mul',finish_word))
            if prior is None:updated=None
            else:
                low,high=prior
                def bound(strict):
                    nonlocal searches
                    lo,hi=low,high+1
                    while lo<hi:
                        mid=(lo+hi)//2;value=evaluate(mid);searches+=1
                        if value<expected or strict and value==expected:lo=mid+1
                        else:hi=mid
                    return lo
                left,right=bound(False),bound(True)-1;updated=(left,right)if left<=right else None
            domains[key]=updated;checked+=1
            if updated is None:
                first=dict(index=index,original=row,expected=answer,loan=loan,previous_initial_coordinates=prior,
                    previous_endpoint_predictions=None if prior is None else ['0x'+struct.pack('>d',float(evaluate(k))).hex()for k in prior]);break
        identity=len(programs)
        if first is None:survivors.append(identity)
        programs.append(dict(id=identity,family=family,advance=advance,update_word=update_word,sum_word=sum_word,
            finish_word=finish_word,checked_until_first_contradiction=checked,excluded=excluded,exact_search_evaluations=searches,
            first_contradiction=first,necessary_initial_domains=[dict(loan=k[0],**({'end':k[1]}if backward else{'start':k[1]}),coordinates=v)for k,v in domains.items()]if first is None else None))
    result=dict(complete=True,original_rows=len(answers),programs=programs,survivors=survivors,
        source_sha256={str(source/'questions.json'):hashlib.sha256((source/'questions.json').read_bytes()).hexdigest(),str(source/'answers.json'):hashlib.sha256((source/'answers.json').read_bytes()).hexdigest(),str(Path(__file__).resolve()):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        public_payment_leaf=balance,limits=('Interior-only necessary loop constraints: every range beginning at period1 is excluded explicitly, without asserting that Excel uses a separate branch. No excluded range counts as a pass. 'if interior else'')+('Backward principal/component recurrences: one arbitrary Ext80 last term per(loan,end), constrained by its exact singleton, then divided by stored1+r or multiplied by its stored reciprocal while accumulating down to each start. No last-term producer is assumed. 'if backward else'')+'Arbitrary extended initial values are necessary constraints only, never selected operands. First beginning-payment ranges are outside each recurrence scope. Exact monotone bisection intersects complete publication cells across the range; no first-term producer or complete CUMPRINC kernel is identified by a nonempty domain.')
    (destination/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(programs=len(programs),survivors=survivors)))


def tagged(s):
    if s.startswith('binary64:'):return Q(number(s.split(':')[1]))
    assert s.startswith('ext80-le:')
    raw=bytes.fromhex(s.split(':')[1]);m=int.from_bytes(raw[:8],'little');e=int.from_bytes(raw[8:],'little')
    assert e&0x7fff!=0x7fff
    return (-1 if e>>15 else 1)*m*pow2(((e&0x7fff)or 1)-16383-63)


def group(row):
    return json.dumps([row['args'][0],row['args'][1],row['source'],row['reading']],sort_keys=True,separators=(',',':'))


def run(trace,growth,out,companions=None):
    started=time.monotonic();out=Path(out);out.mkdir(parents=True,exist_ok=False)
    rows=[json.loads(line)for line in open(trace)]  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
    groups=defaultdict(list)
    for i,row in enumerate(rows):
        a=row['args']
        if (companions or number(a[5])==0) and number(a[3])==number(a[4]):
            groups[group(row)+':'+a[2]+':'+a[5]].append(i)
    payment_map={}
    if companions:
        for row in map(json.loads,open(companions)):  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
            a=row['args'];assert row['function']=='PMT'and len(a)==5 and number(a[3])==0
            key=tuple(a[:3]+[a[4],json.dumps(row['reading'],sort_keys=True)])
            assert key not in payment_map;payment_map[key]=row
    helpers={}
    for row in map(json.loads,open(growth)):  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
        # Canonicalize the JSON group key without changing its typed contents.
        key=json.dumps(json.loads(row['group']),sort_keys=True,separators=(',',':'))
        helpers[key,row['power'],row['growth_model']]=row
    profiles=[];records=[];survivors=[];kill_counts=Counter()
    for pid,(power,gm,consumer)in enumerate(itertools.product(range(2),MODELS,MODELS)):
        profile=dict(id=pid,power=power,growth_model=gm,consumer_model=consumer,
            family='shared_payment'if companions else'shared_initial_principal')
        profiles.append(profile);killed=False
        for loan,indices in sorted(groups.items()):
            current=None;preimages=[];payment=None
            if companions:
                raw=rows[indices[0]];a=raw['args']
                payment=payment_map[tuple(a[:3]+[a[5],json.dumps(raw['reading'],sort_keys=True)])]
                y=Q(number(payment['expected']));current=B64.preimage(Interval(y,y))
            for i in indices:
                raw=rows[i];period=int(number(raw['args'][3]));h=helpers[group(raw),power,gm]
                degree=int(number(raw['args'][1]))-period+1 if companions else period-1
                g=tagged(h['factors'][degree]);y=Q(number(raw['expected']))
                published=B64.preimage(Interval(y,y))
                required=(published if companions and number(raw['args'][5])==1 and period==1
                    else inverse_op(published,MODELS.index(consumer),'div'if companions else'mul',g,True))
                current=(required if not preimages and not companions
                    else current.intersect(required)if current else None)
                preimages.append(dict(row_index=i,period=period,growth=growth,
                    required=required.json()if required else None,
                    intersection=current.json()if current else None))
                if current is None or E64.preimage(current)is None:break
            feasible=current is not None and E64.preimage(current)is not None
            record=dict(profile_id=pid,loan=loan,preimages=preimages,
                payment_observation=payment,feasible=feasible,
                reason=('nonempty'if feasible else 'empty_real_intersection'if current is None else 'empty_ext80_lattice'))
            records.append(record)
            if not feasible:
                killed=True;kill_counts[record['reason']]+=1;break
        if not killed:survivors.append(pid)
    (out/'profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
    with (out/'constraints.jsonl').open('x')as f:
        for record in records:f.write(json.dumps(record,separators=(',',':'))+'\n')
    (out/'survivors.json').write_text(json.dumps(survivors)+'\n')
    summary=dict(profiles=50,survivors=len(survivors),original_available_rows=len(rows),
        scoped_singletons=sum(map(len,groups.values())),loans=len(groups),kill_counts=kill_counts,
        seconds=time.monotonic()-started,
        assumptions=('One arbitrary common Ext80 payment per loan, constrained by its independently recorded PMT publication; singleton=RN53(model(P/G(nper-period+1))), with first beginning-timing payment returned directly. No PMT-producing program assumed.'if companions else
            'One arbitrary common Ext80 A per exact loan/environment: singleton output RN53(model(A*G(period-1))). G is the explicitly calculated integer-growth hypothesis. No PMT-producing program or fixed A is assumed.'),
        claim='Necessary shared-latent constraints only; no coefficients fitted, no complete numerical program or Excel helper identified.',
        inputs={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest()for p in ([trace,growth]+([companions]if companions else[]))})
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


def questions(trace,out):
    """Freeze a coupled surface/transport audit after zero complete survivors."""
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    original=[json.loads(line)for line in open(trace)]  # ubs:ignore — finite one-shot iterator; fully consumed and released before process exit
    bits=lambda x:'0x'+struct.pack('>d',x).hex()
    first=original[0]['args']
    loans=[first[:3],[bits(0.05/12),bits(360.0),bits(200000.0)]]
    batches={f:[]for f in ('CUMPRINC','PMT','PPMT','CUMIPMT')}
    selection=[];relations=[]
    def add(function,args,role):
        qid=f'r11-{function.lower()}-{len(batches[function]):02}'
        batches[function].append(dict(id=qid,args=args,input_transport='cell_value2_bulk'))
        known=[dict(trace_index=i,witness=r)for i,r in enumerate(original)
               if r['function']==function and r['args']==args]
        selection.append(dict(id=qid,function=function,role=role,known_bank_readings=known,
            validation_exposure='all old banks and corpus are discovery; no frozen complete candidate exists'))
        return qid
    for loan_index,(rate,nper,pv)in enumerate(loans):
        n=int(number(nper));middle=n//2
        for timing in (0,1):
            t=bits(float(timing));payment=add('PMT',[rate,nper,pv,bits(0.0),t],f'loan{loan_index} timing{timing}: shared-payment publication')
            cumulative={}
            for start,end in ((1,1),(2,2),(middle,middle),(n,n),(1,middle),(middle+1,n),(1,n),(1,2)):
                cumulative[start,end]=add('CUMPRINC',[rate,nper,pv,bits(float(start)),bits(float(end)),t],
                    f'loan{loan_index} timing{timing}: singleton, prefix, tail or whole interval')
            for period in (1,2,middle,n):
                principal=add('PPMT',[rate,bits(float(period)),nper,pv,bits(0.0),t],
                    f'loan{loan_index} timing{timing}: corresponding independently published principal')
                relations.append(dict(kind='singleton CUMPRINC versus standalone PPMT',
                    CUMPRINC=cumulative[period,period],PPMT=principal,
                    question='Do their publications coincide? Different results reject a simple sum of standalone published PPMT; equal results do not identify internal sharing.'))
            if timing==1:
                relations.append(dict(kind='first beginning-timing principal/payment control',
                    CUMPRINC=cumulative[1,1],PMT=payment,
                    question='Do the publications coincide, as in all30old paired loans?'))
            for start,end in ((1,middle),(middle+1,n),(1,n)):
                interest=add('CUMIPMT',[rate,nper,pv,bits(float(start)),bits(float(end)),t],
                    f'loan{loan_index} timing{timing}: payments less principal relation')
                relations.append(dict(kind='conditional arithmetic relation, not a whole-function prediction',
                    CUMPRINC=cumulative[start,end],CUMIPMT=interest,PMT=payment,count=end-start+1,
                    rule='For each of25fixed product/subtraction model pairs, intersect publication preimages of one shared Ext80 payment and principal with interest=RN53(subtract_model(product_model(count*payment),principal)). No coefficient is selected or fitted.',
                    models=MODELS))
    assert {k:len(v)for k,v in batches.items()}==dict(CUMPRINC=32,PMT=4,PPMT=16,CUMIPMT=12)
    for function,rows in batches.items():
        assert len({tuple(r['args'])for r in rows})==len(rows)
        (out/f'requests-{function}.json').write_text(json.dumps(dict(function=function,witnesses=rows),indent=2)+'\n')
    (out/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    (out/'conditional-relations.json').write_text(json.dumps(relations,indent=2)+'\n')
    readme='''Round 11: 64 coupled annuity observations. No new Excel answers have been read.
Requests contain IDs/arguments/transport only. Application predictions and
conditional diagnostic relations are separate research files, never answers.
The 90,000 complete programs and both 50-profile singleton projections have no
survivors. These questions audit helper/surface/transport assumptions; they
are not discriminators for surviving complete programs or fresh validation of one.
Use all four batches below on Windows; do not use --new-only (old controls repeat).
Score the immutable application predictions first when answers arrive, then
evaluate the frozen coupled relations. Keep requested/submitted/read-back inputs,
argument types, settings, raw results and failed obligations. No predictions are
to be installed into an evidence bank.

cargo xtask witness read --func CUMPRINC --inputs requests-CUMPRINC.json
cargo xtask witness read --func PMT --inputs requests-PMT.json
cargo xtask witness read --func PPMT --inputs requests-PPMT.json
cargo xtask witness read --func CUMIPMT --inputs requests-CUMIPMT.json
'''
    (out/'RUN-ROUND11.txt').write_text(readme)
    print(json.dumps({k:len(v)for k,v in batches.items()}))


def coupled(packet, trace, out):
    """Apply the frozen R11 payment/principal/interest relation, without a producer."""
    packet, trace, out = map(Path, (packet, trace, out))
    out.mkdir(parents=True, exist_ok=False)
    rows = {r['id']:r for r in map(json.loads, trace.read_text().splitlines())}
    relations = json.loads((packet/'conditional-relations.json').read_text())
    r15 = json.loads((packet/'manifest.json').read_text()).get('calls') == 236
    loans = defaultdict(list); identities = []
    for relation in relations:
        if relation.get('kind') == 'published payment partition audit':
            continue  # The component publication cells are checked separately.
        c = rows[relation['CUMPRINC']]
        assert c['function'] == 'CUMPRINC'
        if 'PPMT' in relation:
            p = rows[relation['PPMT']]
            a, b = c['args'], p['args']
            assert p['function'] == 'PPMT' and a[3] == a[4] == b[1]
            assert a[:3] == [b[0], b[2], b[3]] and b[4] == '0x0000000000000000' and a[5] == b[5]
            assert c['reading'] == p['reading']
            identities.append(dict(relation=relation, equal=c['expected']==p['expected'],
                                   cumulative=c, standalone=p))
        elif 'CUMIPMT' in relation:
            p, i = rows[relation['PMT']], rows[relation['CUMIPMT']]
            assert p['function']=='PMT' and i['function']=='CUMIPMT'
            assert c['args']==i['args'] and c['args'][:3]==p['args'][:3]
            assert p['args'][3]=='0x0000000000000000' and c['args'][5]==p['args'][4]
            assert c['reading']==i['reading']==p['reading']
            assert relation['count']==int(number(c['args'][4])-number(c['args'][3])+1)
            loans[p['id']].append((relation,p,c,i))
        else:
            p=rows[relation['PMT']]
            assert c['args'][3:6]==['0x3ff0000000000000']*3 and c['args'][:3]==p['args'][:3]
            assert c['reading']==p['reading']
            identities.append(dict(relation=relation,equal=c['expected']==p['expected'],
                                   cumulative=c,standalone=p))
    assert (len(loans), sum(map(len,loans.values())), len(identities)) == ((8,38,20) if r15 else (4,12,18))
    def rounded(x,m):
        if m==0:return B64.rn(x)
        z=(E64 if m in (1,2) else E53).rn(x)
        return B64.rn(z) if m in (1,3) else z
    results=[];survivors=[]
    for pid,(pm,sm) in enumerate(itertools.product(range(5),repeat=2)):
        feasible=True
        for loan, triplets in sorted(loans.items()):
            payment=triplets[0][1];y=Q(number(payment['expected']));cell=B64.preimage(Interval(y,y))
            low,high=E64.ceil(cell.lo),E64.floor(cell.hi)
            if low==cell.lo and not cell.lc:low=E64.next(low)
            if high==cell.hi and not cell.hc:high=E64.prev(high)
            assert low<=high<0 and log2floor(-low)==log2floor(-high)
            step=pow2(log2floor(-low)-63);count=int((high-low)/step)+1
            assert 1<=count<=4096 and low+(count-1)*step==high
            admitted=[]
            for k in range(count):
                latent=low+k*step;possible=True
                for relation,p,c,i in triplets:
                    product=rounded(relation['count']*(B64.rn(latent)if pm==0 else latent),pm)
                    principal=Q(number(c['expected']));interest=Q(number(i['expected']))
                    required=inverse_add(B64.preimage(Interval(interest,interest)),sm,product)
                    required=required.affine(Q(-1))if required else None
                    wanted=B64.preimage(Interval(principal,principal))
                    overlap=wanted.intersect(required) if required else None
                    if overlap is None:possible=False;break
                    first,last=E64.ceil(overlap.lo),E64.floor(overlap.hi)
                    if first==overlap.lo and not overlap.lc:first=E64.next(first)
                    if last==overlap.hi and not overlap.hc:last=E64.prev(last)
                    if first is None or last is None or first>last:possible=False;break
                admitted.append(possible)
            feasible &= any(admitted)
            results.append(dict(profile_id=pid,product_model=MODELS[pm],subtraction_model=MODELS[sm],
                payment=payment,relations=[r[0]for r in triplets],domain_count=count,
                admissible_payment_count=sum(admitted),
                admission_sha256=hashlib.sha256(bytes(admitted)).hexdigest()))
        if feasible:survivors.append(pid)
    (out/'relations.json').write_text(json.dumps(results,indent=2)+'\n')
    (out/'identities.json').write_text(json.dumps(identities,indent=2)+'\n')
    summary=dict(complete=True,profiles=25,survivors=survivors,loans=len(loans),intervals=sum(map(len,loans.values())),
        singleton_pairs=sum('PPMT'in r['relation']for r in identities),singleton_equal=sum(r['equal']for r in identities if 'PPMT'in r['relation']),
        first_beginning_equal=sum(r['equal']for r in identities if 'PPMT'not in r['relation']),
        inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in (trace,packet/'conditional-relations.json',Path(__file__))},
        limits='Frozen necessary relation only: one shared Ext80 payment per loan, and independent '
               'Ext80 principal per interval, constrained by all three publications. No latent '
               'value is selected or installed; nonempty domains do not identify a complete program.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


def balance_domains(trace,structure,out,variant=None):
    """Complete first-singleton preimages; no published-PMT assumption."""
    from intervals import inverse_add,log2floor,E53
    assert variant in (None,'interest');interest_sum=variant=='interest'
    rows=[json.loads(line)for line in Path(trace).read_text().splitlines()]
    groups=json.loads(Path(structure).read_text())['groups'];domains=[]
    def encoding(x):
        assert x<0
        e=log2floor(-x);m=(-x)/pow2(e-63)
        assert m.denominator==1 and 1<<63<=m<1<<64 and -16382<=e<=16383
        return 'ext80-le:'+(int(m).to_bytes(8,'little')+(e+16383+32768).to_bytes(2,'little')).hex()
    def rounded(x,model):
        if model==0:return B64.rn(x)
        y=(E64 if model in (1,2)else E53).rn(x)
        return B64.rn(y)if model in (1,3)else y
    for loan,group in enumerate(groups):
        index=group['row_indices'][0];row=rows[index];a=list(map(lambda x:Q(number(x)),row['args']))
        assert a[3]==a[4]==1
        for store,im,pm,sm,order in itertools.product((False,True),range(5),range(5),range(5),range(3)if interest_sum else(0,)):
            y=Q(number(row['expected']));required=B64.preimage(Interval(y,y))
            # For this finite bank, the +0 and direct pairwise paths have equal
            # applicable Ext80 lattices; their continuous real preimages can
            # differ. The independent checker explicitly verifies lattice equality.
            interest=rounded(a[0]*a[2],im)
            if interest_sum:
                interest=Q(0)if a[5]else interest
                if order!=2:interest=rounded(interest,sm)
                required=inverse_add(required,pm,interest)
                required=inverse_op(required,pm,'mul',Q(1),True)
            else:
                required=inverse_add(required,sm,Q(0))
                if a[5]==0:required=inverse_add(required,pm,interest)
            if store:required=B64.preimage(required)
            lo=E64.ceil(required.lo)if required else None;hi=E64.floor(required.hi)if required else None
            if lo is not None and lo==required.lo and not required.lc:lo=E64.next(lo)
            if hi is not None and hi==required.hi and not required.hc:hi=E64.prev(hi)
            empty=lo is None or hi is None or lo>hi
            if not empty:
                assert lo<=hi<0 and log2floor(-lo)==log2floor(-hi),'domain crosses an Ext80 binade'
                step=pow2(log2floor(-lo)-63);count=(hi-lo)/step+1
                assert count.denominator==1 and 1<=count<=100000
            else:count=0
            domains.append(dict(loan_index=loan,row_index=index,payment_store=store,interest_model=MODELS[im],
                principal_model=MODELS[pm],sum_model=MODELS[sm],required=required.json()if required else None,
                lower=encoding(lo)if not empty else None,upper=encoding(hi)if not empty else None,count=int(count)))
            if interest_sum:domains[-1]['sum_order']=order
    doc=dict(kind='complete first-singleton Ext80 domains; external PMT not a constraint',domains=domains,
        trace_sha256=hashlib.sha256(Path(trace).read_bytes()).hexdigest(),
        structure_sha256=hashlib.sha256(Path(structure).read_bytes()).hexdigest())
    if interest_sum:doc['interest_sum']=True
    with Path(out).open('x')as f:json.dump(doc,f,separators=(',',':'));f.write('\n')
    print(json.dumps(dict(domains=len(domains),minimum=min(d['count']for d in domains),maximum=max(d['count']for d in domains))))


def questions_r15(trace, executable, out):
    """Extend the coupled audit to all nine annuity surfaces using fixed anchors."""
    import subprocess
    trace, executable, out = Path(trace), Path(executable).resolve(), Path(out)
    out.mkdir(parents=True,exist_ok=False)
    observations=list(map(json.loads,trace.read_text().splitlines()))
    byid={r['id']:r for r in observations}
    anchors=('r11-pmt-00','r11-pmt-02','r12-pmt-group-00','r12-pmt-group-05')
    bits=lambda x:'0x'+struct.pack('>d',x).hex()
    zero,one,minusone=bits(0.0),bits(1.0),bits(-1.0)
    functions=('PMT','IPMT','PPMT','CUMPRINC','CUMIPMT','FV','PV','NPER','RATE')
    batches={f:[]for f in functions};seen={f:{}for f in functions};selection=[];relations=[]
    def add(f,args,anchor,timing,role,**metadata):
        key=tuple(args)
        if key not in seen[f]:
            qid=f'r15-{f.lower()}-{len(batches[f]):03d}'
            seen[f][key]=qid
            batches[f].append(dict(id=qid,args=args,input_transport='cell_value2_bulk'))
        qid=seen[f][key]
        selection.append(dict(id=qid,function=f,anchor_id=anchor,timing=timing,role=role,**metadata))
        return qid
    for anchor in anchors:
        source=byid[anchor];r,n,pv,fv,_=source['args'];periods=int(number(n))
        assert source['function']=='PMT'and number(fv)==0 and periods>=1
        selected=sorted(set(k for k in (1,2,periods//2,periods)if 1<=k<=periods))
        intervals=sorted(set([(k,k)for k in (1,2,periods)if k<=periods]+[(1,periods)]+
            ([(1,periods//2),(periods//2+1,periods)]if periods>1 else[])))
        for t in (0,1):
            tb=bits(float(t));matches=[x for x in observations if x['function']=='PMT'and x['args']==[r,n,pv,zero,tb]]
            assert matches and len({x['expected']for x in matches})==1
            payment=matches[0];p=payment['expected'];assert number(p)<0
            ids={}
            for scale in (1,2,-1):
                ids[scale]=add('PMT',[r,n,bits(number(pv)*scale),zero,tb],anchor,t,
                               'cashflow scaling control',scale=scale)
            for k in selected:
                args=[r,bits(float(k)),n,pv,zero,tb]
                i=add('IPMT',args,anchor,t,'standalone interest',period=k)
                q=add('PPMT',args,anchor,t,'standalone principal',period=k)
                relations.append(dict(kind='published payment partition audit',PMT=ids[1],IPMT=i,PPMT=q,
                    limitation='A failed publication identity does not exclude a shared unrounded component.'))
            for a,b in intervals:
                args=[r,n,pv,bits(float(a)),bits(float(b)),tb]
                c=add('CUMPRINC',args,anchor,t,'cumulative principal',start=a,end=b)
                i=add('CUMIPMT',args,anchor,t,'cumulative interest',start=a,end=b)
                relations.append(dict(kind='coupled payment/principal/interest preimages',PMT=ids[1],
                    CUMPRINC=c,CUMIPMT=i,count=b-a+1,
                    limitation='Compare fixed store/precision profiles with one shared latent payment, never fit it.'))
                if a==b:
                    pp=seen['PPMT'][tuple([r,bits(float(a)),n,pv,zero,tb])]
                    ip=seen['IPMT'][tuple([r,bits(float(a)),n,pv,zero,tb])]
                    relations.append(dict(kind='singleton surface audit',CUMPRINC=c,PPMT=pp,CUMIPMT=i,IPMT=ip))
            for f in ('FV','PV'):
                add(f,[r,n,zero,minusone,tb],anchor,t,'unit growth or discount control')
                add(f,[r,n,minusone,zero,tb],anchor,t,'unit annuity control')
                add(f,[r,n,p,pv if f=='FV'else zero,tb],anchor,t,'published-payment cross function control',
                    payment_source_id=payment['id'])
                if periods>1:
                    add(f,[r,bits(float(periods-1)),zero,minusone,tb],anchor,t,'adjacent period growth control')
            add('NPER',[r,p,pv,zero,tb],anchor,t,'inverse with original published payment',payment_source_id=payment['id'])
            for guess in (r,bits(0.1)):
                add('RATE',[n,p,pv,zero,tb,guess],anchor,t,'explicit guess sensitivity control',
                    payment_source_id=payment['id'],guess='original rate'if guess==r else'mathematical 0.1')
    batchlist=[]
    for f in functions:
        request=out/f'requests-{f}.json'
        request.write_text(json.dumps(dict(function=f,witnesses=batches[f]),indent=2)+'\n')
        prediction=out/f'application-predictions-{f}.jsonl'
        subprocess.run([str(executable),'invoke-inputs',str(request),str(prediction)],check=True,
                       timeout=300)  # ubs:ignore — explicit local executable/argv, no shell; manifest binds executable hash
        batchlist.append(dict(function=f,path=request.name,calls=len(batches[f])))
    (out/'batches.json').write_text(json.dumps(batchlist,indent=2)+'\n')
    (out/'selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    (out/'conditional-relations.json').write_text(json.dumps(relations,indent=2)+'\n')
    (out/'RUN-R15.ps1').write_text('''param([Parameter(Mandatory=$true)][string]$DnaCalcCheckout)
$ErrorActionPreference = 'Stop'
$packet = $PSScriptRoot
$manifest = Get-Content -LiteralPath (Join-Path $packet 'manifest.json') -Raw | ConvertFrom-Json
foreach ($entry in $manifest.files_sha256.PSObject.Properties) {
    $got = (Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $packet $entry.Name)).Hash.ToLowerInvariant()
    if ($got -ne $entry.Value) { throw "Packet hash differs: $($entry.Name)" }
}
Push-Location $DnaCalcCheckout
try {
    $batches = Get-Content -LiteralPath (Join-Path $packet 'batches.json') -Raw | ConvertFrom-Json
    foreach ($batch in $batches) {
        & cargo xtask witness read --func $batch.function --inputs (Join-Path $packet $batch.path)
        if ($LASTEXITCODE -ne 0) { throw "Reading failed: $($batch.path). Preserve the raw worker output." }
    }
} finally { Pop-Location }
''')
    (out/'RUN-R15.txt').write_text('''R15 annuity structural audit; nine witness-reader batches.
Run RUN-R15.ps1 -DnaCalcCheckout <your master checkout>. Use the existing witness reader;
retain requested/submitted/read-back inputs, calculation settings, environment, raw
outputs and failed/unread obligations. Do not use --new-only: repeated controls matter.
Requests contain arguments only. Score immutable application predictions first;
they are never Excel answers. Then check conditional-relations.json.
The four anchor loans were fixed before this packet: two ordinary R11 loans and
the first tiny-rate 360-period and one-period R12 groups. Cashflow scalings are
exact factors 1, 2 and -1. Inverse-call payment arguments are original recorded
PMT publications, with their source IDs retained; no private helper is inferred.
Questions compare singleton cumulative/standalone surfaces, common unrounded
cashflows, unit growth/annuity calls, and explicit RATE guesses. A nonempty
necessary domain or equal component result never identifies a complete program.
All prior observations are discovery. There is no surviving complete candidate
to validate and no new Excel answer in this packet. Preserve the original R13/R14 packets.
''')
    manifest=dict(calls=sum(map(len,batches.values())),new_excel_answers=0,
        calls_by_function={f:len(batches[f])for f in functions},anchors=anchors,
        source_trace_sha256=hashlib.sha256(trace.read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        application_executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        files_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in out.iterdir()if p.is_file()},
        limitation='Unanswered structural questions; no complete candidate, selected helper or fresh validation claim.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))


if __name__=='__main__':
    if sys.argv[1]=='balance-sums':balance_sum_questions(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=='long-finish':long_finish_questions(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=='paired-finish':paired_finish(*sys.argv[2:]);raise SystemExit(0)
    if sys.argv[1]=='interior-backward-sum-domains':recurrence_sum_domains(*sys.argv[2:],backward=True,interior=True)
    elif sys.argv[1]=='interior-sum-domains':recurrence_sum_domains(*sys.argv[2:],interior=True)
    elif sys.argv[1]=='backward-sum-domains':recurrence_sum_domains(*sys.argv[2:],backward=True)
    elif sys.argv[1]=='affine-sum-domains':recurrence_sum_domains(*sys.argv[2:],balance=True)
    elif sys.argv[1]=='recurrence-sum-domains':recurrence_sum_domains(*sys.argv[2:])
    elif sys.argv[1]=='stored-component-sums':stored_component_sums(*sys.argv[2:])
    elif sys.argv[1]=='sum-factor-constraints':sum_factor_constraints(*sys.argv[2:])
    elif sys.argv[1]=='paired-sum-inverse':paired_sum_inverse(*sys.argv[2:])
    elif sys.argv[1]=='summation-analyze':summation_analyze(*sys.argv[2:])
    elif sys.argv[1]=='summation-questions':summation_questions(*sys.argv[2:])
    elif sys.argv[1]=='direct':direct_run(*sys.argv[2:])
    elif sys.argv[1]=='questions':questions(*sys.argv[2:])
    elif sys.argv[1]=='questions-r15':questions_r15(*sys.argv[2:])
    elif sys.argv[1]=='coupled':coupled(*sys.argv[2:])
    elif sys.argv[1]=='balance-domains':balance_domains(*sys.argv[2:])
    else:run(*sys.argv[1:])
