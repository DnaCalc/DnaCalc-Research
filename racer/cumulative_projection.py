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

from intervals import B64, E64, E53, Interval, inverse_op, inverse_add, log2floor, number, pow2

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
    if sys.argv[1]=='direct':direct_run(*sys.argv[2:])
    elif sys.argv[1]=='questions':questions(*sys.argv[2:])
    elif sys.argv[1]=='questions-r15':questions_r15(*sys.argv[2:])
    elif sys.argv[1]=='coupled':coupled(*sys.argv[2:])
    elif sys.argv[1]=='balance-domains':balance_domains(*sys.argv[2:])
    else:run(*sys.argv[1:])
