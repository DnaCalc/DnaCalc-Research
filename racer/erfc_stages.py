"""Invert the upper ERFC product, and pair the normal wrapper with public ERFC.

Original MIT tooling. No rational coefficients are fitted, and no match count
chooses a numerical program. Scalar receipts include every original readback.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import subprocess

from fingerprint import CachedService,load_service,Native,bits,value,sha,write
from exact import BINARY64 as B,Span,binary64,inverse
from x87 import load64

SOURCE_BINDINGS={str(p.resolve()):sha(p) for p in (
    Path(__file__),Path(__file__).with_name('x87_programs.py'),
    Path(__file__).with_name('x87_native.c'),
    *(Path(__file__).parents[1]/'checker'/name for name in ('geomean.py','exact.py','check.py','x87.py')))}


def complements(checkout,study,output,cache,library,offline=False):
    """Test the upper and negative wrappers against independent public ERFC."""
    from geomean import verified_study_answers
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete upper-ERFC study')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    frozen={**SOURCE_BINDINGS,**bindings,str(source.resolve()):sha(source),
            str(Path(library).resolve()):sha(library)}
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    xs=[row['input'] for row in doc['records']]
    erf=service.ask('ERF.PRECISE',[[client.number(sign*value(x))] for x in xs for sign in (1.,-1.)])['rows']
    negative=service.ask('ERFC.PRECISE',[[client.number(-value(x))] for x in xs])['rows']
    counts=Counter();records=[]
    for i,row in enumerate(doc['records']):
        q=row['erfc'];x=row['input']
        if q['type']!='number':raise ValueError('upper ERFC must be numeric')
        qx=load64(int(q['bits'],16))
        strict=bits(1.-value(q['bits']))
        stored=f'{native.publish(native.leaf("sub",load64(0x3ff0000000000000),qx,0x133f)):016x}'
        neg_strict=bits(2.-value(q['bits']))
        neg_stored=f'{native.publish(native.leaf("sub",load64(0x4000000000000000),qx,0x133f)):016x}'
        predictions=dict(erf_strict=strict,erf_stored_pc64=stored,
                         negative_erfc_strict=neg_strict,negative_erfc_stored_pc64=neg_stored)
        counts['positive_inputs']+=1
        for name,prediction in predictions.items():
            answer=erf[2*i] if name.startswith('erf_') else negative[i]
            counts[name+'_match']+=answer==dict(type='number',bits=prediction)
        counts['erf_odd_match']+=erf[2*i+1]==dict(type='number',bits=bits(-value(erf[2*i]['bits'])))
        records.append(dict(input=x,public_erfc=q,positive_erf=erf[2*i],negative_erf=erf[2*i+1],
                            negative_erfc=negative[i],predictions=predictions))
    result=dict(format='dnacalc-erfc-public-complements/1',complete=True,counts=counts,
        records=records,receipts=doc['receipts']+service.receipts,source_bindings=frozen,
        limits='Original upper ERFC answers, independent new positive/negative ERF and negative ERFC calls. '
        'Strict binary64 and stored-PC64 complement predictions are fixed arithmetic programs. '
        'An observed equality identifies a functional wrapper only on this grid; F remains unidentified '
        'and no private polynomial or universal ABI is inferred.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def complement_edges(checkout,output,cache,library,offline=False,count=256):
    """Preselect complement ties from their exponent, then pair public calls."""
    frozen={**SOURCE_BINDINGS,str(Path(library).resolve()):sha(library)}
    sites=[]
    for exponent in (53,54):
        low,high=.5,8.
        target=2.**-exponent
        for _ in range(80):
            middle=(low+high)/2.
            if math.erfc(middle)>target:low=middle
            else:high=middle
        center=(low+high)/2.;first=int(bits(center),16)-count//2
        sites.append(dict(exponent=exponent,center=bits(center),
                          inputs=[f'{first+k:016x}' for k in range(count)]))
    first=int(bits(.5),16)-16
    sites.append(dict(exponent=None,center=bits(.5),inputs=[f'{first+k:016x}' for k in range(33)]))
    rows=[dict(site=i,input=x) for i,s in enumerate(sites) for x in s['inputs']]
    write(str(output)+'.questions.json',dict(format='dnacalc-erfc-complement-edge-questions/1',complete=True,
        sites=sites,source_bindings=frozen,
        selection='Two fixed complement midpoint powers,2^-53 and2^-54, located with host erfc before Excel; adjacent-double windows and33 controls at0.5.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    qs=service.ask('ERFC.PRECISE',[[client.number(r['input'])] for r in rows])['rows']
    predictions=[]
    for row,q in zip(rows,qs):
        if q['type']!='number':raise ValueError('finite positive ERFC must be numeric')
        raw=load64(int(q['bits'],16));one=load64(0x3ff0000000000000);two=load64(0x4000000000000000)
        first=f'{native.publish(native.leaf("sub",one,raw,0x133f)):016x}'
        double=f'{native.publish(native.leaf("sub",two,raw,0x133f)):016x}'
        nested=f'{native.publish(native.leaf("add",one,load64(int(first,16)),0x133f)):016x}'
        predictions.append(dict(row,public_erfc=q,positive_erf_strict=bits(1.-value(q['bits'])),
                                positive_erf_stored_pc64=first,negative_erfc_strict=bits(2.-value(q['bits'])),
                                negative_erfc_stored_pc64=double,negative_erfc_two_stores=nested))
    # Retain predictions before observing the wrappers.
    write(str(output)+'.predictions.json',dict(format='dnacalc-erfc-complement-edge-predictions/1',complete=True,
        records=predictions,receipts=service.receipts,source_bindings=frozen))
    observations={}
    for name in ('ERF.PRECISE','ERF','ERFC.PRECISE','ERFC'):
        sign=1. if name.startswith('ERF') and not name.startswith('ERFC') else -1.
        observations[name]=service.ask(name,[[client.number(sign*value(r['input']))] for r in rows])['rows']
    counts=Counter()
    for i,row in enumerate(predictions):
        row['answers']={name:a[i] for name,a in observations.items()}
        counts['rows']+=1
        for name,answer in row['answers'].items():
            family='negative_erfc_' if name.startswith('ERFC') else 'positive_erf_'
            for model in ('strict','stored_pc64')+ (('two_stores',) if family=='negative_erfc_' else ()):
                counts[name+'/'+model]+=answer==dict(type='number',bits=row[family+model])
    result=dict(format='dnacalc-erfc-complement-exponent-edges/1',complete=True,counts=counts,
        records=predictions,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),predictions_sha256=sha(str(output)+'.predictions.json'),
        limits='One complement stage. Sites come from exponent midpoints, not current match counts. '
        'Original public ERFC leaves determine every prediction before new wrapper answers arrive. '
        'The0.5 controls include the neighboring lower branch. No private F table or correction is inferred.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def classify_losses(bank,bank_source,executable,checkout,output,cache,library):
    """Bind every new lost bank match to its original public ERFC counterpart."""
    import gzip
    bank,bank_source,executable,checkout,output,library=map(Path,(bank,bank_source,executable,checkout,output,library))
    changed=bank/'changed-outcomes.jsonl.gz'
    frozen={**SOURCE_BINDINGS,str(changed.resolve()):sha(changed),
            str(executable.resolve()):sha(executable),str(library.resolve()):sha(library)}
    losses=[json.loads(line) for line in gzip.open(changed,'rt') if json.loads(line)['before_matched']]
    arguments={}
    for row in losses:
        args=[value(x.removeprefix('0x')) for x in row['args']]
        if row['function'] in ('ERF','ERF.PRECISE') and len(args)==1:
            x=args[0]
        elif row['function'] in ('GAMMA.DIST','GAMMADIST') and len(args)==4 and args[1]==.5 and args[3]!=0:
            x=math.sqrt(args[0]/args[2])
        else:
            raise ValueError('unclassified lost function/shape: '+row['function'])
        if not math.isfinite(x) or abs(x)<=.5:raise ValueError('loss outside changed upper complement')
        arguments[bits(abs(x))]=None
        row['upper_erf_argument']=bits(abs(x));row['negative']=x<0
    counterpart=bank_source/'cache/oracle/build-20131/ERFC.PRECISE.jsonl.gz'
    frozen[str(counterpart.resolve())]=sha(counterpart)
    for line in gzip.open(counterpart,'rt'):
        raw=json.loads(line);args=raw['args']
        if len(args)==1 and len(args[0])==1:
            arg=args[0][0]
            if arg.get('kind')=='number':
                b=arg['bits_hex'].removeprefix('0x').lower()
                if b in arguments:arguments[b]=dict(original=raw,answer=dict(type='number',bits=raw['outcome']['bits_hex'].removeprefix('0x').lower()))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache)
    missing=[b for b,v in arguments.items() if v is None]
    if missing:
        answers=service.ask('ERFC.PRECISE',[[client.number(b)] for b in missing])['rows']
        for b,a in zip(missing,answers):arguments[b]=dict(original=None,answer=a)
    inputs=Path(str(output)+'.inputs.json');actuals=Path(str(output)+'.actuals.jsonl')
    witnesses=[dict(id='upper-'+str(i),args=['0x'+b]) for i,b in enumerate(arguments)]
    write(inputs,dict(function='ERFC.PRECISE',witnesses=witnesses))
    # ubs:ignore — fixed argv to the bound local application executable; no shell, bounded execution.
    subprocess.run([str(executable),'invoke-inputs',str(inputs),str(actuals)],check=True,timeout=300)
    computed=[json.loads(line) for line in actuals.read_text().splitlines()]
    if len(computed)!=len(witnesses):raise ValueError('incomplete public ERFC replay')
    native=Native(library)
    for (b,entry),wanted,actual in zip(arguments.items(),witnesses,computed):
        if any(actual[k]!=wanted[k] for k in ('id','args')):raise ValueError('upper ERFC replay identity changed')
        expected=entry['answer']
        if expected['type']!='number':raise ValueError('upper ERFC counterpart is not numeric')
        if not actual['calc_func_prediction'].startswith('0x'):raise ValueError('upper ERFC product did not calculate')
        calculated=actual['calc_func_prediction'].removeprefix('0x')
        def complement(q):
            return f'{native.publish(native.leaf("sub",load64(0x3ff0000000000000),load64(int(q,16)),0x133f)):016x}'
        entry.update(calculated=calculated,original_complement=complement(expected['bits']),
                     calculated_complement=complement(calculated))
    for row in losses:
        entry=arguments[row['upper_erf_argument']]
        sign=(1<<63) if row['negative'] else 0
        old=f'0x{int(entry["original_complement"],16)^sign:016x}'
        new=f'0x{int(entry["calculated_complement"],16)^sign:016x}'
        if row['expected']!=old or row['actual']!=new or entry['calculated']==entry['answer']['bits']:
            raise ValueError('lost row is not explained by its upper ERFC counterpart')
        if entry['original'] is not None:
            environment=entry['original']['environment']
            if str(environment['excel_build'])!=row['reading']['excel_build'] or environment['cpu_id']!=row['reading']['cpu']:
                raise ValueError('original upper counterpart context differs')
        row['classification']='unchanged_positive_upper_ERFC_body_gap_exposed_by_measured_complement'
    result=dict(format='dnacalc-erf-complement-loss-classification/1',complete=True,
        lost_bank_rows=len(losses),unique_upper_arguments=len(arguments),original_bank_counterparts=len(arguments)-len(missing),
        losses=losses,upper_counterparts=arguments,receipts=service.receipts,source_bindings=frozen,
        application_replay=dict(inputs_sha256=sha(inputs),actuals_sha256=sha(actuals)),
        limits='Every lost original bank match is retained with its input, expected and calculated bits. '
        'The old expected result is the stored-PC64 complement of the original public ERFC reading; '
        'the new result is the same complement of the unchanged product upper body. This locates '
        'the exposed discrepancy in that stand-in body, without claiming its private F table, '
        'exponential ABI or underlying coefficient cause is identified. Fresh counterparts, if '
        'needed, retain their own environment and never replace the old answer.')
    write(output,result);print(json.dumps({k:result[k] for k in ('lost_bank_rows','unique_upper_arguments','original_bank_counterparts')}),flush=True)
    return result


def run(checkout,output,cache,library,count=256,offline=False):
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    centers=(.5001,.7,1.,1.11567925,1.1217,1.5,2.5,3.,4.,5.5,8.,16.,26.)
    xs=[f'{int(bits(x),16)-count//2+k:016x}' for x in centers for k in range(count)]
    answers=service.ask('ERFC.PRECISE',[[client.number(x)] for x in xs])['rows']
    square={x:f'{native.publish(native.leaf("mul",load64(int(x,16)),load64(int(x,16)),0x133f)):016x}' for x in xs}
    exponents=sorted({bits(-value(z)) for z in square.values()})
    exponential=dict(zip(exponents,service.ask('EXP',[[client.number(x)] for x in exponents])['rows']))
    counts=Counter();records=[];span_hash=hashlib.sha256()
    for x,answer in zip(xs,answers):
        e=exponential[bits(-value(square[x]))]
        if answer['type']!=e['type'] or answer['type']!='number':
            raise ValueError('normal finite upper ERFC/EXP must be numeric')
        target=binary64(answer['bits']);factor=binary64(e['bits'])
        domain=inverse(Span(target,target),factor,'mul',True,'stored_pc64')
        span_hash.update(json.dumps(domain.json(),sort_keys=True).encode()+b'\n')
        coordinates=B.coordinates(domain)
        n=0 if coordinates is None else coordinates[1]-coordinates[0]+1
        if n>4:
            raise ValueError('unexpectedly broad F preimage')
        row=dict(input=x,square=square[x],exponential=e,erfc=answer,binary64_F_count=n,
            binary64_F_values=[] if coordinates is None else [bits(float(B.value(k))) for k in range(coordinates[0],coordinates[1]+1)])
        if not n:
            row['adjacent_F_bits']=[bits(float(B.value(k))) for k in
                (B.floor_coordinate(domain.low),B.ceiling_coordinate(domain.high))]
        counts['rows']+=1;counts['binary64_F_'+str(n)]+=1
        records.append(row)
    # Independent public wrapper calls use their own staged z, not an assumed
    # exact inverse of multiplying an ERFC argument by sqrt(2).
    ys=[f'{int(bits(y),16)-64+k:016x}' for y in (.9,1.,2.,4.,8.,16.,32.,37.5) for k in range(128)]
    c=load64(0x3fe6a09e667f3bcd)
    zs={y:f'{native.publish(native.leaf("mul",load64(int(y,16)),c,0x133f)):016x}' for y in ys}
    distinct_z=sorted(set(zs.values()))
    public_erfc=dict(zip(distinct_z,service.ask('ERFC.PRECISE',[[client.number(z)] for z in distinct_z])['rows']))
    normals=service.ask('NORM.S.DIST',[[client.number(sign*value(y)),client.boolean(True)] for y in ys for sign in (-1.,1.)])['rows']
    wrappers=[]
    for i,y in enumerate(ys):
        q=public_erfc[zs[y]]
        if q['type']!='number':
            raise ValueError('normal wrapper ERFC control must be numeric')
        half=.5*value(q['bits'])
        if half<2.**-1022:
            half=0.
        for j,sign in enumerate((-1.,1.)):
            prediction=bits(half if sign<0 else 1.-half)
            row=dict(input=bits(sign*value(y)),z=zs[y],public_erfc=q,
                     answer=normals[2*i+j],prediction=prediction)
            counts['normal_rows']+=1
            counts['normal_wrapper_matched']+=row['answer']==dict(type='number',bits=prediction)
            wrappers.append(row)
    result=dict(format='dnacalc-erfc-product-preimages/1',complete=True,counts=counts,records=records,
        normal_wrappers=wrappers,receipts=service.receipts,F_domain_sha256=span_hash.hexdigest(),
        source_bindings={str(p):sha(p) for p in (Path(__file__),Path(library))},
        selection='Thirteen fixed upper-ERFC sites, adjacent double windows; eight fixed normal-wrapper windows with both signs.',
        limits='The F inverse conditions on z=storedPC64(x*x), independently measured public EXP(-z), and a final storedPC64 E*F. Empty binary64 preimages refute a stored-F composition under these premises; they do not identify a rational, coefficient, ABI, or the alternative association E*N/D. Normal wrapper predictions use separately observed ERFC values at staged z. No numerical kernel is selected by counts.')
    write(output,result);print(json.dumps(counts,indent=2),flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout');p.add_argument('output')
    p.add_argument('--cache',required=True);p.add_argument('--library',required=True)
    p.add_argument('--count',type=int,default=256);p.add_argument('--offline',action='store_true')
    p.add_argument('--complements')
    p.add_argument('--complement-edges',action='store_true')
    p.add_argument('--classify-losses',nargs=3,metavar=('BANK','BANK_SOURCE','EXECUTABLE'))
    a=p.parse_args()
    if not 2<=a.count<=512:
        raise ValueError('small windows only')
    if a.classify_losses:
        classify_losses(*a.classify_losses,a.checkout,a.output,a.cache,a.library)
    elif a.complement_edges:
        complement_edges(a.checkout,a.output,a.cache,a.library,a.offline,a.count)
    elif a.complements:
        complements(a.checkout,a.complements,a.output,a.cache,a.library,a.offline)
    else:
        run(a.checkout,a.output,a.cache,a.library,a.count,a.offline)
