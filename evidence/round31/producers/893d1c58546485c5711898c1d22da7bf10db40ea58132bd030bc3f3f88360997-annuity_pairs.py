"""Independent native replay of R29's isolated PMT arithmetic stages.

The scalar lineage and exact inverse checker are reused. No producer or
application module is imported. The machine executes each arithmetic leaf
directly, independently of the rational producer and its profile evaluation.
"""
from collections import Counter
import ctypes as C
import json
from pathlib import Path
import struct
import sys

from check import require,digest
from exact import BINARY64 as B,Span,binary64,inverse
from geomean import verified_study_answers
from x87 import load64,decode

SOURCE_BINDINGS={str(p.resolve()):digest(p) for p in (
    Path(__file__),*(Path(__file__).with_name(n) for n in ('exact.py','check.py','geomean.py','x87.py')),
    Path(__file__).parents[1]/'racer/x87_native.c')}


class Register(C.Structure):
    _fields_=[('bytes',C.c_ubyte*10)]


class Arithmetic:
    def __init__(self,library):
        self.lib=C.CDLL(str(Path(library).resolve()))
        ptr=C.POINTER(Register)
        self.lib.x87_leaf.argtypes=[C.c_int,ptr,ptr,C.c_uint16,ptr]
        self.lib.x87_leaf.restype=None
        self.lib.x87_publish.argtypes=[ptr,C.c_uint16]
        self.lib.x87_publish.restype=C.c_uint64
        self.counts=Counter()

    def op(self,kind,a,b):
        aa,bb=Register.from_buffer_copy(a),Register.from_buffer_copy(b)
        raw=Register()
        self.lib.x87_leaf({'add':2,'mul':4,'div':5}[kind],C.byref(aa),C.byref(bb),0x133f,C.byref(raw))
        result=self.lib.x87_publish(C.byref(raw),0x133f)
        self.counts[kind]+=1
        return load64(result)

    def bits(self,a):
        raw=Register.from_buffer_copy(a)
        return f'{self.lib.x87_publish(C.byref(raw),0x133f):016x}'

    def finish(self,a,r,model):
        if model=='strict':
            x,y=(struct.unpack('>d',bytes.fromhex(self.bits(v)))[0] for v in (a,r))
            return struct.pack('>d',x*y).hex()
        require(model=='stored_pc64','unsupported final model')
        return self.bits(self.op('mul',a,r))


def check(study,cache,library,output):
    study,cache,library,output=map(Path,(study,cache,library,output))
    doc=json.loads(study.read_text())
    require(doc.get('complete') is True,'incomplete stage study')
    answers,bindings,ordered=verified_study_answers(doc,cache)
    frozen={str(p):digest(p) for p in (study,library)};frozen.update(SOURCE_BINDINGS);frozen.update(bindings)
    a=Arithmetic(library);one=load64(0x3ff0000000000000)
    model=doc.get('final_model','stored_pc64')
    counts=Counter()
    def observed(args):
        args=[dict(type='number',bits=x.removeprefix('0x')) for x in args]
        return answers[('PMT',json.dumps(args,sort_keys=True))]
    def assert_answer(args,expected):
        got=observed(args)
        require(got==dict(type='number',bits=expected),'native isolated-stage replay differs from original Excel')
        counts['native_stage_matches']+=1
    def infer(rate,periods,observations):
        domain=None;r=binary64(rate)
        for row in observations:
            pv=row['pv'];zero=row['zero']
            require(observed([rate,periods,pv,'0000000000000000','0000000000000000'])==zero,'original type-zero reading changed')
            y=binary64(zero['bits'])
            required=inverse(inverse(Span(y,y),r,'mul',True,model),binary64(pv),'div',False,'stored_pc64')
            domain=required if domain is None else domain.meet(required)
            require(domain is not None,'no common denominator')
        coordinates=B.coordinates(domain)
        require(coordinates is not None and coordinates[0]==coordinates[1],'denominator not unique')
        em=B.value(coordinates[0]);counts['unique_denominators']+=1
        return struct.pack('>d',float(em)).hex()
    def finish(rate,pv,fv,em,timing):
        r,p,f,e=(load64(int(x,16)) for x in (rate,pv,fv,em))
        balance=a.op('add',p,f)
        q=a.op('div',balance,e)
        corrected=a.op('add',q,f)
        factor=r if not timing else a.op('div',one,a.op('add',one,a.op('div',one,r)))
        return a.finish(corrected,factor,model)
    if doc['format']=='dnacalc-annuity-reciprocal-domain/1':
        require(doc['phase']==2,'initial contradicted hypothesis is not an accepted domain rule')
        questions=Path(str(study)+'.questions.json')
        require(digest(questions)==doc['questions_sha256'],'domain question manifest changed')  # ubs:ignore — public evidence checksum, not a secret or authentication token.
        frozen[str(questions)]=digest(questions)
        by_function={}
        poles={}
        for row in doc['records']:
            name,args=row['function'],row['args']
            request=[dict(type='number',bits=x)for x in args]
            original=answers[(name,json.dumps(request,sort_keys=True))]
            require(original==row['answer'],'original domain answer changed')
            rate=binary64(args[0]);timing=binary64(args[-1])
            if args[0] not in poles:
                denominator=a.op('add',one,a.op('div',one,load64(int(args[0],16)))) if rate else one
                poles[args[0]]=binary64(a.bits(denominator))==0
            pole=poles[args[0]]
            require(pole==row['reciprocal_pole'],'native and rational reciprocal poles differ')
            first=name=='IPMT'and binary64(args[1])==1 and timing==1
            rejected=rate<=-1 or (pole and timing==1)
            wanted='not_asserted'if name=='FV'else 'zero'if first else 'error:Num'if rejected else 'not_asserted'
            require(wanted==row['expected_admission'],'declared domain branch differs')
            seen='zero'if original['type']=='number'and binary64(original['bits'])==0 else 'error:Num'if original==dict(type='error',code='#NUM!')else 'other'
            require(seen==row['observed_admission'],'observed domain classification differs')
            tally=by_function.setdefault(name,Counter())
            tally['calls']+=1;tally['constrained']+=wanted!='not_asserted'
            tally['contradictions']+=wanted!='not_asserted'and wanted!=seen
            tally['pole_calls']+=pole
        require(by_function==doc['counts'],'domain counts differ')
        require(all(c['contradictions']==0 for c in by_function.values()),'public domain contradiction')
        counts.update({name+'/'+field:n for name,c in by_function.items()for field,n in c.items()})
    elif doc['format']=='dnacalc-annuity-paired-stages/1':
        for group in doc['groups']:
            r,n=group['rate'],group['periods']
            observations=[dict(pv=row['pv'],zero=row['zero']) for row in group['timing']]
            em=infer(r,n,observations)
            require(group['em_values']==[em],'reported constrained denominator differs')
            for row in group['timing']:
                for timing in (0,1):
                    expected=finish(r,row['pv'],'0000000000000000',em,timing)
                    assert_answer([r,n,row['pv'],'0000000000000000',f'{0x3ff0000000000000 if timing else 0:016x}'],expected)
            for row in group['numerator']:
                expected=finish(r,row['pv'],row['fv'],em,0)
                assert_answer([r,n,row['pv'],row['fv'],'0000000000000000'],expected)
    elif doc['format']=='dnacalc-annuity-stage-discriminants/1':
        questions=Path(str(study)+'.questions.json')
        require(digest(questions)==doc['questions_sha256'],'pre-answer manifest changed')  # ubs:ignore — public file hashes, not secrets or authentication tokens.
        frozen[str(questions)]=digest(questions)
        for row in doc['sums']:
            r,n,p,f,t=row['args']
            require(t=='0000000000000000','sum discriminant changed timing')
            expected=finish(r,p,f,row['em'],0)
            require(expected==row['predictions']['stored_pc64/stored_pc64'],'native and rational sum-stage predictions differ')
            assert_answer(row['args'],expected)
            counts[row['stage']]+=1
        for group in doc['rates']:
            r=group['rate'];n='4028000000000000'
            em=infer(r,n,group['observations'])
            require(em==group['em'],'new-rate denominator changed')
            for row in group['observations']:
                for timing in (0,1):
                    expected=finish(r,row['pv'],'0000000000000000',em,timing)
                    assert_answer([r,n,row['pv'],'0000000000000000',f'{0x3ff0000000000000 if timing else 0:016x}'],expected)
    else:
        raise ValueError('unknown study')
    require(frozen=={p:digest(p) for p in frozen},'checker inputs changed')
    result=dict(complete=True,counts=counts,native_operations=a.counts,inputs_sha256=frozen,
        original_calls=len(ordered),
        limits='Independent direct-instruction arithmetic validates each declared isolated stage against original scalar readings. Denominators remain exact necessary constraints, not adopted helper values. Actual application results and bank/corpus gains require their separate registered-route replays.')
    with output.open('x') as f:
        json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(dict(counts=counts,native_operations=a.counts)),flush=True)
    return result


if __name__=='__main__':
    check(*sys.argv[1:])
