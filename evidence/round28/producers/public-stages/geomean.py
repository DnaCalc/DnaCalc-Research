"""Independent ISA grammar and rational arithmetic checks for GEOMEAN.

Only the earlier independent checker primitives are reused. No search module
or product source is imported. FYL2X/F2XM1/ROM leaves use a separate native
single-instruction replay; they are measurements, not mathematical substitutes.
"""
from collections import Counter
from fractions import Fraction as Q
import itertools
import json
from pathlib import Path
import sys
import time

from exact import binary64,two
from check import records,require,digest
from x87 import Exact,Instructions,published_bits

WORDS=(0x123f,0x133f)


def coordinates(row):
    if row['function']!='GEOMEAN'or row['result_index']is not None:return None,'surface'
    if 'reference_binding' in row:
        verify_static_reference(row)
    def flatten(v):
        if not isinstance(v,list):return[v]
        return list(itertools.chain.from_iterable(flatten(a)for a in v))
    args=flatten(row['args'])
    if not all(isinstance(x,str)and x.startswith('0x')and len(x)==18 for x in args):return None,'non_numeric_argument'
    values=list(map(binary64,args))
    if not(2<=len(values)<=32)or any(x<=0 for x in values):return None,'length_or_domain'
    def product(xs,order):
        def multiply(a,b):
            z=a*b
            require(two(-1000)<=z<=two(1000),'partial_product_bound')
            return z
        if order==2:
            if not xs:return Q(1)
            if len(xs)==1:return xs[0]
            middle=len(xs)//2
            return multiply(product(xs[:middle],2),product(xs[middle:],2))
        result=Q(1)
        for x in xs if order==0 else xs[::-1]:result=multiply(result,x)
        return result
    try:
        for order in range(3):product(values,order)
    except ValueError:return None,'partial_product_bound'
    if not row['expected'].startswith('0x'):return None,'non_numeric_expected'
    return [int(x,16)for x in args],None


def verify_static_reference(row):
    """Re-read the immutable package/capture, rather than trust derived operands."""
    import zipfile
    import xml.etree.ElementTree as ET
    binding=row['reference_binding'];package=Path(binding['package']);capture=Path(row['source'])
    require(digest(package)==binding['package_sha256']and digest(capture)==binding['capture_sha256'],'reference source changed') # ubs:ignore — public file digests, no authentication secret or timing boundary.
    require(row['formula'].lstrip('=')in('GEOMEAN(Data!A1:A3)','_xlfn.GEOMEAN(Data!A1:A3)')and
            binding['sheet']=='Data'and binding['cells']==['A1','A2','A3'],'unsupported reference binding')
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'};rid='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
    with zipfile.ZipFile(package)as z:
        rels={r.attrib['Id']:r.attrib['Target']for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheet=next(s for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets',ns)if s.attrib['name']=='Data')
        cells={c.attrib['r']:c for c in ET.fromstring(z.read('xl/'+rels[sheet.attrib[rid]])).findall('.//m:c',ns)}
        values=[]
        for name in('A1','A2','A3'):
            cell=cells[name];require(cell.find('m:f',ns)is None and cell.attrib.get('t','n')=='n','not a static numeric reference')
            value=Q(cell.find('m:v',ns).text)
            require(value==binary64(binding['values'][len(values)]),'reference numeric value changed')
            values.append(binding['values'][len(values)])
    require(row['args']==[[[v]for v in values]],'reference operands differ')
    doc=json.loads(capture.read_text())
    package_source='sha256:'+binding['package_sha256']
    require(doc['source']==package_source and doc['session']['packages']['input']==package_source,
            'reference capture/package lineage differs')
    f=row['original_formula'];key=(f['sheet'],f['row'],f['column'])
    require(f in doc['formulas']and f['formula']==row['formula'],'reference formula changed')
    original=next(r for r in doc['observations']if(r['sheet'],r['row'],r['column'])==key)
    require(original==row['original_capture']and row['reading']==doc['environment'],'reference observation changed')
    require(row['expected']=='0x'+original['excel_raw']['number'],'reference answer changed')


def grammar(count,path,order):
    nodes=[dict(op='input',a=i)for i in range(count)]
    def const(bits):nodes.append(dict(op='const',bits=bits));return len(nodes)-1
    one=const('3ff0000000000000');zero=const('0000000000000000')
    # A count <=32 is an exact binary64 integer; construct its encoding by
    # integer arithmetic rather than the search's host-float conversion.
    exponent=count.bit_length()-1
    n=const(f'{((exponent+1023)<<52)|((count-(1<<exponent))<<(52-exponent)):016x}')
    phase=0
    def append(kind,a=0,b=0,site=None):
        node=dict(op=kind,a=a,b=b)
        if phase:node['phase']=phase
        if site is not None:node['site']=site
        nodes.append(node);return len(nodes)-1
    def reduce(items,kind,identity):
        if order==2:
            if not items:return identity
            if len(items)==1:return items[0]
            middle=len(items)//2
            return append(kind,reduce(items[:middle],kind,identity),reduce(items[middle:],kind,identity),0)
        result=identity
        for item in items if order==0 else items[::-1]:result=append(kind,result,item,0)
        return result
    if path in (0,1):
        reduced=reduce(list(range(count)),'mul',one);nodes[reduced]['site']=1
        if path==0:
            log=append('fyl2x',one,reduced,2);argument=append('div',log,n,3)
        else:
            coefficient=append('div',one,n,3);argument=append('fyl2x',coefficient,reduced,2)
    else:
        coefficient=one if path==2 else append('fldln2')
        logs=[append('fyl2x',coefficient,i,2)for i in range(count)]
        reduced=reduce(logs,'add',zero);nodes[reduced]['site']=1
        argument=append('div',reduced,n,3)
        if path==3:argument=append('mul',argument,append('fldl2e'),3)
    phase=1
    integer=append('frndint',argument)
    fraction=append('sub',argument,integer,4)
    small=append('f2xm1',fraction,0,5)
    mantissa=append('add',one,small,6)
    result=append('fscale',mantissa,integer,7)
    require(len(nodes)<=128,'native interpreter bound')
    return dict(nodes=nodes,outputs=[result,result],sites=8)


def check(input_path,directory,library,output):
    start=time.monotonic();source=Path(__file__);directory=Path(directory)
    hashes={str(p):digest(p)for p in (source,source.with_name('x87.py'),source.with_name('exact.py'),source.with_name('x87_isa.c'))}
    original=list(records(input_path));scoped=[];excluded=[];scope=Counter()
    for index,raw in enumerate(original):
        inputs,reason=coordinates(raw);row=dict(trace_index=index,raw=raw)
        if reason:row['reason']=reason;excluded.append(row);scope[reason]+=1
        else:row['inputs']=inputs;scoped.append(row);scope['scoped']+=1
    require(scoped==list(records(directory/'scoped.jsonl'))and excluded==list(records(directory/'excluded.jsonl')),'changed input/scope')
    ps=[dict(id=i,path=p,order=o,cw=c,inner_cw=w,mask=m)for i,(p,o,c,w,m)in enumerate(itertools.product(range(4),range(3),WORDS,WORDS,range(256)))]
    require(ps==list(records(directory/'profiles.jsonl'))and len(ps)==12288,'incomplete grammar')
    summary=json.loads((directory/'summary.json').read_text());require(summary['complete']is True,'partial run')
    instructions=Instructions(library);exact=Exact(instructions)
    killed=set();counts=Counter();executions=0;per_witness=Counter()
    def replay(identity,order,traced):
        p=ps[identity];row=scoped[order]
        g=grammar(len(row['inputs']),p['path'],p['order'])
        # Storing each reduction result also stores the final reduction node.
        mask=p['mask']|(2 if p['mask']&1 else 0)
        return exact.graph(g,row['inputs'],p['cw'],mask,inner_cw=p['inner_cw'],trace=traced)
    for record in records(directory/'first-kills.jsonl'):
        identity,order=record['candidate_id'],record['order']
        require(identity not in killed and 0<=identity<len(ps)and 0<=order<len(scoped),'duplicate/unknown program or witness')
        require(record['trace_index']==scoped[order]['trace_index']and record['expected']==scoped[order]['raw']['expected'],'changed observation')
        for previous in range(order):
            root=replay(identity,previous,False)[:10];executions+=1
            require(f'0x{published_bits(root):016x}'==scoped[previous]['raw']['expected'],'not first contradiction')
        roots,raw,values=replay(identity,order,True);executions+=1
        root=roots[:10];got=f'0x{published_bits(root):016x}'
        require(root.hex()==record['root']and raw.hex()==record['raw']and values.hex()==record['values'],'full-width instruction trace differs')
        require(got==record['got']and got!=record['expected'],'invalid first kill')
        counts['full_width_instruction_values']+=len(values)//10
        counts['certified_first_kills']+=1;killed.add(identity);per_witness[order]+=1
    alive=json.loads((directory/'survivors.json').read_text())
    require(len(set(alive))==len(alive)and killed.isdisjoint(alive)and killed|set(alive)==set(range(len(ps))),'unsearched program')
    for identity in alive:
        for order,row in enumerate(scoped):
            root=replay(identity,order,False)[:10];executions+=1
            require(f'0x{published_bits(root):016x}'==row['raw']['expected'],'invalid survivor')
    remaining=len(ps);witnesses=list(records(directory/'witnesses.jsonl'))
    require(len(witnesses)==len(scoped),'missing witness accounting')
    for order,(record,row)in enumerate(zip(witnesses,scoped)):
        remaining-=per_witness[order]
        require(record==dict(order=order,trace_index=row['trace_index'],killed=per_witness[order],remaining=remaining),'incorrect witness accounting')
    require(summary['programs']==len(ps)and summary['recorded']==len(original)and summary['scoped']==len(scoped)and
            summary['scope']==scope and summary['survivors']==len(alive)and summary['executions']==executions,'incorrect summary')
    require(all(digest(p)==h for p,h in hashes.items()),'checker changed during run')  # ubs:ignore — public integrity hashes, not secrets
    transcript=Path(str(output)+'.native-leaves.jsonl')
    with transcript.open('x')as out:
        for (kind,a,b,cw),value in sorted(instructions.transcript.items()):
            out.write(json.dumps(dict(kind=kind,a=a.hex(),b=b.hex(),cw=cw,value=value.hex()))+'\n')
    result=dict(complete=True,programs=len(ps),scoped=len(scoped),survivors=len(alive),executions=executions,
                counts=counts,rational_operations=exact.counts,native_leaves=instructions.counts,
                seconds=time.monotonic()-start,source_sha256=hashes,native_replay_sha256=digest(library),
                input_sha256={str(p):digest(p)for p in (Path(input_path),transcript,*sorted(directory.iterdir()))if p.is_file()},
                limits='Declared x87 grammar on the recorded native CPU only. Microcode leaves are replayed instructions, not correctly-rounded mathematical functions. No general Windows/CPU equivalence or complete Excel kernel identification.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def check_exact_roots(input_path, directory, output):
    """Verify roots by inverse midpoint powers, without the producer's root search."""
    from exact import operation
    directory = Path(directory)
    frozen = {str(p):digest(p) for p in (Path(input_path), *sorted(directory.iterdir()),
              Path(__file__), Path(__file__).with_name('exact.py')) if p.is_file()}
    original = list(records(input_path)); rows = []; scope = Counter()
    for index, raw in enumerate(original):
        inputs, reason = coordinates(raw)
        if reason: scope[reason] += 1
        else: rows.append(dict(trace_index=index, raw=raw, inputs=inputs)); scope['scoped'] += 1
    require(rows, 'empty exact-root scope')
    require(json.loads((directory/'scoped.json').read_text()) == rows, 'changed exact-root scope')
    modes = ('strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53','exact')
    ps = [dict(id=i, order=order, model=model) for i,(order,model) in enumerate(itertools.product(range(3),modes))]
    require(json.loads((directory/'profiles.json').read_text()) == ps, 'changed exact-root grammar')
    def product(profile, row):
        xs = list(map(binary64, row['inputs']))
        def multiply(a,b): return a*b if profile['model']=='exact' else operation(a,b,'mul',profile['model'])
        def balanced(items):
            if len(items)==1:return items[0]
            k=len(items)//2
            return multiply(balanced(items[:k]),balanced(items[k:]))
        if profile['order']==2:return balanced(xs)
        result=Q(1)
        for x in xs if profile['order']==0 else xs[::-1]:result=multiply(result,x)
        return result
    def verifies(value, count, encoded):
        bits=int(encoded,16)
        require(0 < bits < 0x7ff0000000000000, 'nonpositive or nonfinite root')
        center=binary64(bits);low=((binary64(bits-1)+center)/2)**count;high=((binary64(bits+1)+center)/2)**count
        return (low <= value <= high) if bits%2==0 else (low < value < high)
    seen=set();executions=0
    for kill in records(directory/'first-kills.jsonl'):
        identity, at=kill['candidate_id'], kill['order']
        require(identity not in seen and 0<=identity<len(ps) and 0<=at<len(rows), 'bad root kill identity')
        for prior in range(at+1):
            row=rows[prior];value=product(ps[identity],row);executions+=1
            if prior<at:require(verifies(value,len(row['inputs']),row['raw']['expected']), 'not first root contradiction')
            else:
                require(value==Q(*map(int,kill['product'])) and
                        verifies(value,len(row['inputs']),kill['got']), 'incorrect product/root')
                require(kill['expected']==row['raw']['expected']!=kill['got'], 'changed root target')
        seen.add(identity)
    survivors=json.loads((directory/'survivors.json').read_text())
    require(len(set(survivors))==len(survivors) and seen.isdisjoint(survivors) and
            seen|set(survivors)==set(range(len(ps))), 'unsearched root program')
    for identity in survivors:
        for row in rows:
            require(verifies(product(ps[identity],row),len(row['inputs']),row['raw']['expected']), 'invalid root survivor')
            executions+=1
    summary=json.loads((directory/'summary.json').read_text())
    require(summary['scope']==scope,'wrong public root scope') # ubs:ignore — a public numerical-study label, not an OAuth or authentication secret.
    require(summary['complete'] and summary['programs']==len(ps) and summary['survivors']==len(survivors)
            and summary['observations']==len(rows) and summary['executions']==executions,
            'wrong root summary')
    require(all(digest(p)==h for p,h in frozen.items()), 'root proof/source changed') # ubs:ignore — equality of public file digests, not an authentication secret.
    result=dict(complete=True,programs=len(ps),observations=len(rows),survivors=len(survivors),executions=executions,
                input_sha256=frozen,limits='Exact inverse midpoint-power verification of the specified rounded-product/correct-root programs. No native transcendental instruction, fitted answer or Excel internal precision is assumed.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def root_publication_domains(products, count, target, staged=False):
    from exact import BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P, Span
    binary = B.pullback(Span(target, target))
    formats = {'binary64': binary}
    if staged:
        formats['extended64_then_binary64'] = E.pullback(binary)
        formats['extended53_then_binary64'] = P.pullback(binary)
    results = {}
    for name, span in formats.items():
        require(span is not None and span.low > 0, 'positive normal root publication required')
        inverse_root = Span(span.low**count, span.high**count, span.left, span.right)
        results[name] = dict(span=span.json(), reachable=any(inverse_root.contains(v) for v in products))
    return results


def product_result(raw, staged=False):
    """Explicit ordered trees; independent of the producer's subset closure."""
    from exact import operation
    inputs, reason = coordinates(raw)
    require(reason is None and 2 <= len(inputs) <= 3, 'product-closure arity/domain')
    xs = list(map(binary64, inputs))
    require(all(two(-100) <= x <= two(100) for x in xs), 'product-closure magnitude')
    models = ('strict', 'stored_pc64', 'resident_pc64', 'stored_pc53', 'resident_pc53', 'exact')
    def mul(a, b, model):
        return a*b if model == 'exact' else operation(a, b, 'mul', model)
    table = {}
    for mask in range(1, 1 << len(xs)):
        indices = [i for i in range(len(xs)) if mask & (1 << i)]
        values = set()
        for order in itertools.permutations(indices):
            if len(order) == 1:
                values.add(xs[order[0]])
            elif len(order) == 2:
                values.update(mul(xs[order[0]], xs[order[1]], m) for m in models)
            else:
                a, b, c = (xs[i] for i in order)
                for first, last in itertools.product(models, repeat=2):
                    values.add(mul(mul(a, b, first), c, last))
                    values.add(mul(a, mul(b, c, first), last))
        table[str(mask)] = [[str(x.numerator), str(x.denominator)] for x in sorted(values)]
    target = binary64(raw['expected'])
    products = [Q(*map(int, p)) for p in table[str((1 << len(xs))-1)]]
    domains = root_publication_domains(products, len(xs), target, staged)
    result = dict(raw=raw, products=table, root_span=domains['binary64']['span'],
                  reachable=any(d['reachable'] for d in domains.values()))
    if staged:
        result['root_roundings'] = domains
    return result


def check_product_closure(input_path, directory, output, staged=False):
    directory = Path(directory)
    frozen = {str(p):digest(p) for p in (Path(input_path), Path(__file__),
        Path(__file__).with_name('exact.py'), Path(__file__).with_name('check.py'),
        *sorted(directory.iterdir())) if p.is_file()}
    counts = Counter()
    proofs, scope = iter(records(directory/'closures.jsonl')), iter(records(directory/'scope.jsonl'))
    for index, raw in enumerate(records(input_path)):
        inputs, reason = coordinates(raw)
        if reason is None and not 2 <= len(inputs) <= 3:
            reason = 'arity_outside_2_to_3'
        if reason is None and any(not two(-100) <= binary64(x) <= two(100) for x in inputs):
            reason = 'outside_modest_magnitudes'
        if reason is None and binary64(raw['expected']) <= 0:
            reason = 'nonpositive_expected'
        require(next(scope, None) == dict(trace_index=index, reason=reason), 'changed product scope')
        if reason:
            counts[reason] += 1
            continue
        independently = product_result(raw, staged)
        require(next(proofs, None) == dict(trace_index=index, **independently), 'incomplete product closure, answer or source differs')
        counts['scoped'] += 1
        counts['unreachable'] += not independently['reachable']
        counts['root_product_values'] += len(independently['products'][str((1 << len(inputs))-1)])
        if staged:
            for name, result in independently['root_roundings'].items():
                counts[name+'_unreachable'] += not result['reachable']
    require(counts['scoped'] and next(scope, None) is None and next(proofs, None) is None, 'empty or extra product obligations')
    summary = json.loads((directory/'summary.json').read_text())
    require(summary['complete'] is True and summary['counts'] == counts, 'incorrect product closure totals')
    require(summary.get('root_staging', False) == staged, 'incorrect product root-stage grammar')
    require(frozen == {p:digest(p) for p in frozen}, 'product closure inputs/checker changed')
    receipt = dict(complete=True, counts=counts, input_sha256=frozen,
        limits='Independent explicit permutation/two-tree reconstruction of complete product sets and exact midpoint-power root preimages. '
               'Two/three modest positive operands; six choices per multiply (five rounding/store modes or exact). '
               'Reference bindings are re-read from original workbook/capture bytes. No native power leaf or complete Excel algorithm is inferred.')
    if staged:
        receipt.update(root_staging=True,
            limits='Independent complete product-tree sets and exact inverse publication/nth-power boundaries for direct binary64, extended64-then-binary64, or extended53-then-binary64 real roots. '
                   'Every intermediate root is assumed correctly rounded in its declared format; native pow/log approximations and other programs remain open. '
                   'The direct-root family is included, not added. Original static reference lineage is revalidated.')
    with open(output, 'x') as out:
        json.dump(receipt, out, indent=2)
        out.write('\n')
    print(json.dumps(receipt))


def verified_service_batch(cache,receipt,function='GEOMEAN'):
    """Read the one existing scalar-worker lineage, with exact input echoes."""
    import gzip
    path=Path(cache)/receipt['file']
    require(Path(receipt['file']).name==receipt['file']and digest(path)==receipt['sha256'],'service reply bytes changed') # ubs:ignore — public evidence digest and filename, not secrets.
    with gzip.open(path,'rt')as stream:raw=json.load(stream)
    reply=raw['reply'];request=raw['request']
    require(request['function']==function and reply['ok']and reply['version']==1,'wrong/failed service reply')
    expected=[[{k:v for k,v in a.items()if k!='mode'}for a in row]for row in request['rows']]
    require(reply['inputs']==expected and len(reply['rows'])==len(request['rows']),'service operands/answers incomplete')
    require(reply['automation_security_before']==reply['automation_security_after'],'service worker changed security')
    require(all(reply[k]==receipt[k]for k in ('request_sha256','excel_build','machine','timestamp')),'service provenance changed')
    return raw,path


def verified_study_answers(doc,cache):
    """Collect one study without discarding duplicate original observations."""
    answers={};bindings={};ordered=[]
    for receipt in doc['receipts']:
        raw,path=verified_service_batch(cache,receipt,receipt['function'])
        bindings[str(path)]=digest(path)
        request,reply=raw['request'],raw['reply']
        require(len(request['rows'])==receipt['calls'],'study receipt count changed')
        for row,a in zip(request['rows'],reply['rows']):
            key=request['function'],json.dumps(row,sort_keys=True)
            require(key not in answers or answers[key]==a,'contradicting original duplicate readings')
            answers[key]=a;ordered.append((request['function'],row,a))
    return answers,bindings,ordered


def check_stored_products(study,cache,before,after,originals,output):
    """Independent product/store proof and exhaustive changed-result audit.

    Exact rounding is the checker's integer/rational arithmetic. The root
    controls are independent Excel observations. Neither a native root nor
    a mathematical root supplies their expected values.
    """
    import gzip
    from exact import BINARY64 as B,EXTENDED64 as E
    from x87 import bits64
    study,cache,before,after,originals,output=map(Path,(study,cache,before,after,originals,output))
    frozen={str(p):digest(p)for p in (study,before,after,originals,Path(__file__),Path(__file__).with_name('exact.py'))}
    doc=json.loads(study.read_text())
    require(doc['complete']and len(doc['records'])==61,'incomplete product study')
    require([r['scale_exponent']for r in doc['records']]==list(range(-30,31)),'changed predeclared scale grid')
    require(len(doc['receipts'])==1,'unexpected product-control batch inventory')
    receipt=doc['receipts'][0];raw,path=verified_service_batch(cache,receipt)
    frozen[str(path)]=digest(path)
    require(raw['request']['function']=='GEOMEAN'and len(raw['request']['rows'])==183,'wrong product questions')
    reply=raw['reply']
    require(reply['ok']and reply['inputs']==raw['request']['rows']and len(reply['rows'])==183,'product operands/answers incomplete')
    require(reply['automation_security_before']==reply['automation_security_after'],'product worker changed security')
    require(all(reply[k]==receipt[k]for k in ('request_sha256','excel_build','machine','timestamp')),'product provenance changed')
    key=lambda args:tuple(a.removeprefix('0x')for a in args)
    controls={};counts=Counter()
    for i,r in enumerate(doc['records']):
        require(r['args'][:2]==['401f9ff8b1ae9622']*2,'changed independently known square input')
        require(binary64(r['args'][2])==two(r['scale_exponent']),'scale is not the declared exact power of two')
        strict=stored=Q(1)
        for arg in r['args']:
            x=binary64(arg);strict=B.round(strict*x);stored=B.round(E.round(stored*x))
        require(f'{bits64(strict):016x}'==r['strict']and f'{bits64(stored):016x}'==r['stored_pc64'],'incorrect independently rounded products')
        expected_questions=[[dict(type='number',bits=b)for b in r['args']]]
        expected_questions += [[dict(type='number',bits=r[k]),dict(type='number',bits='3ff0000000000000'),dict(type='number',bits='3ff0000000000000')]for k in ('strict','stored_pc64')]
        require(raw['request']['rows'][3*i:3*i+3]==expected_questions,'product/root question operands changed')
        target,old_root,new_root=reply['rows'][3*i:3*i+3]
        require([target,old_root,new_root]==[r['original'],r['strict_root_control'],r['stored_root_control']],'reported root controls differ from worker')
        require(all(a['type']=='number'for a in (target,old_root,new_root)),'unread product control')
        outcome=('root publications indistinguishable'if old_root==new_root else'strict product only'if target==old_root
                 else'storedPC64 product only'if target==new_root else'neither product root control')
        require(outcome==r['outcome'],'incorrect product classification');counts[outcome]+=1
        controls[key(r['args'])]=dict(original=target,stored_root=new_root,
            stored_root_args=key([r['stored_pc64'],'3ff0000000000000','3ff0000000000000']),record=i)
    require(dict(counts)==doc['counts']and counts['storedPC64 product only']and not counts['strict product only']
            and not counts['neither product root control'],'product evidence does not isolate stored-PC64')
    old_rows,new_rows=list(records(before)),list(records(after))
    original_rows=json.loads(originals.read_text())
    require(len(old_rows)==len(new_rows)==len(original_rows),'missing before/after live predictions')
    source_answers={}
    for name in sorted({s for r in original_rows for s in r['sources']}):
        source=Path(name).resolve()
        require(source.is_relative_to(cache.resolve()),'live source escapes the retained service cache')
        frozen[str(source)]=digest(source)
        with gzip.open(source,'rt')as stream:reading=json.load(stream)
        request,response=reading['request'],reading['reply']
        require(request['function']=='GEOMEAN'and response['ok']and response['inputs']==request['rows']
                and len(response['rows'])==len(request['rows']),'live source operands/answers changed')
        require(response['automation_security_before']==response['automation_security_after'],'live source security differs')
        values={}
        for args,answer in zip(request['rows'],response['rows']):
            require(all(a['type']=='number'for a in args),'non-numeric live replay source')
            k=key([a['bits']for a in args])
            require(k not in values or values[k]==answer,'conflicting repeated source answer')
            values[k]=answer
        source_answers[name]=values
    predictions={key(r['args']):r['calc_func_prediction']for r in new_rows}
    require(len(predictions)==len(new_rows),'duplicate replay arguments')
    totals=Counter();losses=[];changed=[]
    for old,new,original in zip(old_rows,new_rows,original_rows):
        request=original['request'];args=key(request['args']);answer=original['answer']
        require(old['id']==new['id']==request['id']and old['args']==new['args']==request['args'],'changed original replay identity')
        require(original['sources']and all(source_answers[s].get(args)==answer for s in original['sources']),
                'live expectation differs from its exact raw Excel question/answer')
        require(answer['type']=='number','unsupported live result type')
        expected='0x'+answer['bits'];was=old['calc_func_prediction']==expected;now=new['calc_func_prediction']==expected
        for name,v in dict(rows=1,before_matched=was,after_matched=now,gains=now and not was,losses=was and not now,
                           changed_outputs=old['calc_func_prediction']!=new['calc_func_prediction']).items():totals[name]+=v
        if old['calc_func_prediction']!=new['calc_func_prediction']:
            require(args in controls,'changed output outside the isolated product study')
            control=controls[args]
            require(answer==control['original']==control['stored_root'],'original/root-control answer changed')
            require(new['calc_func_prediction']==predictions[control['stored_root_args']],
                    'new original call differs from its stored-product root control')
            entry=dict(id=request['id'],args=request['args'],expected=expected,before=old['calc_func_prediction'],after=new['calc_func_prediction'],
                       product_study_record=control['record'])
            changed.append(entry)
            if was and not now:
                entry=dict(entry,classification='unresolved root helper at the independently measured stored product',
                           measured_stored_product_root=control['stored_root'])
                losses.append(entry)
    require(totals['losses']==len(losses),'unclassified live loss')
    require(frozen=={p:digest(Path(p))for p in frozen},'product evidence/checker changed during verification')
    result=dict(complete=True,product_control_counts=counts,live_counts=totals,losses=losses,changed=changed,input_sha256=frozen,
        limits='Sixty-one normal scalar square-tie controls isolate only the product loop; no private root is identified. Every changed result is tied to an independent actual root control. Other recorded bank/corpus scopes need their own comparison. No coefficient, tolerance, expected answer or selector is fitted.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(product_control_counts=counts,live_counts=totals,classified_losses=len(losses))))
    return result


def check_intermediate_stores(study,cache,output):
    """Independently distinguish per-step stores from one final product store."""
    from exact import BINARY64 as B,EXTENDED64 as E
    from x87 import bits64
    study,cache,output=map(Path,(study,cache,output))
    frozen={str(p):digest(p)for p in (study,Path(__file__),Path(__file__).with_name('exact.py'))}
    doc=json.loads(study.read_text());models=['strict','stored_each','resident_final']
    require(doc['complete'] and doc['models']==models and len(doc['records'])==512,'incomplete intermediate-store study')
    require(doc['consecutive_inputs_per_center']==256 and len(doc['receipts'])==1,'changed predeclared window size')
    raw,path=verified_service_batch(cache,doc['receipts'][0]);frozen[str(path)]=digest(path)
    require(len(raw['request']['rows'])==2048,'incomplete intermediate-store questions')
    centers=[bits64(B.round(Q(13,10))),bits64(B.round(Q(37,10)))]
    counts=Counter();product_differences=Counter()
    for i,r in enumerate(doc['records']):
        require(r['args'][:2]==['401f9ff8b1ae9622']*2,'changed square witness')
        require(r['args'][2]==f'{centers[i//256]-128+i%256:016x}','changed consecutive input window')
        calculated={}
        for model in models:
            p=Q(1)
            for arg in r['args']:
                x=binary64(arg)
                p=(E.round(p*x)if model=='resident_final'else B.round(E.round(p*x))if model=='stored_each'else B.round(p*x))
            calculated[model]=f'{bits64(B.round(p)):016x}'
        require(r['products']==calculated,'incorrect independently computed intermediate products')
        product_differences['stored_vs_resident']+=calculated['stored_each']!=calculated['resident_final']
        product_differences['stored_vs_strict']+=calculated['stored_each']!=calculated['strict']
        questions=[[dict(type='number',bits=a)for a in r['args']]]
        questions += [[dict(type='number',bits=calculated[m]),dict(type='number',bits='3ff0000000000000'),dict(type='number',bits='3ff0000000000000')]for m in models]
        require(raw['request']['rows'][4*i:4*i+4]==questions,'changed candidate root operands')
        actual,*roots=raw['reply']['rows'][4*i:4*i+4]
        require(actual==r['original']and dict(zip(models,roots))==r['root_controls'],'reported controls differ from raw worker')
        require(all(x['type']=='number'for x in [actual,*roots]),'non-numeric intermediate-store control')
        compatible=[m for m,root in zip(models,roots)if root==actual]
        require(compatible==r['compatible'],'incorrect model compatibility')
        require('stored_each'in compatible,'measured contradiction to stored-each program')
        label='only '+compatible[0] if len(compatible)==1 else 'ambiguous '+','.join(compatible)
        counts[label]+=1
    require(dict(counts)==doc['counts']and counts['only stored_each']>0,'no distinguishing per-step-store observation')
    require(frozen=={p:digest(Path(p))for p in frozen},'intermediate-store evidence changed')
    result=dict(complete=True,counts=counts,product_differences=product_differences,input_sha256=frozen,
        limits='Three declared complete product programs, three positive normal factors, one known square tie, fixed512 third-factor inputs. The per-step stored program is consistent;55 sole-model publications distinguish it from the alternatives. Controls assume binary64 publication before the unknown root; an unpublished resident extended-root producer remains outside this certificate.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(counts=counts,product_differences=product_differences)))
    return result


def check_stage_replay(study,cache,before,after,output):
    """Count every stage question and classify all changed original results."""
    study,cache,before,after,output=map(Path,(study,cache,before,after,output))
    proof=check_intermediate_stores(study,cache,output.with_suffix('.products.json'))
    doc=json.loads(study.read_text());raw,path=verified_service_batch(cache,doc['receipts'][0])
    old_rows,new_rows=list(records(before)),list(records(after))
    require(len(old_rows)==len(new_rows)==2048,'missing full stage replay')
    frozen=dict(proof['input_sha256'],**{str(p):digest(p)for p in (before,after)})
    totals=Counter();original_counts=Counter();losses=[];changes=[]
    for i,(old,new,question,answer)in enumerate(zip(old_rows,new_rows,raw['request']['rows'],raw['reply']['rows'])):
        args=['0x'+a['bits']for a in question]
        require(old['id']==new['id']==f'stage-full/{i}'and old['args']==new['args']==args,'stage replay identity changed')
        expected='0x'+answer['bits'];a,b=old['calc_func_prediction'],new['calc_func_prediction']
        was,now=a==expected,b==expected
        delta=dict(rows=1,before_matched=was,after_matched=now,gains=now and not was,losses=was and not now,changed_outputs=a!=b)
        for k,v in delta.items():totals[k]+=v
        if i%4==0:
            for k,v in delta.items():original_counts[k]+=v
        else:require(a==b,'unit-factor root control changed with the product step')
        if a!=b:
            require(i%4==0,'changed non-original stage control')
            root_index=i+2
            require(answer==raw['reply']['rows'][root_index],'original stage answer differs from measured stored-product root')
            require(b==new_rows[root_index]['calc_func_prediction'],'new stage result differs from the independently invoked stored-product root')
            entry=dict(index=i//4,args=args,expected=expected,before=a,after=b,
                       stored_product=doc['records'][i//4]['products']['stored_each'])
            changes.append(entry)
            if was and not now:
                losses.append(dict(entry,classification='unresolved root helper at independently measured stored product',
                                   measured_stored_product_root=raw['reply']['rows'][root_index]))
    require(len(losses)==totals['losses'],'unclassified stage loss')
    require(frozen=={p:digest(Path(p))for p in frozen},'stage replay inputs changed')
    result=dict(complete=True,all_reading_counts=totals,original_counts=original_counts,losses=losses,changes=changes,input_sha256=frozen,
        limits='All512 original triples and1536 independently measured root controls counted. Controls include duplicate arguments and are not added to other cohorts as unique calls. Every changed original result agrees with the corresponding current stored-product root invocation; the root remains unresolved.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_reading_counts=totals,original_counts=original_counts,classified_losses=len(losses))))
    return result


def check_discrete_seed_stages(study,cache,output):
    """Check original discrete observations and one-node inverse constraints."""
    from exact import Span,BINARY64 as B,pull_rounding
    study,cache,output=map(Path,(study,cache,output))
    doc=json.loads(study.read_text());require(doc['complete'],'incomplete discrete stage study')
    frozen={str(p):digest(p)for p in (study,Path(__file__),Path(__file__).with_name('exact.py'))}
    answers,bindings,ordered=verified_study_answers(doc,cache);frozen.update(bindings)
    used=set();counts=Counter();aliases=Counter()
    def number(bits):return dict(type='number',bits=bits)
    def scalar(x):return number(f'{bits64(Q(x)):016x}')
    def logical(x):return dict(type='boolean',value=x,mode='cell')
    def observed(name,row):
        key=name,json.dumps(row,sort_keys=True);used.add(key);return answers[key]
    from x87 import bits64
    if doc['format']=='dnacalc-hypergeom-identical-factor-return-stage/1':
        groups={};seen=set();one=scalar(1)
        for r in doc['records']:
            d,n,k=r['draws'],r['total'],r['population_successes'];key=d,n,k
            require(key not in seen and k in (0,n)and 0<d<n and n&(n-1)==0,'wrong hypergeometric stage domain')
            seen.add(key)
            s=0 if k==0 else d;args=[scalar(v)for v in (s,d,k,n)]
            pmf=observed('HYPGEOM.DIST',args+[logical(False)])
            cdf=observed('HYPGEOM.DIST',args+[logical(True)])
            legacy=observed('HYPGEOMDIST',args)
            b=observed('BINOM.DIST',[scalar(d),scalar(n),scalar(Q(d,n)),logical(False)])
            require((pmf,cdf,legacy,b)==(r['pmf'],r['cdf'],r['legacy'],r['public_binom']),'discrete record differs from original readings')
            groups.setdefault((d,n),{})[k]=pmf
            aliases['equal'if pmf==legacy else'different']+=1
            counts['pmf_mathematical_one'if pmf==one else'pmf_differs_from_one']+=1
            counts['cdf_mathematical_one'if cdf==one else'cdf_differs_from_one']+=1
            require(pmf['type']==b['type']=='number','unhandled nonnumeric return control')
            stored=binary64(b['bits']);target=binary64(pmf['bits'])
            stored_cell=B.pullback(Span(stored,stored));cell=pull_rounding(Span(target,target),'stored_pc64')
            domains={}
            for orientation in ('resident_over_stored','stored_over_resident'):
                allowed=cell.linear(stored)if orientation=='resident_over_stored'else cell.reciprocal(stored)
                intersection=allowed.meet(stored_cell)
                domains[orientation]=intersection.json()if intersection is not None else None
                counts[orientation+'_possible'if intersection is not None else orientation+'_impossible']+=1
            require(domains==r['conditional_resident_mass_intervals'],'wrong conditional mass-return inverse')
        for (d,n),group in groups.items():
            require(set(group)=={0,n},'missing zero/all population pair')
            counts['zero_all_same'if group[0]==group[n]else'zero_all_differ']+=1
    elif doc['format']=='dnacalc-poisson-resident-seed-multiply-stage/1':
        inputs={r['input']for r in doc['records']};seen=set()
        for x in inputs:
            rows=[[scalar(k),number(x),logical(False)]for k in (0,1,2)]
            ps=[observed('POISSON.DIST',row)for row in rows]
            for row,a in zip(rows,ps):aliases['equal'if observed('POISSON',row)==a else'different']+=1
            exponential=observed('EXP',[scalar(-binary64(x))])
            counts['zero_equals_EXP'if ps[0]==exponential else'zero_differs_from_EXP']+=1
        for r in doc['records']:
            x,stage=r['input'],r['stage'];require((x,stage)not in seen,'duplicate Poisson stage');seen.add((x,stage))
            index={'zero_to_one':0,'one_to_two':1}[stage]
            a=observed('POISSON.DIST',[scalar(index),number(x),logical(False)])
            b=observed('POISSON.DIST',[scalar(index+1),number(x),logical(False)])
            require((a,b)==(r['published_seed'],r['published_output']),'Poisson stage target differs from original reading')
            require(a['type']==b['type']=='number','unhandled nonnumeric seed stage')
            seed,target=binary64(a['bits']),binary64(b['bits']);multiplier=binary64(x)/(index+1)
            source=B.pullback(Span(seed,seed));target_cell=pull_rounding(Span(target,target),'stored_pc64')
            domain=source.meet(target_cell.linear(1/multiplier))
            require((domain.json()if domain is not None else None)==r['resident_seed_interval'],'Poisson seed inverse differs')
            counts[stage+'_possible'if domain is not None else stage+'_impossible']+=1
        require(seen=={(x,k)for x in inputs for k in ('zero_to_one','one_to_two')},'missing Poisson seed pair')
        require([r['mean']for r in doc['mode_records']]==[*range(1,65),128,256,500,1000],'missing original mode observation')
        for r in doc['mode_records']:
            mean=r['mean'];require(r['previous_count']==mean-1,'mode count changed')
            a,b=[observed('POISSON.DIST',[scalar(k),scalar(mean),logical(False)])for k in (mean-1,mean)]
            require((a,b)==(r['previous'],r['current']),'mode answer differs from original reading')
            counts['exact_mode_step_equal'if a==b else'exact_mode_step_different']+=1
    else:raise ValueError('unknown discrete stage format')
    require(used==set(answers),'unaccounted original discrete observation')
    require(counts==doc['counts']and aliases==doc['alias_counts'],'discrete stage totals differ')
    require(frozen=={p:digest(p)for p in frozen},'discrete stage source or evidence changed')
    result=dict(complete=True,original_calls=len(ordered),unique_original_calls=len(answers),counts=counts,
                alias_counts=aliases,inputs_sha256=frozen,
                limits='Independent rational inversion of one stored-PC64 node using the complete original scalar-worker readings. Every question and alias observation accounted for. Conditional seed intervals do not identify a private seed or authorize a correction.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items()if k!='inputs_sha256'}));return result


def quotient_summary_coverage(models,groups):
    require([(m['quotient'],m['published'])for m in models]==
        [(q,p)for q in ('divide','stored_reciprocal')for p in (False,True)],'quotient model order changed')
    for model in models:
        q,p=model['quotient'],model['published']
        require([g['product_power']for g in model['groups']]==sorted(k for a,b,k in groups if(a,b)==(q,p)),
                'missing, repeated or unordered quotient summary')


def check_root_public_stages(study,original,original_cache,public_cache,output):
    """Verify the measured public leaves and exact stored argument independently."""
    from exact import BINARY64 as B, EXTENDED64 as E
    study,original,original_cache,public_cache,output=map(Path,(study,original,original_cache,public_cache,output))
    doc=json.loads(study.read_text());prior=json.loads(original.read_text())
    require(doc['complete']and prior['complete']and doc['original_sha256']==digest(original),'public-stage original changed')
    require(doc['original_receipts']==prior['receipts'],'public-stage receipt scope changed')
    roots,root_files,_=verified_study_answers(prior,original_cache)
    public,public_files,ordered=verified_study_answers(doc,public_cache)
    frozen={str(p):digest(p)for p in (study,original,Path(__file__))};frozen.update(root_files);frozen.update(public_files)
    number=lambda b:dict(type='number',bits=b)
    key=lambda name,b:(name,json.dumps([number(b)],sort_keys=True))
    expected={}
    for row in prior['records']:
        n=row['count'];b=int(row['input'],16);e=((b>>52)&2047)-1023
        q=(abs(e)//n)*(-1 if e<0 else 1);shift=e-q*n+1023
        require(0<shift<2047,'public-stage input not normal')
        x=f'{(b&((1<<52)-1))|(shift<<52):016x}'
        args=[number(x)]+[number('3ff0000000000000')]*(n-1)
        expected[n,x]=roots['GEOMEAN',json.dumps(args,sort_keys=True)]
    require([(r['count'],r['input'])for r in doc['records']]==sorted(expected),'missing, repeated or unordered reduced root')
    used=set();counts=Counter()
    for row in doc['records']:
        n,x=row['count'],row['input'];lkey=key('LN',x);log=public[lkey];used.add(lkey)
        require(log['type']=='number','public-stage log not numeric')
        reciprocal=B.round(E.round(Q(1,n)));argument=B.round(E.round(binary64(log['bits'])*reciprocal))
        # bits64 accepts an exact binary64 value; it does not round again.
        from x87 import bits64
        rbits=f'{bits64(reciprocal):016x}';abits=f'{bits64(argument):016x}'
        ekey=key('EXP',abits);exp=public[ekey];used.add(ekey);match=exp==expected[n,x]
        require((row['excel'],row['public_LN'],row['reciprocal'],row['stored_argument'],row['public_EXP'],row['matched'])==
            (expected[n,x],log,rbits,abits,exp,match),'public-stage leaf or stored argument changed')
        counts['matched'if match else'mismatched']+=1
    inputs=sorted({x for _,x in expected});singletons=Counter()
    require([r['input']for r in doc['singleton_controls']]==inputs,'singleton scope changed')
    for row in doc['singleton_controls']:
        x=row['input'];skey=key('GEOMEAN',x);used.add(skey);answer=public[skey]
        require(row['excel']==answer,'singleton answer changed')
        singletons['identity'if answer==number(x)else'different']+=1
    require(used==set(public),'unused public primitive reading')
    require(counts==doc['counts']and singletons==doc['singleton_counts'],'public-stage totals changed')
    require(doc['public_log_calls']==len(inputs)and doc['public_exp_calls']==len({r['stored_argument']for r in doc['records']}),'public-stage call count changed')
    require(frozen=={p:digest(p)for p in frozen},'public-stage source or reply changed')
    result=dict(complete=True,reduced_roots=len(expected),public_calls=len(ordered),counts=counts,singleton_counts=singletons,inputs_sha256=frozen,
        limits='Complete declared reduced-root subset and independently measured LN/EXP/singleton calls. Exact PC64 then binary64 assignments. This identifies the sampled composition only; private implementation internals and fallback remain open.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items()if k!='inputs_sha256'}));return result


def check_root_log_stage(study,original,original_cache,control_cache,library,output):
    """Verify one quotient's inverse using independent rational/ISA operations."""
    from x87 import encode,decode,load64
    from exact import Span,BINARY64 as B
    study,original,original_cache,control_cache,output=map(Path,(study,original,original_cache,control_cache,output))
    doc=json.loads(study.read_text());prior=json.loads(original.read_text())
    require(doc['complete']and prior['complete']and doc['original_sha256']==digest(original),'quotient source changed')
    require(doc['original_receipts']==prior['receipts'],'quotient original scope changed')
    old,old_files,_=verified_study_answers(prior,original_cache)
    public,public_files,public_order=verified_study_answers(doc,control_cache)
    frozen={str(p):digest(p)for p in (study,original,Path(__file__),Path(library))};frozen.update(old_files);frozen.update(public_files)
    number=lambda b:dict(type='number',bits=b)
    source=[r for r in prior['records']if r['role']=='fractional_log2_root'and r['product_power']>0]
    expected={(q,p,r['count'],r['product_power'])for q in ('divide','stored_reciprocal')for p in (False,True)for r in source}
    seen=set();groups={};used=Counter();boundaries=0
    exact=Exact(Instructions(library));cw=0x133f;one=encode(Q(1));zero=encode(Q(0))
    def op(kind,a=zero,b=zero,stored=False):return exact.operation(kind,a,b,cw,stored)
    for row in doc['records']:
        q,p,n,k=row['quotient'],row['published'],row['count'],row['product_power'];key=q,p,n,k
        require(key in expected and key not in seen,'quotient missing, extra or duplicate record');seen.add(key)
        args=[number(row['input'])]+[number('3ff0000000000000')]*(n-1)
        require(binary64(row['input'])==two(k)and row['excel']==old['GEOMEAN',json.dumps(args,sort_keys=True)],'quotient original target changed')
        require(0<k<n,'anchor normalization is not q=0')
        operand=Q(n)if q=='divide'else B.round(Q(1,n))
        quotient_node=dict(op='div'if q=='divide'else'mul',a=0,b=2)
        if p:quotient_node['site']=0
        nodes=[dict(op='input',a=0),dict(op='const',bits='3ff0000000000000'),
            dict(op='const',bits=f'{published_bits(encode(operand)):016x}'),quotient_node,
            dict(op='fldl2e',a=0,b=0),dict(op='mul',a=3,b=4),
            dict(op='frndint',a=5,b=0),dict(op='sub',a=5,b=6),
            dict(op='f2xm1',a=7,b=0),dict(op='add',a=1,b=8),dict(op='fscale',a=9,b=6)]
        require(row['graph']==dict(nodes=nodes,outputs=[10,10],sites=int(p)),'quotient graph metadata changed')
        def evaluate(x):
            t=op('div'if q=='divide'else'mul',x,encode(operand),p)
            t=op('mul',t,op('fldl2e'));whole=op('frndint',t)
            small=op('f2xm1',op('sub',t,whole))
            return f'{published_bits(op("fscale",op("add",one,small),whole)):016x}'
        cell=row['cell'];require(cell is not None,'empty quotient inverse needs a separate certificate')
        lo,hi=(bytes.fromhex(cell[z])for z in ('low','high'))
        def neighbor(raw,delta):
            sig=int.from_bytes(raw[:8],'little')+delta
            require(1<<63<=sig<1<<64,'quotient inverse crosses binade')
            return sig.to_bytes(8,'little')+raw[8:]
        got=[evaluate(x)for x in (neighbor(lo,-1),lo,hi,neighbor(hi,1))];boundaries+=4
        require(got[1]==got[2]==row['excel']['bits']and int(got[0],16)<int(got[1],16)<int(got[3],16),'quotient inverse boundaries changed')
        require((cell['predecessor_output'],cell['successor_output'])==(got[0],got[3]),'quotient neighbor metadata changed')
        domain=Span(decode(lo),decode(hi));g=q,p,k
        groups[g]=domain if g not in groups else groups[g].meet(domain)if groups[g]is not None else None
        used[g]+=1
    require(seen==expected,'incomplete quotient original subset')
    keys={('LN',json.dumps([number(f'{published_bits(encode(two(k))):016x}')],sort_keys=True))for _,_,k in groups}
    require(keys==set(public),'missing or unused public LN observation')
    counts={};models_seen=set()
    quotient_summary_coverage(doc['models'],groups)
    for model in doc['models']:
        q,p=model['quotient'],model['published'];require((q,p)not in models_seen,'duplicate quotient model');models_seen.add((q,p))
        tally=Counter()
        for item in model['groups']:
            k=item['product_power'];domain=groups[q,p,k];answer=public['LN',json.dumps([number(f'{published_bits(encode(two(k))):016x}')],sort_keys=True)]
            possible=domain is not None and answer['type']=='number'and domain.contains(binary64(answer['bits']))
            require(item['observations']==used[q,p,k]and item['common_private_log']==(domain.json()if domain else None)and item['public_LN']==answer and item['public_LN_possible']==possible,'quotient group summary changed')
            tally['groups']+=1;tally['impossible']+=domain is None;tally['public_LN_excluded']+=not possible
        counts[f'{q}/{p}']=tally
    require(models_seen=={(q,p)for q in ('divide','stored_reciprocal')for p in (False,True)},'missing quotient model')
    require(frozen=={p:digest(p)for p in frozen},'quotient sources changed')
    result=dict(complete=True,original_subset=len(source),public_calls=len(public_order),boundary_checks=boundaries,counts=counts,inputs_sha256=frozen,
        limits='Independent exact quotient operations and separate ISA leaves. The public EXP primitive is conditional; no private log or exp implementation is identified.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(counts));return result


def check_root_normalization(study,original,original_cache,control_cache,output):
    """Integer-only exponent reconstruction using independent Excel root controls."""
    study,original,original_cache,control_cache,output=map(Path,(study,original,original_cache,control_cache,output))
    doc=json.loads(study.read_text());prior=json.loads(original.read_text())
    require(doc['complete']and prior['complete']and doc['original_study_sha256']==digest(original),
            'normalization source binding changed')
    require(doc['original_receipts']==prior['receipts'],'normalization original receipt scope changed')
    originals,old_files,ordered=verified_study_answers(prior,original_cache)
    controls,control_files,control_order=verified_study_answers(doc,control_cache)
    frozen={str(p):digest(p)for p in (study,original,Path(__file__))};frozen.update(old_files);frozen.update(control_files)
    require(len(doc['records'])==len(ordered),'missing original normalization result')
    names=['frexp_truncation','adjusted_frexp_truncation','adjusted_frexp_floor']
    require(doc['convention_order']==names,'normalization convention order changed')
    counts={name:Counter()for name in names};used=set();first_contradictions={}
    for index,(r,(function,args,a))in enumerate(zip(doc['records'],ordered)):
        require(function=='GEOMEAN'and r['index']==index and r['count']==len(args)and
                args[0]==dict(type='number',bits=r['input'])and r['excel']==a and
                all(v==dict(type='number',bits='3ff0000000000000')for v in args[1:]),
                'original normalization input or answer changed')
        n=len(args);b=int(r['input'],16);e=((b>>52)&2047)-1023
        require(0<((b>>52)&2047)<2047 and b>>63==0,'nonpositive or nonnormal original product')
        require([v['convention']for v in r['conventions']]==names,'missing normalization convention')
        for item in r['conventions']:
            name=item['convention'];exponent=e+(name=='frexp_truncation')
            quotient=(exponent//n if name.endswith('_floor')else
                      (abs(exponent)//n)*(-1 if exponent<0 else 1))
            shifted=e-quotient*n+1023
            require(0<shifted<2047,'normalized product leaves normal range')
            reduced=f'{(b&((1<<52)-1))|(shifted<<52):016x}'
            require((item['exponent'],item['quotient'],item['reduced'])==(exponent,quotient,reduced),
                    'normalization exponent/quotient differs from exact integer arithmetic')
            key='GEOMEAN',json.dumps([dict(type='number',bits=reduced)]+[dict(type='number',bits='3ff0000000000000')]*(n-1),sort_keys=True)
            used.add(key);measured=controls[key]
            require(measured['type']=='number'and item['measured_normalized_root']==measured,'normalized root is not the original control')
            bits=int(measured['bits'],16);scale_exponent=((bits>>52)&2047)+quotient
            require(0<scale_exponent<2047 and bits>>63==0,'reconstructed root leaves positive normal range')
            prediction=dict(type='number',bits=f'{(bits&((1<<52)-1))|(scale_exponent<<52):016x}')
            require(item['prediction']==prediction,'final exact exponent restoration differs')
            match=prediction==a;counts[name]['matched'if match else'mismatched']+=1
            if not match and name not in first_contradictions:first_contradictions[name]=dict(index=index,input=r['input'],count=n,excel=a,prediction=prediction)
    require(used==set(controls),'unaccounted normalized-root control')
    require(counts==doc['counts'],'normalization totals differ')
    require(frozen=={p:digest(p)for p in frozen},'normalization source or reading changed')
    result=dict(complete=True,original_calls=len(ordered),original_unique_calls=len(originals),
        normalized_calls=len(control_order),normalized_unique_calls=len(controls),counts=counts,
        first_contradictions=first_contradictions,inputs_sha256=frozen,
        limits='Exact integer exponent arithmetic and complete original GEOMEAN controls. Tests only the declared normalization/reconstruction stage; it neither identifies nor replaces the private root. A surviving convention is scoped to these normal products, arities and input route.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in result.items()if k!='inputs_sha256'}));return result


def primitive_stage_grammar(kind,constant53=False,stored_reduction=False,spill=None):
    """Independent fixed ISA metadata, separate from the producer's Builder."""
    nodes=[dict(op='input',a=0),dict(op='const',bits='3ff0000000000000')]
    def node(op,a=0,b=0,store=False):
        item=dict(op=op,a=a,b=b)
        if store:item['site']=0
        nodes.append(item);return len(nodes)-1
    def const(bits):nodes.append(dict(op='const',bits=bits));return len(nodes)-1
    if kind in ('LN','LOG2'):
        c=1 if kind=='LOG2'else const('3fe62e42fefa39ef')if constant53 else node('fldln2')
        final=node('fyl2x',c,0)
    else:
        require(kind in ('EXP','EXP2'),'unknown primitive grammar')
        if kind=='EXP':
            c=const('3ff71547652b82fe')if constant53 else node('fldl2e')
            argument=node('mul',0,c,stored_reduction)
        else:argument=0
        integer=node('frndint',argument)
        fraction=node('sub',argument,integer,spill=='fraction')
        small=node('f2xm1',fraction,store=spill=='f2xm1_return')
        mantissa=node('add',1,small)
        final=node('fscale',mantissa,integer)
    return dict(nodes=nodes,outputs=[final,final],sites=int(stored_reduction or spill is not None))


def check_stage_graph_metadata(graphs):
    expected=dict(LN_ROM80=primitive_stage_grammar('LN'),LN_constant53=primitive_stage_grammar('LN',True),
        EXP_ROM80_resident_reduction=primitive_stage_grammar('EXP'),
        EXP_constant53_resident_reduction=primitive_stage_grammar('EXP',True),
        EXP_ROM80_stored_reduction=primitive_stage_grammar('EXP',stored_reduction=True),
        EXP2=primitive_stage_grammar('EXP2'),LOG2=primitive_stage_grammar('LOG2'))
    require(graphs==expected,'primitive graph metadata differs from the independently declared ISA grammar')


def stage_observation_coverage(doc, answers, programs):
    """Every inverse obligation comes from the complete original answer set."""
    require([p['program']for p in doc['primitive_checks']]==list(programs),
            'missing, duplicate or reordered primitive programmes')
    original_root_inputs={}
    one=dict(type='number',bits='3ff0000000000000')
    for name,encoded in answers:
        if name!='GEOMEAN':continue
        row=json.loads(encoded)
        require(len(row)in (2,4,8)and all(a==one for a in row[1:]),
                'root observation is outside the exact unit-factor premise')
        original_root_inputs.setdefault(row[0]['bits'],set()).add(len(row))
    require(original_root_inputs and all(ns=={2,4,8}for ns in original_root_inputs.values()),
            'missing original root arity')
    expected_pairs={(x,k)for x in original_root_inputs for k in ('EXP','EXP2')}
    listed_pairs=[(r['input'],r['primitive'])for r in doc['root_intervals']]
    require(len(listed_pairs)==len(expected_pairs)and set(listed_pairs)==expected_pairs,
            'missing, extra or duplicated original root inverse input')


def check_resident_stage_inverses(document, cache, library, output, constants=None, refinement=None, refinement_cache=None, spills=None):
    """Independent rational arithmetic and one-instruction inverse-cell audit.

    Reconstructs the declared primitives directly, without importing the
    stage generator, its graph interpreter or the application. Local
    monotonicity of each declared primitive remains an explicit premise.
    """
    import gzip
    from x87 import encode,decode,load64,published_bits
    from exact import BINARY64 as B, EXTENDED64 as E, Span
    document,cache,output=map(Path,(document,cache,output))
    doc=json.loads(document.read_text())
    require(doc['complete']and doc['format']=='dnacalc-resident-primitive-stage-inverses/1','partial or unknown inverse study')
    check_stage_graph_metadata(doc['graphs'])
    frozen={str(p):digest(p)for p in (document,Path(__file__),Path(__file__).with_name('x87.py'),Path(__file__).with_name('exact.py'),Path(library))}
    answers,bindings,ordered=verified_study_answers(doc,cache)
    frozen.update(bindings);calls=len(ordered)
    def observed(name,row):
        return answers[(name,json.dumps(row,sort_keys=True))]
    number=lambda b:dict(type='number',bits=b)
    exact=Exact(Instructions(library));cw=0x133f
    zero,one=encode(Q(0)),encode(Q(1))
    def op(kind,a=zero,b=zero,stored=False):return exact.operation(kind,a,b,cw,stored)
    def primitive(kind,x,constant53=False,stored=False,spill=None):
        if kind in ('LN','LOG2'):
            c=one if kind=='LOG2'else load64(0x3fe62e42fefa39ef)if constant53 else op('fldln2')
            return op('fyl2x',c,x)
        if kind=='EXP':
            c=load64(0x3ff71547652b82fe)if constant53 else op('fldl2e')
            t=op('mul',x,c,stored)
        else:require(kind=='EXP2','unknown inverse primitive');t=x
        whole=op('frndint',t)
        fraction=op('sub',t,whole,spill=='fraction')
        small=op('f2xm1',fraction,stored=spill=='f2xm1_return')
        mantissa=op('add',one,small)
        return op('fscale',mantissa,whole)
    programs={
        'EXP_ROM80_resident_reduction':('EXP',False,False),
        'EXP_constant53_resident_reduction':('EXP',True,False),
        'EXP_ROM80_stored_reduction':('EXP',False,True),
        'LN_ROM80':('LN',False,False),'LN_constant53':('LN',True,False)}
    stage_observation_coverage(doc,answers,programs)
    primitive_counts={p:Counter()for p in programs};root_counts={k:Counter()for k in ('EXP','EXP2')}
    for (name,encoded),a in answers.items():
        if name not in ('LN','EXP'):continue
        args=json.loads(encoded);require(len(args)==1,'primitive arity changed')
        for p,(kind,constant53,stored)in programs.items():
            if name!=kind:continue
            got=f'{published_bits(primitive(kind,load64(int(args[0]["bits"],16)),constant53,stored)):016x}'
            primitive_counts[p]['matched'if a==number(got)else'mismatched']+=1
    for record in doc['primitive_checks']:
        require(primitive_counts[record['program']]==record['counts'],'wrong primitive stage totals')
    seen=set();boundary_checks=0;contradictions=[]
    for r in doc['root_intervals']:
        key=r['input'],r['primitive'];require(key not in seen,'duplicate root inverse');seen.add(key)
        lows,highs=[],[]
        require([c['count']for c in r['cells']]==[2,4,8],'missing exact-scale root arity')
        for cell in r['cells']:
            n=cell['count'];a=cell['excel'];span=cell['cell']
            require(observed('GEOMEAN',[number(r['input'])]+[number('3ff0000000000000')]*(n-1))==a,
                    'inverse target is not the original Excel observation')
            require(span is not None,'empty-cell certificate needs both neighboring input witnesses')
            lo,hi=(bytes.fromhex(span[k])for k in ('low','high'))
            require(decode(lo)<=decode(hi),'reversed inverse cell')
            def neighbor(raw,delta):
                sig=int.from_bytes(raw[:8],'little')+delta
                require(1<<63<=sig<1<<64,'inverse neighbor crosses binade')
                return sig.to_bytes(8,'little')+raw[8:]
            got=[f'{published_bits(primitive(r["primitive"],v)):016x}'for v in (neighbor(lo,-1),lo,hi,neighbor(hi,1))]
            boundary_checks+=4
            require(got[1]==got[2]==a['bits']and int(got[0],16)<int(a['bits'],16)<int(got[3],16),
                    'inverse boundary does not contain exactly the requested publication cell')
            require((span['predecessor_output'],span['successor_output'])==(got[0],got[3]),'changed neighboring publications')
            lows.append(decode(lo)*n);highs.append(decode(hi)*n)
        lo,hi=max(lows),min(highs);feasible=lo<=hi
        require(feasible==r['shared_log_possible'],'wrong shared-log feasibility')
        tally=root_counts[r['primitive']]
        tally['shared_log_possible'if feasible else'shared_log_impossible']+=1
        if feasible:
            require(r['common_log_interval']==[str(lo),str(hi)],'changed shared-log interval')
            kind='LN'if r['primitive']=='EXP'else'LOG2'
            resident=decode(primitive(kind,load64(int(r['input'],16))))
            tally['registered_resident_log_in_intersection'if lo<=resident<=hi else'registered_resident_log_excluded']+=1
            tally['binary64_log_in_intersection'if lo<=B.round(resident)<=hi else'binary64_log_excluded']+=1
        else:contradictions.append(dict(input=r['input'],primitive=r['primitive'],gap=str(lo-hi)))
    require(root_counts==doc['root_stage_counts'],'inverse stage counts do not reconcile')
    require(all({(x,k)for k in ('EXP','EXP2')}<=seen for x,_ in seen),'missing paired inverse primitive')
    constant_families=[];constant_boundary_checks=0
    if constants is not None:
        constants=Path(constants);frozen[str(constants)]=digest(constants)
        cdoc=json.loads(constants.read_text());require(cdoc['complete'],'partial constant study')
        require([f['function']for f in cdoc['families']]==['LN','EXP'],
                'missing, duplicate or reordered constant family')
        common={};constant_seen=set()
        for r in cdoc['records']:
            name,x=r['function'],r['input'];key=name,x
            require(name in ('LN','EXP')and key not in constant_seen,'duplicate or unknown constant input')
            constant_seen.add(key)
            require(observed(name,[number(x)])==r['excel'],'constant target differs from original reading')
            span=r['cell'];lo,hi=(bytes.fromhex(span[k])for k in ('low','high'))
            raw_x=load64(int(x,16));q=decode(raw_x)
            def evaluate(v):
                raw=op('fyl2x',v,raw_x)if name=='LN'else primitive('EXP2',v)
                return f'{published_bits(raw):016x}'
            got=[evaluate(v)for v in (neighbor(lo,-1),lo,hi,neighbor(hi,1))]
            constant_boundary_checks+=4
            target=r['excel']['bits']
            require(got[1]==got[2]==target and int(got[0],16)<int(target,16)<int(got[3],16),
                    'constant-stage inverse boundary invalid')
            require((span['predecessor_output'],span['successor_output'])==(got[0],got[3]),'constant-stage neighbor changed')
            if name=='LN':domain=Span(decode(lo),decode(hi))
            else:
                # Union of all RN64 cells from the first through last allowed
                # product. Division by the exact positive input is monotone.
                lower=(decode(neighbor(lo,-1))+decode(lo))/(2*q)
                upper=(decode(hi)+decode(neighbor(hi,1)))/(2*q)
                domain=Span(lower,upper,int.from_bytes(lo[:8],'little')%2==0,int.from_bytes(hi[:8],'little')%2==0)
            require(domain.json()==r['coefficient_interval'],'wrong multiplication/coefficient inverse')
            common[name]=domain if name not in common else common[name].meet(domain)if common[name]is not None else None
        expected_inputs={(name,json.loads(row)[0]['bits'])for name,row in answers if name in ('LN','EXP')}
        require(constant_seen==expected_inputs,'missing constant input or extra constant observation')
        for family in cdoc['families']:
            name=family['function'];span=common[name]
            require((span.json()if span is not None else None)==family['common_coefficient_interval'],'constant intersection changed')
            rom=op('fldln2'if name=='LN'else'fldl2e')
            require(family['published_rom']==rom.hex()and family['rom_in_intersection']==(span is not None and span.contains(decode(rom))),
                    'ROM comparison differs from independently measured ISA constant')
            candidates=[]
            if span is not None:
                begin=E.ceiling_coordinate(span.low);end=E.floor_coordinate(span.high)
                if not span.contains(E.value(begin)):begin+=1
                if not span.contains(E.value(end)):end-=1
                total=max(0,end-begin+1)
                candidates=[encode(E.value(k)).hex()for k in range(begin,min(end+1,begin+16))]
            else:total=0
            require(total==family['pc64_candidate_count']and candidates==family['first_pc64_candidates'],'wrong constant lattice cardinality')
            constant_families.append(family)
    refined_constants=[];refinement_calls=0
    if refinement is not None:
        require(constants is not None and refinement_cache is not None,'refinement needs the original constant proof and raw answers')
        refinement=Path(refinement);refinement_cache=Path(refinement_cache)
        selection=Path(str(refinement)+'.selection.json')
        for p in (refinement,selection):frozen[str(p)]=digest(p)
        rdoc=json.loads(refinement.read_text());selected=json.loads(selection.read_text())
        # ubs:ignore — public immutable evidence identities, no secret or authentication boundary.
        require(rdoc['complete']and rdoc['prior_sha256']==digest(constants)==selected['prior_sha256']and
                rdoc['selection_sha256']==digest(selection),'refinement prior/selection binding changed')
        prior=next(f for f in constant_families if f['function']=='LN')
        candidates=prior['first_pc64_candidates']
        require(candidates==rdoc['original_candidates']and len(candidates)==prior['pc64_candidate_count'],
                'refinement silently drops a prior coefficient')
        require(selected['coefficient_endpoints']==[candidates[0],candidates[-1]],'refinement endpoints changed')
        _,refinement_bindings,refinement_order=verified_study_answers(rdoc,refinement_cache)
        frozen.update(refinement_bindings)
        require(all(name=='LN'for name,_,_ in refinement_order),'wrong refinement function')
        original=[(row,a)for _,row,a in refinement_order];refinement_calls=len(original)
        require(len(original)==len(rdoc['records'])==len(selected['selected']),'unaccounted refinement observations')
        prior_span=prior['common_coefficient_interval']
        intersection=Span(Q(*map(int,prior_span['lower'])),Q(*map(int,prior_span['upper'])),
                          prior_span['lower_closed'],prior_span['upper_closed'])
        for r,s,(row,a)in zip(rdoc['records'],selected['selected'],original):
            require({k:r[k]for k in s}==s and row==[number(r['input'])]and r['excel']==a,
                    'refinement record differs from selected question or original answer')
            x=load64(int(r['input'],16))
            predictions=[f'{published_bits(op("fyl2x",bytes.fromhex(c),x)):016x}'for c in (candidates[0],candidates[-1])]
            require(predictions[0]!=predictions[1]and predictions==[r['lower_prediction'],r['upper_prediction']],
                    'refinement question did not distinguish the original endpoint coefficients')
            span=r['cell'];lo,hi=(bytes.fromhex(span[k])for k in ('low','high'))
            got=[f'{published_bits(op("fyl2x",v,x)):016x}'for v in (neighbor(lo,-1),lo,hi,neighbor(hi,1))]
            require(got[1]==got[2]==a['bits']and int(got[0],16)<int(a['bits'],16)<int(got[3],16),
                    'refined coefficient output cell is invalid')
            domain=Span(decode(lo),decode(hi))
            intersection=intersection.meet(domain)if intersection is not None else None
        remaining=[c for c in candidates if intersection is not None and intersection.contains(decode(bytes.fromhex(c)))]
        require(remaining==rdoc['remaining_pc64_candidates']and
                (intersection.json()if intersection is not None else None)==rdoc['common_coefficient_interval'],
                'refined coefficient result changed')
        refined_constants=remaining
    spill_counts={};spill_boundary_checks=0
    if spills is not None:
        spills=Path(spills);frozen[str(spills)]=digest(spills)
        sdoc=json.loads(spills.read_text());require(sdoc['complete'],'partial primitive-store study')
        models=[(k,s)for k in ('EXP','EXP2')for s in ('fraction','f2xm1_return')]
        require([(m['primitive'],m['store'])for m in sdoc['models']]==models,'missing, duplicate or reordered primitive-store model')
        for m in sdoc['models']:
            require(m['graph']==primitive_stage_grammar(m['primitive'],spill=m['store']),
                    'primitive-store graph metadata differs from independent grammar')
        inputs={r['input']for r in doc['root_intervals']}
        expected={(x,k,s)for x in inputs for k,s in models}
        keys=[(r['input'],r['primitive'],r['store'])for r in sdoc['records']]
        require(len(keys)==len(expected)and set(keys)==expected,'missing original primitive-store inverse')
        for r in sdoc['records']:
            key=r['primitive'],r['store'];tally=spill_counts.setdefault(key,Counter())
            require([c['count']for c in r['cells']]==[2,4,8],'missing original primitive-store arity')
            lows,highs=[],[]
            for c in r['cells']:
                n=c['count'];a=c['excel'];span=c['cell']
                require(observed('GEOMEAN',[number(r['input'])]+[number('3ff0000000000000')]*(n-1))==a,
                        'primitive-store target differs from original observation')
                require(span is not None,'empty-cell store proof needs neighboring witnesses')
                lo,hi=(bytes.fromhex(span[k])for k in ('low','high'))
                got=[f'{published_bits(primitive(r["primitive"],v,spill=r["store"])):016x}'for v in
                     (neighbor(lo,-1),lo,hi,neighbor(hi,1))]
                spill_boundary_checks+=4
                require(got[1]==got[2]==a['bits']and int(got[0],16)<int(a['bits'],16)<int(got[3],16),
                        'primitive-store inverse boundary changed')
                require((span['predecessor_output'],span['successor_output'])==(got[0],got[3]),'primitive-store adjacent output changed')
                lows.append(decode(lo)*n);highs.append(decode(hi)*n)
            lower,upper=max(lows),min(highs);possible=lower<=upper
            require(possible==r['shared_log_possible']and
                    (not possible or r['common_log_interval']==[str(lower),str(upper)]),'primitive-store shared-log interval changed')
            tally['shared_log_possible'if possible else'shared_log_impossible']+=1
            if possible:
                center=decode(primitive('LN'if r['primitive']=='EXP'else'LOG2',load64(int(r['input'],16))))
                tally['registered_log_possible'if lower<=center<=upper else'registered_log_excluded']+=1
        for m in sdoc['models']:
            require(spill_counts[m['primitive'],m['store']]==m['counts'],'primitive-store counts differ')
    require(frozen=={p:digest(p)for p in frozen},'inverse source or observation changed')
    result=dict(complete=True,original_calls=calls,primitive_counts=primitive_counts,
                inverse_boundary_checks=boundary_checks,root_stage_counts=root_counts,
                constant_boundary_checks=constant_boundary_checks,constant_families=constant_families,
                refinement_calls=refinement_calls,refined_log_constants=refined_constants,
                spill_boundary_checks=spill_boundary_checks,
                spill_counts=[dict(primitive=k,store=s,counts=c)for (k,s),c in spill_counts.items()],
                contradictions=contradictions,input_sha256=frozen,
                limits='Independent rational PC64 arithmetic and native single-instruction leaves. Exact output-cell endpoints and their immediate neighbors verified; local monotonicity is conditional. The intersections bound an unknown shared logarithm, never fit one. No private algorithm or global table/polynomial identification follows.')
    with output.open('x')as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(original_calls=calls,inverse_boundary_checks=boundary_checks,root_stage_counts=root_counts)))
    return result


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='root-public-stages':
        check_root_public_stages(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='root-log-stage':
        check_root_log_stage(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='root-normalization':
        check_root_normalization(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='discrete-stages':
        check_discrete_seed_stages(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='resident-stages':
        check_resident_stage_inverses(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='stage-replay':
        check_stage_replay(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='intermediate-stores':
        check_intermediate_stores(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='stored-products':
        check_stored_products(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='product-closure-staged':
        check_product_closure(*sys.argv[2:], staged=True);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='product-closure':
        check_product_closure(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='exact-root':check_exact_roots(*sys.argv[2:])
    else:check(*sys.argv[1:])
