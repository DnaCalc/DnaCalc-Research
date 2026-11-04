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


def verified_service_batch(cache,receipt):
    """Read the one existing scalar-worker lineage, with exact input echoes."""
    import gzip
    path=Path(cache)/receipt['file']
    require(Path(receipt['file']).name==receipt['file']and digest(path)==receipt['sha256'],'service reply bytes changed') # ubs:ignore — public evidence digest and filename, not secrets.
    with gzip.open(path,'rt')as stream:raw=json.load(stream)
    reply=raw['reply'];request=raw['request']
    require(request['function']=='GEOMEAN'and reply['ok']and reply['version']==1,'wrong/failed service reply')
    require(reply['inputs']==request['rows']and len(reply['rows'])==len(request['rows']),'service operands/answers incomplete')
    require(reply['automation_security_before']==reply['automation_security_after'],'service worker changed security')
    require(all(reply[k]==receipt[k]for k in ('request_sha256','excel_build','machine','timestamp')),'service provenance changed')
    return raw,path


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


if __name__=='__main__':
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
