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
    require(digest(package)==binding['package_sha256']and digest(capture)==binding['capture_sha256'],'reference source changed')
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
    require(summary['complete'] and summary['programs']==len(ps) and summary['survivors']==len(survivors)
            and summary['observations']==len(rows) and summary['scope']==scope and summary['executions']==executions,
            'wrong root summary')
    require(all(digest(p)==h for p,h in frozen.items()), 'root proof/source changed')
    result=dict(complete=True,programs=len(ps),observations=len(rows),survivors=len(survivors),executions=executions,
                input_sha256=frozen,limits='Exact inverse midpoint-power verification of the specified rounded-product/correct-root programs. No native transcendental instruction, fitted answer or Excel internal precision is assumed.')
    with open(output,'x')as out:json.dump(result,out,indent=2);out.write('\n')
    print(json.dumps(result))


def product_result(raw):
    """Explicit ordered trees; independent of the producer's subset closure."""
    from exact import operation, BINARY64 as B, Span
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
    span = B.pullback(Span(target, target))
    interval = Span(span.low**len(xs), span.high**len(xs), span.left, span.right)
    products = (Q(*map(int, p)) for p in table[str((1 << len(xs))-1)])
    return dict(raw=raw, products=table, root_span=span.json(), reachable=any(interval.contains(p) for p in products))


def check_product_closure(input_path, directory, output):
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
        independently = product_result(raw)
        require(next(proofs, None) == dict(trace_index=index, **independently), 'incomplete product closure, answer or source differs')
        counts['scoped'] += 1
        counts['unreachable'] += not independently['reachable']
        counts['root_product_values'] += len(independently['products'][str((1 << len(inputs))-1)])
    require(counts['scoped'] and next(scope, None) is None and next(proofs, None) is None, 'empty or extra product obligations')
    summary = json.loads((directory/'summary.json').read_text())
    require(summary['complete'] is True and summary['counts'] == counts, 'incorrect product closure totals')
    require(frozen == {p:digest(p) for p in frozen}, 'product closure inputs/checker changed')
    receipt = dict(complete=True, counts=counts, input_sha256=frozen,
        limits='Independent explicit permutation/two-tree reconstruction of complete product sets and exact midpoint-power root preimages. '
               'Two/three modest positive operands; six choices per multiply (five rounding/store modes or exact). '
               'Reference bindings are re-read from original workbook/capture bytes. No native power leaf or complete Excel algorithm is inferred.')
    with open(output, 'x') as out:
        json.dump(receipt, out, indent=2)
        out.write('\n')
    print(json.dumps(receipt))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='product-closure':
        check_product_closure(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='exact-root':check_exact_roots(*sys.argv[2:])
    else:check(*sys.argv[1:])
