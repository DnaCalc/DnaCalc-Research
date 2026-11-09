"""Paired public gamma calls, with one prefactor and one series stage.

Original MIT research code. Reuses the scalar recorder, ISA interpreter and
exact inverse arithmetic. Candidate labels are fixed before reading answers;
counts report contradictions and never select coefficients or a correction.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys

try:
    from .fingerprint import (CachedService, bits, value, load_service, Native,
                             primitive_stage_graph, positive_extended_cell, Builder, sha, write)
except ImportError:
    from fingerprint import (CachedService, bits, value, load_service, Native,
                             primitive_stage_graph, positive_extended_cell, Builder, sha, write)
from exact import BINARY64, EXTENDED64, Span, binary64, inverse
from x87 import load64, decode

STAGE_BINDINGS = {str(p.resolve()):sha(p) for p in (
    Path(__file__), Path(__file__).with_name('x87_programs.py'),
    Path(__file__).with_name('x87_native.c'),
    *(Path(__file__).parents[1]/'checker'/name for name in
      ('geomean.py', 'exact.py', 'check.py', 'x87.py')))}


def seed_phase_tracks(checkout, output, cache, retained, offline=False, radius=.001, count=1024):
    """Continue the observed GAMMA(1.3) lines, without a candidate table.

    The retained 20k window and neighboring windows supplied these frequencies
    and slopes. They are observations, not a search over named approximations.
    New calls form a fixed geometric ladder; cached replies retain readbacks.
    """
    from fingerprint import PhaseTrack, window_residual
    if not math.isfinite(radius) or not 0 < radius <= .125 or count < 16:
        raise ValueError('phase study needs 16 or more samples and a radius in (0,.125]')
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    retained = Path(retained).resolve()
    rows = [json.loads(line) for line in retained.read_text().splitlines()][:4096]
    xs = [value(r['a'][0]) for r in rows]
    if len(rows) != 4096 or any(int(rows[i]['a'][0],16) != int(rows[0]['a'][0],16)+i for i in range(len(rows))):
        raise ValueError('retained consecutive window changed')
    x0 = xs[0]
    specifications = ((.03475,-.29), (.44946,-.49), (.14258,.50),
                      (.34679,.21), (.47787,-.04), (.20379,.275), (.04187,-1.10))
    tracks = [PhaseTrack('L'+str(i+1),x0,f,s,degree=12) for i,(f,s) in enumerate(specifications)]
    residual = window_residual([r['r'] for r in rows])
    for t in tracks:
        t.observe(xs,residual)
    centers=[]
    step=6000*math.ulp(x0)
    distance=step
    while distance<radius:
        centers.extend((x0+distance,x0-distance))
        distance=distance*1.25+step
    windows=[]
    # Independent questions may be packed; phase fits still consume windows
    # sequentially. The queue's existing 20k-call limit remains in force.
    for offset in range(0,len(centers),16):
        batch=centers[offset:offset+16]
        inputs=[f'{int(bits(c),16)-count//2+k:016x}' for c in batch for k in range(count)]
        reply=service.ask('GAMMA',[[client.number(x)] for x in inputs])
        for i,c in enumerate(batch):
            answers=reply['rows'][i*count:(i+1)*count]
            ws=[value(x) for x in inputs[i*count:(i+1)*count]]
            residual=window_residual(answers)
            readings=[t.observe(ws,residual) for t in tracks]
            windows.append(dict(center=bits(c),readings=readings))
            print(json.dumps(dict(center=c,readings=readings)),flush=True)
        result=dict(format='dnacalc-gamma-seed-phase-reader/1',complete=offset+16>=len(centers),
            retained_window=dict(path=str(retained),sha256=sha(retained),rows=len(rows),
                                 limits='Retained scalar rows omit worker metadata; bound to their original file, not promoted to a new capture.'),
            specifications=[dict(frequency=f,slope=s) for f,s in specifications],
            radius=radius,count=count,windows=windows,tracks=[t.record() for t in tracks],
            receipts=service.receipts,source_bindings={str(p):sha(p) for p in (Path(__file__),Path(__file__).with_name('fingerprint.py'))},
            limits='Diagnostic phase polynomials only. No coefficient or store is installed; aliases, absolute stage scales and additive constants remain open. No bank match counts select this study.')
        write(output,result)
    return result


def seed_phase_validation(study, output, cache, retained_dir=None):
    """Check a read phase and one Horner link against retained coarse windows.

    The coarse windows selected initial frequencies/slopes, but never entered
    the phase polynomial fit. They are validation of that fit, not blind data.
    Only L5's constant and slope are calibrated; its higher derivatives come
    from the fixed L5 = 2*(x-1)*L1 + constant hypothesis.
    """
    import gzip
    from decimal import Decimal as D, localcontext
    from fingerprint import PhaseTrack, window_residual
    study=Path(study).resolve()
    doc=json.loads(study.read_text())
    if not doc['complete']:
        raise ValueError('incomplete phase ladder')
    retained=Path(retained_dir)/'q21.jsonl' if retained_dir else Path(doc['retained_window']['path'])
    if sha(retained)!=doc['retained_window']['sha256']:
        raise ValueError('retained seed window changed')
    rows=[json.loads(line) for line in retained.read_text().splitlines()][:4096]
    xs=[value(r['a'][0]) for r in rows]
    noise=window_residual([r['r'] for r in rows])
    def calibrated(index):
        spec=doc['specifications'][index]
        track=PhaseTrack('L'+str(index+1),xs[0],spec['frequency'],spec['slope'],degree=12)
        track.observe(xs,noise)
        return track
    root=calibrated(0)
    root.coefficients=[D(c) for c in doc['tracks'][0]['coefficients_cycles']]
    root.active_degree=doc['tracks'][0]['active_degree']
    root.observations=[(value(r['x']),D(r['cycles']),r['amplitude']) for r in doc['tracks'][0]['observations']]
    root.misses=doc['tracks'][0]['misses']
    linked=calibrated(4)
    linked.observations=[]
    linked.coefficients=[D(0)]*len(root.coefficients)
    linked.coefficients[1]=D(doc['tracks'][4]['coefficients_cycles'][1])
    with localcontext() as ctx:
        ctx.prec=70
        for k in range(2,len(linked.coefficients)):
            linked.coefficients[k]=2*((D(xs[0])-1)*root.coefficients[k]+root.coefficients[k-1])
    # Fit a constant and slope only, preserving every derived higher term.
    linked.degree=1
    linked.active_degree=1
    linked.observe(xs,noise)
    training=[]
    raw_bindings={str(retained):sha(retained)}
    for receipt in doc['receipts']:
        path=Path(cache)/receipt['file']
        if sha(path)!=receipt['sha256']:
            raise ValueError('phase reply changed')
        raw_bindings[str(path.resolve())]=sha(path)
        with gzip.open(path,'rt') as stream:
            original=json.load(stream)
        if original['request']['function']!='GAMMA':
            raise ValueError('wrong phase function')
        for i in range(0,len(original['reply']['rows']),doc['count']):
            xx=[value(r[0]['bits']) for r in original['request']['rows'][i:i+doc['count']]]
            residual=window_residual(original['reply']['rows'][i:i+doc['count']])
            training.append(linked.observe(xx,residual))
    validations=[]
    directory=Path(retained_dir) if retained_dir else retained.parent
    for tag in 'abcdefg':
        path=directory/('q24'+tag+'.jsonl')
        raw_bindings[str(path.resolve())]=sha(path)
        rows=[json.loads(line) for line in path.read_text().splitlines()][:4096]
        xx=[value(r['a'][0]) for r in rows]
        residual=window_residual([r['r'] for r in rows])
        results=[]
        for label,track in (('read_L1',root),('linked_L5',linked)):
            _,_,amplitude,phase=track.demodulate(xx,residual)
            innovation=(phase-track.reference_phase+.5)%1-.5
            results.append(dict(model=label,amplitude=amplitude,innovation_cycles=innovation))
        validations.append(dict(path=str(path.resolve()),first=bits(xx[0]),count=len(xx),results=results))
    summary={label:max(abs(r['innovation_cycles']) for w in validations for r in w['results'] if r['model']==label)
             for label in ('read_L1','linked_L5')}
    result=dict(format='dnacalc-gamma-seed-phase-validation/1',complete=True,
        study=dict(path=str(study),sha256=sha(study)),raw_bindings=raw_bindings,
        root_phase=root.record(),linked_phase=linked.record(),training=training,validation=validations,
        summary=dict(max_abs_validation_phase_cycles=summary,
            linked_training_windows=len(training),
            linked_large_innovations=sum(abs(r.get('innovation_cycles',0))>.2 for r in training)),
        hypothesis='L5 = 2*(x-1)*L1 + constant, modulo the arbitrary phase offset and integer/sign slope aliases; higher terms are inherited, constant and slope alone are calibrated.',
        source_bindings={str(p):sha(p) for p in (Path(__file__),Path(__file__).with_name('fingerprint.py'))},
        limits='A supported smooth phase relation, not a literal private-store proof. The other five tracks, absolute stage scale, polynomial offsets and remaining table coefficients are unidentified. Coarse windows informed initial slopes and are not blind observations. No application value changes or bank gains are claimed.')
    write(output,result)
    print(json.dumps(result['summary']),flush=True)
    return result


def seed_phase_absolute(study, output):
    """Expose the remaining affine/scale ambiguity of the R32 Horner lead.

    This is algebra on already measured phase derivatives, not a recovered
    numerical table. Each phase is measured modulo one on a binary64 grid.
    """
    from decimal import Decimal as D, localcontext
    doc = json.loads(Path(study).read_text())
    if doc['format'] != 'dnacalc-gamma-seed-phase-validation/1':
        raise ValueError('validated phase relation required')
    root, linked = doc['root_phase'], doc['linked_phase']
    with localcontext() as ctx:
        ctx.prec = 70
        x0 = D(value(root['origin']))
        dx = D(math.ulp(float(x0)))
        p, q = ([D(c) for c in track['coefficients_cycles']] for track in (root, linked))
        # This is an identity imposed by the linked R32 model; its
        # inherited second derivative is not another observation.
        alias = (q[2]/2-(x0-1)*p[2]-p[1])*dx
        root_value = q[1]/2-(x0-1)*p[1]
        choices = []
        for integer in range(-8, 9):
            normalized = root_value+D(integer)/(2*dx)
            if D(2)**52 <= abs(normalized) < D(2)**53:
                choices.append(dict(linked_slope_integer_alias=integer,
                    root_value_in_own_ulps=str(normalized),
                    root_value_in_binades=str(normalized*D(2)**-52)))
        higher = {str(j):str(p[j]*dx) for j in range(2,root['active_degree']+1)}
    result = dict(format='dnacalc-gamma-seed-phase-absolute/1', complete=True,
        study=dict(path=str(Path(study).resolve()),sha256=sha(study)),
        relation='Q = 2*(x-1)*P + constant in normalized stage cycles',
        inherited_second_derivative_identity=str(alias),
        conditional_premise='Chosen root slope branch (integer alias zero) and observed spectral signs; linked higher derivatives were inherited, not independently measured.',
        conditional_root_value_choices=choices,
        higher_taylor_coefficients_in_input_grid_units=higher,
        source_bindings={str(Path(__file__).resolve()):sha(__file__)},
        limits='Conditional algebra using the R32 smooth-phase model, not rigorous hidden-stage enclosures. '
               'The linked slope has integer aliases, both spectral signs need interpretation, the root phase '
               'has an arbitrary constant, and a common power-of-two scale is unobserved. The remaining '
               'phase tracks do not yet supply validated adjacent Horner links. These observations do not '
               'identify absolute binary64 coefficients; no table or application result changes.')
    write(output,result)
    print(json.dumps(result,indent=2),flush=True)
    return result


def low_gamma_recurrence(checkout, output, cache, library, offline=False, count=32):
    """Separate simple quotient loops using independently observed seed calls."""
    from fractions import Fraction as Q
    from exact import operation
    native=Native(library)
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    seeds=sorted({f'{int(bits(center),16)-count//2+i:016x}'
                  for center in (1.125,1.25,1.3,1.375,1.5,1.625,1.7,1.75,1.875)
                  for i in range(count)})
    cases=[]
    for seed in seeds:
        y=value(seed)-1.
        for shift in (0,-1,-2):
            x=y+shift
            if not -2<x<1 or x==0 or Q(x)!=Q(y)+shift:
                continue
            factors=[y-i for i in range(1-shift)]
            if any(Q(f)!=Q(y)-i for i,f in enumerate(factors)):
                raise ValueError('recurrence factor rounded')
            cases.append(dict(x=bits(x),seed=seed,factors=[bits(f)for f in factors]))
    labels=[(order,kind) for order in ('top_down','bottom_up')
            for kind in ('resident_divisions','stored_divisions','resident_product','stored_product')]
    questions=dict(format='dnacalc-gamma-low-recurrence-questions/1',complete=True,
        seeds=seeds,cases=cases,program_order=labels,
        selection='Fixed dyadic/ordinary centers and adjacent seed doubles. Exact rational operands '
                  'must preserve the same fractional seed. No output selects the inputs or program order.',
        source_bindings={**STAGE_BINDINGS,str(Path(library).resolve()):sha(library)})
    write(str(output)+'.questions.json',questions)
    answers=service.ask('GAMMA',[[client.number(s)]for s in seeds])['rows']
    if any(a['type']!='number'for a in answers):
        raise ValueError('positive seed must be numeric')
    anchors={s:a['bits']for s,a in zip(seeds,answers)}
    targets=service.ask('GAMMA',[[client.number(c['x'])]for c in cases])['rows']
    counts=Counter();records=[]
    for c,target in zip(cases,targets):
        if target['type']!='number':
            raise ValueError('non-pole low recurrence must be numeric')
        predictions={}
        for order,kind in labels:
            factors=c['factors'] if order=='top_down'else list(reversed(c['factors']))
            operand=load64(int(anchors[c['seed']],16))
            rational=binary64(anchors[c['seed']])
            if kind.endswith('divisions'):
                for f in factors:
                    operand=native.leaf('div',operand,load64(int(f,16)),0x133f)
                    rational=EXTENDED64.round(rational/binary64(f))
                    if kind.startswith('stored'):
                        raw=native.publish(operand)
                        operand=load64(raw)
                        rational=BINARY64.round(rational)
            else:
                product=load64(int(bits(1.),16));rproduct=Q(1)
                for f in factors:
                    product=native.leaf('mul',product,load64(int(f,16)),0x133f)
                    rproduct=EXTENDED64.round(rproduct*binary64(f))
                    if kind.startswith('stored'):
                        product=load64(native.publish(product))
                        rproduct=BINARY64.round(rproduct)
                operand=native.leaf('div',operand,product,0x133f)
                rational=EXTENDED64.round(rational/rproduct)
            raw=f'{native.publish(operand):016x}'
            if raw!=bits(float(BINARY64.round(rational))):
                raise ValueError('native recurrence differs from exact rational checker')
            label=order+'_'+kind
            predictions[label]=raw
            counts[label+'_matches']+=raw==target['bits']
        counts['rows']+=1
        counts['all_models_refuted']+=target['bits']not in predictions.values()
        records.append(dict(c,seed_answer=anchors[c['seed']],answer=target,predictions=predictions))
    result=dict(format='dnacalc-gamma-low-recurrence/1',complete=True,counts=counts,
        records=records,receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        source_bindings=questions['source_bindings'],
        limits='Paired public seeds, exact factor identities and independently checked PC64 arithmetic. '
               'A match identifies only a sampled recurrence under the public-seed premise; it does not '
               'identify the seed table. Every contradiction and original return is retained.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def low_gamma_shared_seed(study, output):
    """Invert each low recurrence before assuming the seed is stored."""
    from collections import defaultdict
    from fractions import Fraction as Q
    from exact import operation
    doc=json.loads(Path(study).read_text())
    groups=defaultdict(list)
    for r in doc['records']:
        groups[r['seed']].append(r)
    programs=[(order,kind)for order in ('top_down','bottom_up')
              for kind in ('resident_divisions','stored_divisions','resident_product','stored_product')]
    counts=Counter();records=[]
    for seed,rows in sorted(groups.items()):
        anchor=binary64(rows[0]['seed_answer'])
        anchor_cell=BINARY64.pullback(Span(anchor,anchor))
        models=[]
        for order,kind in programs:
            common=anchor_cell
            constraints=[]
            for r in rows:
                factors=r['factors']if order=='top_down'else list(reversed(r['factors']))
                v=binary64(r['answer']['bits'])
                span=BINARY64.pullback(Span(v,v))
                if kind.endswith('divisions'):
                    for f in reversed(factors):
                        span=inverse(span,binary64(f),'div',True,
                                     'stored_pc64'if kind.startswith('stored')else'resident_pc64')
                else:
                    product=Q(1)
                    for f in factors:
                        product=operation(product,binary64(f),'mul',
                                          'stored_pc64'if kind.startswith('stored')else'resident_pc64')
                    span=inverse(span,product,'div',True,'resident_pc64')
                constraints.append(dict(x=r['x'],seed_preimage=span.json()if span else None))
                common=common.meet(span)if common else None
            coordinates=EXTENDED64.coordinates(common)
            label=order+'_'+kind
            counts[label+'_feasible_groups']+=coordinates is not None
            models.append(dict(program=label,shared_real_interval=common.json()if common else None,
                               resident_count=coordinates[1]-coordinates[0]+1 if coordinates else 0,
                               constraints=constraints))
        counts['seed_groups']+=1
        counts['all_models_refuted_groups']+=all(m['resident_count']==0 for m in models)
        records.append(dict(seed=seed,seed_answer=rows[0]['seed_answer'],models=models))
    result=dict(format='dnacalc-gamma-low-shared-seed/1',complete=True,counts=counts,records=records,
        study=dict(path=str(Path(study).resolve()),sha256=sha(study)),
        source_bindings={**STAGE_BINDINGS},
        limits='Necessary exact inverse constraints for a common arbitrary resident PC64 seed in the '
               'observed public seed rounding cell. An empty intersection refutes the declared recurrence '
               'even without assuming binary64 seed storage. A nonempty interval is possibility, not '
               'identification of a seed producer or its final operation. No output changes.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def gammaln_low_composition(checkout, output, cache, library, offline=False, count=64):
    """Observe the below-.7 recurrence with public upper-band and LN calls."""
    from exact import operation
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    centers=(2.**-100,2.**-53,2.**-52,2.**-20,.125,.25,.375,.5,.625,
             math.nextafter(.7,0.))
    inputs=sorted({f'{int(bits(c),16)-count//2+i:016x}'for c in centers for i in range(count)
                   if 0<value(f'{int(bits(c),16)-count//2+i:016x}')<.7})
    # Exponent edges separate x+1 publication from a hypothetical unrounded
    # addition; each upper call independently observes that rounded argument.
    shifts=sorted({bits(value(x)+1.)for x in inputs})
    frozen={**STAGE_BINDINGS,str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-gammaln-low-composition-questions/1',
        complete=True,inputs=inputs,upper_arguments=shifts,source_bindings=frozen,
        hypothesis='GAMMALN(x) = GAMMALN(RN53(1+x)) - LN(x) below .7; test strict and stored-PC64 subtraction using independent public calls.'))
    upper=service.ask('GAMMALN',[[client.number(x)]for x in shifts])['rows']
    if any(r['type']!='number'for r in upper):raise ValueError('numeric upper band required')
    uppers={x:r['bits']for x,r in zip(shifts,upper)}
    logs=service.ask('LN',[[client.number(x)]for x in inputs])['rows']
    target=service.ask('GAMMALN',[[client.number(x)]for x in inputs])['rows']
    precise=service.ask('GAMMALN.PRECISE',[[client.number(x)]for x in inputs])['rows']
    counts=Counter();records=[]
    for x,log,answer,other in zip(inputs,logs,target,precise):
        if any(r['type']!='number'for r in (log,answer,other)):raise ValueError('numeric low band required')
        a=uppers[bits(value(x)+1.)];b=log['bits']
        predictions={model:bits(float(operation(binary64(a),-binary64(b),'add',model)))
                     for model in ('strict','stored_pc64')}
        native_result=f'{native.publish(native.leaf("sub",load64(int(a,16)),load64(int(b,16)),0x133f)):016x}'
        if native_result!=predictions['stored_pc64']:raise ValueError('low subtraction checker disagrees')
        for name,result in predictions.items():counts[name+'_matches']+=result==answer['bits']
        counts['rows']+=1;counts['aliases_equal']+=answer==other
        counts['both_refuted']+=answer['bits']not in predictions.values()
        records.append(dict(x=x,upper_argument=bits(value(x)+1.),upper_answer=a,log=log,
                            answer=answer,precise=other,predictions=predictions))
    result=dict(format='dnacalc-gammaln-low-composition/1',complete=True,counts=counts,records=records,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Independent public calls test one recurrence and its stored argument. A contradiction '
               'does not identify a private low-band approximation, logarithm, or an unstored upper-band '
               'value. Exact aliases and all mismatches are retained; no table is fitted.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def gammaln_low_tables(checkout, output, cache, library, offline=False, count=128):
    """Challenge existing low-band tables without moving any coefficient."""
    import re
    path=Path(checkout)/'crates/calc-func/src/ported/excel_numeric/gammaln.rs'
    source=path.read_text()
    tables={}
    for name in ('P1','Q1'):
        block=source.split('const '+name+':',1)[1].split('];',1)[0]
        tables[name]=[b.replace('_','').lower()for b in re.findall(r'0x([0-9a-fA-F_]+)',block)]
        if len(tables[name])!=8:raise ValueError('low table shape changed')
    for name in ('P2','Q2'):
        block=source.split('const '+name+':',1)[1].split('];',1)[0].split('= [',1)[1]
        tables[name]=[bits(float(v.strip()))for v in block.split(',')if v.strip()]
        if len(tables[name])!=8:raise ValueError('shifted low table shape changed')
    native=Native(library)
    programs={}
    for band in ('B1','B2'):
        variants=[('continuous',False,False,False,False,0x133f),
                  ('C_Horner',False,True,True,True,0x133f),
                  ('strict',True,True,True,True,0x123f)]
        variants += [('products_final_sums_'+str(int(p))+str(int(q)),True,True,p,q,0x133f)
                     for p,q in ((False,False),(True,False),(False,True),(True,True))]
        for label,products,sums,p_final,q_final,cw in variants:
            g=Builder(1)
            def op(kind,a,b,stored):
                n=g.op(kind,a,b,store=False)
                if stored:g.nodes[n]['site']=0;g.sites=1
                return n
            t=0 if band=='B1'else op('sub',0,g.constant(bits(2.)),False)
            def chain(name,initial,last):
                if band=='B1':
                    v=g.constant(tables[name][-1])
                    coefficients=list(reversed(tables[name][:-1]))
                else:
                    v=g.constant(bits(initial));coefficients=tables[name]
                for i,c in enumerate(coefficients):
                    v=op('mul',v,t,products)
                    v=op('add',v,g.constant(c),last if i==len(coefficients)-1 else sums)
                return v
            n=chain('P'+band[-1],0.,p_final);d=chain('Q'+band[-1],1.,q_final)
            ratio=op('div',n,d,products)
            if band=='B1':
                delta=op('sub',0,g.one,False)
                result=op('mul',delta,ratio,False)
            else:
                inner=op('mul',t,ratio,False)
                inner=op('add',g.constant('3fdb0ee6072093ce'),inner,False)
                result=op('mul',t,inner,False)
            programs[band+'_'+label]=dict(graph=g.graph(result,result),cw=cw)
    centers=(.7+2.**-12,.75,.875,1.+2.**-12,1.125,1.375,1.5-2.**-12,1.625,1.75,1.875,2.-2.**-12)
    windows=[dict(center=bits(c),band='B1'if c<1.5 else'B2',
                  inputs=[f'{int(bits(c),16)-count//2+i:016x}'for i in range(count)])for c in centers]
    frozen={**STAGE_BINDINGS,str(path.resolve()):sha(path),str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-gammaln-low-table-questions/1',complete=True,
        tables=tables,windows=windows,programs=programs,source_bindings=frozen,
        selection='Fixed ordinary centers and root/band edges. Existing table bits are frozen; '
                  'programs follow continuous Horner, C assignments, strict operations, then stored products and final-sum variants.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    questions=[[client.number(x)]for w in windows for x in w['inputs']]
    answers=service.ask('GAMMALN',questions)['rows']
    precise=service.ask('GAMMALN.PRECISE',questions)['rows']
    counts=Counter();records=[];offset=0
    for window in windows:
        local=Counter()
        for x,answer,other in zip(window['inputs'],answers[offset:offset+count],precise[offset:offset+count]):
            if answer['type']!='number'or other['type']!='number':raise ValueError('numeric low table window required')
            predictions={name:f'{native.publish(native.evaluate(p["graph"],[int(x,16)],p["cw"],mask=1)[:10],p["cw"]):016x}'
                         for name,p in programs.items()if name.startswith(window['band'])}
            for name,p in predictions.items():local[name+'_matches']+=p==answer['bits']
            local['rows']+=1;local['aliases_equal']+=answer==other
            local['all_models_refuted']+=answer['bits']not in predictions.values()
            records.append(dict(x=x,band=window['band'],answer=answer,precise=other,predictions=predictions))
        offset+=count;counts.update(local)
        window['counts']=local
    result=dict(format='dnacalc-gammaln-low-tables/1',complete=True,counts=counts,windows=windows,
        records=records,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Frozen existing published tables under declared Horner stores only. These contradictions '
               'cannot distinguish missing table bits from an unknown private approximation. No coefficient, '
               'branch, or kernel is selected by match count; alias calls remain separate observations.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def beta_shape_two(checkout, output, cache, library, offline=False, count=64):
    """Isolate the series-free I_x(2,1) exponent before general beta work.

    In the retained TOMS708 bpser form, the series is zero at b=1 and
    betaln(2,1) is -ln(2). No unknown gamma table enters this control.
    """
    from exact import operation
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    centers=(2.**-100,2.**-30,2.**-10,.03125,.125,.25,.3,.375,.5-2.**-12)
    xs=sorted({f'{int(bits(c),16)-count//2+i:016x}'for c in centers for i in range(count)})
    programs={}
    for stored_log in (True,False):
        for stored_sum in (True,False):
            g=Builder(1)
            def op(kind,a,b,stored=False):
                n=g.op(kind,a,b,store=False)
                if stored:g.nodes[n]['site']=0;g.sites=1
                return n
            logarithm=op('fyl2x',op('fldln2',0,0),0,stored_log)
            exponent=op('mul',g.constant(bits(2.)),logarithm)
            exponent=op('add',exponent,g.constant(bits(math.log(2.))),stored_sum)
            base2=op('mul',exponent,op('fldl2e',0,0))
            whole=op('frndint',base2,0)
            fraction=op('sub',base2,whole)
            small=op('f2xm1',fraction,0)
            exp=op('fscale',op('add',g.one,small),whole)
            result=op('div',exp,g.constant(bits(2.)))
            label='log_'+('stored'if stored_log else'resident')+'_sum_'+('stored'if stored_sum else'resident')
            programs[label]=g.graph(result,result)
    frozen={**STAGE_BINDINGS,str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-beta-shape-two-questions/1',complete=True,
        inputs=xs,programs=programs,source_bindings=frozen,
        premise='Retained TOMS708 series-free b=1, a=2 arm: EXP(2*LN(x)+LN(2))/2. '
                'Fixed exponent/ordinary adjacent windows below .5; no fit, coefficient, or count selection.'))
    logs=service.ask('LN',[[client.number(x)]for x in xs])['rows']
    if any(r['type']!='number'for r in logs):raise ValueError('positive LN controls must be numeric')
    arguments=[bits(float(operation(2*binary64(r['bits']),binary64(bits(math.log(2.))),'add','stored_pc64')))for r in logs]
    exp_controls=service.ask('EXP',[[client.number(x)]for x in arguments])['rows']
    targets=service.ask('BETA.DIST',[[client.number(x),client.number(2.),client.number(1.),client.boolean(True)]for x in xs])['rows']
    legacy=service.ask('BETADIST',[[client.number(x),client.number(2.),client.number(1.)]for x in xs])['rows']
    records=[];counts=Counter()
    ln_graph=primitive_stage_graph('LN');exp_graph=primitive_stage_graph('EXP')
    for x,log,arg,exp,target,alias in zip(xs,logs,arguments,exp_controls,targets,legacy):
        if any(r['type']!='number'for r in (exp,target,alias)):raise ValueError('positive beta controls must be numeric')
        ln=f'{native.publish(native.evaluate(ln_graph,[int(x,16)],0x133f)[:10]):016x}'
        ee=f'{native.publish(native.evaluate(exp_graph,[int(arg,16)],0x133f)[:10]):016x}'
        if ln!=log['bits']or ee!=exp['bits']:raise ValueError('independent public primitive differs from native instruction model')
        predictions={name:f'{native.publish(native.evaluate(graph,[int(x,16)],0x133f,mask=1)[:10]):016x}'for name,graph in programs.items()}
        paired_prediction=bits(value(exp['bits'])/2.)
        if predictions['log_stored_sum_stored']!=paired_prediction:raise ValueError('paired public beta exponent differs from native')
        predictions['stored_square']=bits(float(operation(binary64(x),binary64(x),'mul','stored_pc64')))
        for label,p in predictions.items():counts[label+'_matches']+=p==target['bits']
        counts['rows']+=1;counts['aliases_equal']+=target==alias
        counts['all_models_refuted']+=target['bits']not in predictions.values()
        records.append(dict(x=x,log=log,exponent_argument=arg,exp=exp,answer=target,legacy=alias,predictions=predictions))
    result=dict(format='dnacalc-beta-shape-two/1',complete=True,counts=counts,records=records,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Conditional series-free exponent controls and independently verified public primitives. '
               'The square is a mathematical control, not a proposed special case. A complete match '
               'would identify only this sampled arm; other beta normalizers, series and callers remain separate.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def gammaln_log_controls(checkout, replay, output, cache, offline=False):
    """Use primitive-log disagreements before asking the low-band target."""
    from exact import operation
    doc=json.loads(Path(replay).read_text())
    if not doc['complete']or doc['format']!='dnacalc-gammaln-low-stage-replay/1':
        raise ValueError('complete current stage replay required')
    selected=[r for r in doc['records']if r['runtime_log']!=r['worksheet_log']]
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    uppers=service.ask('GAMMALN',[[client.number(r['upper_argument'])]for r in selected])['rows']
    logs=service.ask('LN',[[client.number(r['original']['x'])]for r in selected])['rows']
    controls=[];indistinguishable=[]
    for r,upper,log in zip(selected,uppers,logs):
        if upper['type']!='number'or log['type']!='number':raise ValueError('primitive controls must be numeric')
        if log['bits']!=r['worksheet_log']:raise ValueError('public logarithm differs from current stage replay')
        predictions={name+'_'+model:bits(float(operation(binary64(upper['bits']),-binary64(b),'add',model)))
                     for name,b in (('runtime',r['runtime_log']),('worksheet',r['worksheet_log']))
                     for model in ('strict','stored_pc64')}
        record=dict(x=r['original']['x'],upper=upper,log=log,replay=r,predictions=predictions)
        if len(set(predictions.values()))>1:controls.append(record)
        else:indistinguishable.append(record)
    frozen={**STAGE_BINDINGS,str(Path(replay).resolve()):sha(replay)}
    questions=dict(format='dnacalc-gammaln-log-controls-questions/1',complete=True,
        replay_rows=len(doc['records']),primitive_log_disagreements=len(selected),controls=controls,
        indistinguishable=indistinguishable,source_bindings=frozen,
        selection='All primitive runtime/worksheet LN disagreements in the fixed prior RNG input set; '
                  'independent public upper-band values distinguish subtraction predictions before any low-band target is requested.')
    write(str(output)+'.questions.json',questions)
    answers=service.ask('GAMMALN',[[client.number(r['x'])]for r in controls])['rows']
    aliases=service.ask('GAMMALN.PRECISE',[[client.number(r['x'])]for r in controls])['rows']
    records=[];counts=Counter()
    for r,answer,alias in zip(controls,answers,aliases):
        if answer['type']!='number'or alias['type']!='number':raise ValueError('finite low band target required')
        for name,p in r['predictions'].items():counts[name+'_matches']+=p==answer['bits']
        counts['rows']+=1;counts['aliases_equal']+=answer==alias
        counts['all_models_refuted']+=answer['bits']not in r['predictions'].values()
        records.append(dict(r,answer=answer,alias=alias))
    result=dict(format='dnacalc-gammaln-log-controls/1',complete=True,counts=counts,records=records,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='A sampled logarithm consumer discrimination with independent public upper-band and LN '
               'observations. The upper approximation remains unresolved; this does not recover its table. '
               'Strict and PC64 subtraction are reported separately even when indistinguishable.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def beta_log_replay(study, replay, output, library):
    """Separate the stored runtime LN leaf from exponent publication."""
    import copy
    doc=json.loads(Path(study).read_text())
    questions_path=str(study)+'.questions.json'
    if sha(questions_path)!=doc['questions_sha256']:raise ValueError('beta questions changed')
    questions=json.loads(Path(questions_path).read_text())
    leaves=json.loads(Path(replay).read_text())
    if not leaves['complete']:raise ValueError('incomplete primitive replay')
    runtime={r['original']['x']:r['runtime_log']for r in leaves['records']}
    native=Native(library);graphs={}
    for name,original in questions['programs'].items():
        if not name.startswith('log_stored_'):continue
        graph=copy.deepcopy(original)
        found=[n for n in graph['nodes']if n['op']=='fyl2x']
        if len(found)!=1:raise ValueError('changed beta logarithm graph')
        found[0].update(op='input',a=1)
        graphs['runtime_'+name]=graph
        resident_normalizer=copy.deepcopy(graph)
        constants=[n for n in resident_normalizer['nodes']if n['op']=='const'and n.get('bits')==bits(math.log(2.))]
        if len(constants)!=1:raise ValueError('changed beta normalizer graph')
        constants[0]['op']='fldln2'
        graphs['runtime_'+name+'_resident_ln2']=resident_normalizer
    records=[];counts=Counter()
    for r in doc['records']:
        if r['x']not in runtime:raise ValueError('missing runtime logarithm')
        predictions={name:f'{native.publish(native.evaluate(graph,[int(r["x"],16),int(runtime[r["x"]],16)],0x133f,mask=1)[:10]):016x}'
                     for name,graph in graphs.items()}
        for name,p in predictions.items():counts[name+'_matches']+=p==r['answer']['bits']
        counts['rows']+=1
        counts['runtime_public_log_different']+=runtime[r['x']]!=r['log']['bits']
        counts['all_models_refuted']+=r['answer']['bits']not in predictions.values()
        records.append(dict(original=r,runtime_log=runtime[r['x']],predictions=predictions))
    result=dict(format='dnacalc-beta-runtime-log-replay/1',complete=True,counts=counts,records=records,
        graphs=graphs,source_bindings={**STAGE_BINDINGS,str(Path(study).resolve()):sha(study),
            str(Path(replay).resolve()):sha(replay),str(Path(library).resolve()):sha(library)},
        limits='The existing cleared runtime logarithm is a declared hypothesis, not an Excel '
               'observation. Fixed series-free beta exponent programs and every contradiction remain '
               'visible. No normalizer or coefficient is fitted and no kernel is selected by counts.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def beta_integer_unit_shape(checkout, output, cache, library, offline=False, count=64, prior=None):
    """Separate the published b=1 normalizer loop from its EXP consumer.

    Integer shapes make the other gamma helpers zero. The retained TOMS708
    reduction is a short product of b/(1+b), then a logarithm. This is an
    algorithm-derived control, never an integer-shape replacement kernel.
    """
    import itertools
    from fractions import Fraction as Q
    from exact import operation
    native=Native(library);client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    shapes=(2,3,4,5,6,7)
    centers=(2.**-10,.03125,.125,.25,.3,.375,.5-2.**-12)
    xs=sorted({f'{int(bits(c),16)-count//2+i:016x}'for c in centers for i in range(count)})
    loops={}
    for a,ratio_store,product_store in itertools.product(shapes,(True,False),(True,False)):
        product=load64(int(bits(1.),16));rproduct=Q(1)
        for b in range(a-1,0,-1):
            ratio=native.leaf('div',load64(int(bits(b),16)),load64(int(bits(1+b),16)),0x133f)
            rr=operation(Q(b),Q(1+b),'div','stored_pc64'if ratio_store else'resident_pc64')
            if ratio_store:ratio=load64(native.publish(ratio))
            product=native.leaf('mul',product,ratio,0x133f)
            rproduct=operation(rproduct,rr,'mul','stored_pc64'if product_store else'resident_pc64')
            if product_store:product=load64(native.publish(product))
        if decode(product)!=rproduct:raise ValueError('beta reduction differs from rational checker')
        loops[(a,ratio_store,product_store)]=product
    graphs={}
    for a,mul_store,norm_store,sum_store,exp_store,reciprocal in itertools.product(shapes,*([(True,False)]*5)):
        g=Builder(2)
        def op(kind,a,b=0,stored=False):
            n=g.op(kind,a,b,store=False)
            if stored:g.nodes[n]['site']=0;g.sites=1
            return n
        normalizer=op('fyl2x',op('fldln2',0),1,norm_store)
        term=op('mul',g.constant(bits(a)),0,mul_store)
        exponent=op('sub',term,normalizer,sum_store)
        base2=op('mul',exponent,op('fldl2e',0))
        if reciprocal:base2=op('neg',base2)
        whole=op('frndint',base2)
        fraction=op('sub',base2,whole)
        exp=op('fscale',op('add',g.one,op('f2xm1',fraction)),whole)
        if reciprocal:exp=op('div',g.one,exp)
        if exp_store:exp=op('add',exp,g.constant(bits(0.)),True)
        result=op('div',exp,g.constant(bits(a)))
        graphs[(a,mul_store,norm_store,sum_store,exp_store,reciprocal)]=g.graph(result,result)
    labels=list(itertools.product((True,False),repeat=7))
    def label(p):
        return '_'.join(n+('S'if v else'R')for n,v in zip(('ratio','product','mul','normalizer','sum','exp'),p[:6]))+('_coreReciprocal'if p[6]else'_coreDirect')
    frozen={**STAGE_BINDINGS,str(Path(library).resolve()):sha(library)}
    manifest=dict(format='dnacalc-beta-integer-unit-shape-questions/1',complete=True,
        shapes=shapes,inputs=xs,program_order=[label(p)for p in labels],
        normalizer_arguments={','.join(map(str,k)):v.hex()for k,v in loops.items()},
        graphs={','.join(map(str,k)):v for k,v in graphs.items()},source_bindings=frozen,
        scope='The retained betaln reduction arm has min(a,b)=1 and max(a,b)<8. Earlier exploratory shapes 8,16,31 leave that arm and cannot refute its reduction loop.',
        selection='Fixed integer shapes and adjacent windows below .5. Published TOMS708 reduction '
                  'and finite store/EXP programs frozen before target calls. No integer shortcut is proposed.')
    if prior:manifest['reused_study']=dict(path=str(Path(prior).resolve()),sha256=sha(prior))
    write(str(output)+'.questions.json',manifest)
    logs=service.ask('LN',[[client.number(x)]for x in xs])['rows']
    existing={}
    if prior:
        previous=json.loads(Path(prior).read_text())
        if sha(str(prior)+'.questions.json')!=previous['questions_sha256']:raise ValueError('reused beta questions changed')
        existing={(r['a'],r['x']):r['answer']for r in previous['records']}
    pending=[(a,x)for a in shapes for x in xs if(a,x)not in existing]
    if pending:
        questions=[[client.number(x),client.number(a),client.number(1.),client.boolean(True)]for a,x in pending]
        original=service.ask('BETA.DIST',questions)['rows']
        existing.update(zip(pending,original))
    answers=[existing[(a,x)]for a in shapes for x in xs]
    counts={label(p):Counter()for p in labels};records=[];matrix=hashlib.sha256()
    ln_graph=primitive_stage_graph('LN')
    for index,(a,x,log,answer)in enumerate((a,x,logs[j],answers[i*len(xs)+j])for i,a in enumerate(shapes)for j,x in enumerate(xs)):
        if log['type']!='number'or answer['type']!='number':raise ValueError('finite beta controls required')
        actual_log=f'{native.publish(native.evaluate(ln_graph,[int(x,16)],0x133f)[:10]):016x}'
        if actual_log!=log['bits']:raise ValueError('beta control LN differs from independent primitive')
        mask=0
        for i,p in enumerate(labels):
            ratio_store,product_store,*outer=p
            raw=f'{native.publish(native.evaluate(graphs[(a,*outer)],[int(log["bits"],16)],0x133f,mask=1,extended=loops[(a,ratio_store,product_store)])[:10]):016x}'
            matrix.update(bytes.fromhex(raw))
            tally=counts[label(p)];tally['rows']+=1
            if raw==answer['bits']:tally['matched']+=1;mask|=1<<i
            else:
                tally['failed']+=1
                # One first contradiction suffices to reject a complete tree;
                # all answers and the prediction matrix digest remain bound.
                if 'first_contradiction'not in tally:tally['first_contradiction']=dict(index=index,x=x,a=a,prediction=raw,answer=answer)
        records.append(dict(x=x,a=a,answer=answer,matching_program_mask=f'{mask:032x}'))
    result=dict(format='dnacalc-beta-integer-unit-shape/1',complete=True,counts=counts,records=records,
        prediction_matrix_sha256=matrix.hexdigest(),receipts=service.receipts,source_bindings=frozen,
        reused_study=manifest.get('reused_study'),new_target_calls=len(pending),reused_target_calls=len(answers)-len(pending),
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Fixed published-loop and arithmetic programs. Complete models require every observation; '
               'partial scores do not choose a kernel or a per-row selector. The exponent, normalizer '
               'and EXP publication are distinct hypotheses. General nonwhole shapes and the series remain open.')
    write(output,result)
    print(json.dumps(dict(rows=len(records),complete_programs=[k for k,c in counts.items()if not c['failed']],
                          no_program_matches=sum(r['matching_program_mask']=='0'*32 for r in records))),flush=True)
    return result


def beta_f_unit_prefix(checkout, prior, output, cache, library, offline=False, power_pairs=False):
    """Reuse the established b=1 body to question the F input transform."""
    import itertools
    from fractions import Fraction as Q
    from exact import operation
    from x87 import encode
    prior=Path(prior);study=json.loads(prior.read_text())
    question_path=Path(str(prior)+'.questions.json')
    if sha(question_path)!=study['questions_sha256']:
        raise ValueError('established beta questions changed')
    questions=json.loads(question_path.read_text())
    for core in ('Reciprocal','Direct'):
        name='ratioR_productR_mulS_normalizerR_sumR_expR_core'+core
        if study['counts'][name].get('failed')or study['counts'][name]['matched']!=len(study['records']):
            raise ValueError('complete prior beta body required')
    native=Native(library);client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    cases=[]
    for a,scale in itertools.product(range(2,8),(1.125,2.,4.,16.,256.)):
        center=a*scale
        cases.extend(dict(a=a,x=f'{int(bits(center),16)-32+i:016x}')for i in range(64))
    profiles=list(itertools.product(('strict','stored_pc64','resident_pc64'),repeat=2))
    profiles=list(itertools.product(profiles,(True,False),(True,False)))
    def label(p):return p[0][0]+'_'+p[0][1]+('_logS'if p[1]else'_logR')+('_coreReciprocal'if p[2]else'_coreDirect')
    frozen={**STAGE_BINDINGS,str(prior.resolve()):sha(prior),str(question_path.resolve()):sha(question_path),str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-beta-f-unit-prefix-questions/1',complete=True,
        cases=cases,programs=[label(p)for p in profiles],source_bindings=frozen,power_pairs=power_pairs,
        premise='For F(x,2,2a) right tail the b=1 beta input is z=a/(a+x). '
                'Reuse only the previously established resident normalizer and stored a*LN stage. '
                'Question denominator, quotient and LN lifetime independently; no partial score selects a model. '
                'Optional public POWER pairs test the elementary ratio^a and (1+x/a)^(-a) forms on these exposed controls.'))
    args=[[client.number(c['x']),client.number(2.),client.number(2*c['a'])]for c in cases]
    answers=service.ask_rows('FDIST',args);aliases=service.ask_rows('F.DIST.RT',args)
    powers={}
    if power_pairs:
        for form,word in itertools.product(('ratio','one_plus'),('strict','stored_pc64')):
            arguments=[]
            for case in cases:
                a=Q(case['a']);x=binary64(case['x'])
                if form=='ratio':base=operation(a,operation(a,x,'add',word),'div',word);exponent=a
                else:base=operation(Q(1),operation(x,a,'div',word),'add',word);exponent=-a
                arguments.append([client.number(bits(float(base))),client.number(float(exponent))])
            powers[form+'_'+word]=service.ask_rows('POWER',arguments)
    one=load64(int(bits(1.),16));ln2=native.leaf('fldln2',one,one,0x133f)
    counts={label(p):Counter()for p in profiles};records=[]
    for index,(case,answer,alias)in enumerate(zip(cases,answers,aliases)):
        if answer.get('type')!='number'or alias!=answer:
            raise ValueError('finite matching F-tail aliases required')
        a=case['a'];x=binary64(case['x']);predictions={}
        normalizer=bytes.fromhex(questions['normalizer_arguments'][f'{a},False,False'])
        for profile in profiles:
            (dw,qw),log_store,reciprocal=profile
            den=operation(Q(a),x,'add',dw);z=operation(Q(a),den,'div',qw)
            logarithm=native.leaf('fyl2x',ln2,encode(z),0x133f)
            if log_store:logarithm=load64(native.publish(logarithm))
            graph=questions['graphs'][f'{a},True,False,False,False,{reciprocal}']
            raw=native.evaluate(graph,[],0x133f,mask=1,extended=logarithm+normalizer)
            predicted=f'{native.publish(raw[:10]):016x}';name=label(profile)
            predictions[name]=predicted;tally=counts[name];tally['rows']+=1
            if predicted==answer['bits']:tally['matched']+=1
            else:
                tally['failed']+=1
                if 'first_contradiction'not in tally:tally['first_contradiction']=dict(index=index,case=case,answer=answer,prediction=predicted)
        records.append(dict(case=case,answer=answer,alias=alias,predictions=predictions))
    power_counts={}
    for name,replies in powers.items():
        if any(reply.get('type')!='number'for reply in replies):raise ValueError('finite public powers required')
        power_counts[name]=sum(reply==answer for reply,answer in zip(replies,answers))
        for record,reply in zip(records,replies):record.setdefault('public_powers',{})[name]=reply
    result=dict(format='dnacalc-beta-f-unit-prefix/1',complete=True,counts=counts,records=records,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        public_power_counts=power_counts,
        complete_programs=[k for k,v in counts.items()if not v['failed']],
        limits='Conditional declared F transform and previously identified b=1 beta body only. '
               'These controls identify no general beta series, caller or per-row choice.')
    write(output,result);print(json.dumps(dict(rows=len(records),complete_programs=result['complete_programs'],public_power_counts=power_counts)),flush=True)
    return result


def beta_public_f_tails(checkout, trace, output, cache, offline=False):
    """Pair the public F surfaces before attributing their failures to beta."""
    trace=Path(trace)
    original=[json.loads(line)for line in trace.read_text().splitlines()]
    selected=[r for r in original if r['actual']!=r['expected']]
    inputs=[tuple(r['args'])for r in selected]
    for center,d1,d2 in ((.22,2.,19.),(.22,5.,11.),(.5,1.,1.),(1.,2.,2.),(1.,4.,4.),(4.,7.,12.)):
        inputs.extend(('0x'+f'{int(bits(center),16)-32+i:016x}','0x'+bits(d1),'0x'+bits(d2))for i in range(64))
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    for row in inputs:
        if len(row)!=3 or any(not math.isfinite(value(x.removeprefix('0x')))or value(x.removeprefix('0x'))<0 for x in row):
            raise ValueError('finite positive F controls required')
    frozen={**STAGE_BINDINGS,str(trace.resolve()):sha(trace)}
    write(str(output)+'.questions.json',dict(format='dnacalc-beta-public-f-tail-questions/1',complete=True,
        inputs=inputs,selected_bank_rows=len(selected),source_bindings=frozen,
        selection='Every previously observed FDIST bank failure, then six fixed64-adjacent windows. '
                  'These are exposed diagnostic controls, not blind validation. Compare aliases, historical '
                  'expectations and one-minus-public-CDF before changing the shared beta evaluator.'))
    args=[[client.number(x.removeprefix('0x'))for x in row]for row in inputs]
    old=service.ask_rows('FDIST',args)
    right=service.ask_rows('F.DIST.RT',args)
    cdf=service.ask_rows('F.DIST',[row+[client.boolean(True)]for row in args])
    counts=Counter();records=[]
    for i,(row,a,b,c)in enumerate(zip(inputs,old,right,cdf)):
        if any(v.get('type')!='number'for v in(a,b,c)):
            raise ValueError('numeric F-tail controls required')
        counts['rows']+=1;counts['public_aliases_equal']+=a==b
        complement=bits(1.-value(c['bits']))
        counts['one_minus_cdf_matches']+=a['bits']==complement
        record=dict(args=row,fdist=a,right_tail=b,cdf=c,one_minus_cdf=complement)
        if i<len(selected):
            r=selected[i];record['original_bank']=r
            repeat='0x'+a['bits']
            counts['bank_rows']+=1
            counts['old_expected_repeated']+=repeat==r['expected']
            counts['old_actual_repeated']+=repeat==r['actual']
            counts['neither_repeated']+=repeat not in(r['expected'],r['actual'])
        records.append(record)
    result=dict(format='dnacalc-beta-public-f-tails/1',complete=True,counts=counts,records=records,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Original public pairs and exposed historical failures only. Alias agreement does not '
               'identify a private beta producer, and a changed historical answer does not score as a '
               'canonical bank pass. No context/profile is inferred or installed from these observations.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def beta_series_normalizer(checkout, output, cache, library, offline=False, count=128, resident_log=False, strict_final=False):
    """Invert I_x(2,2)'s one-term series into a shared normalizer interval."""
    import itertools
    from fractions import Fraction as Q
    from exact import operation
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library)
    centers=(2.**-10,.03125,.125,.25,.3)
    xs=sorted({f'{int(bits(c),16)-count//2+i:016x}'for c in centers for i in range(count)})
    graphs={}
    for reciprocal in (True,False):
        g=Builder(1)
        # The positive inverse coordinate is e=-t, since all these
        # prefactors have negative exponents. Output remains resident.
        t=g.op('mul',0,g.op('fldl2e',store=False),store=False)
        if not reciprocal:t=g.op('neg',t,store=False)
        whole=g.op('frndint',t,store=False)
        fraction=g.op('sub',t,whole,store=False)
        small=g.op('f2xm1',fraction,store=False)
        exp=g.op('fscale',g.op('add',g.one,small,store=False),whole,store=False)
        if reciprocal:exp=g.op('div',g.one,exp,store=False)
        graphs[reciprocal]=g.graph(exp,exp)
    specifications=list(itertools.product((True,False),repeat=5))
    def label(p):return '_'.join(n+('S'if s else'R')for n,s in zip(('series','product','factor','prefactor'),p[:4]))+('_coreReciprocal'if p[4]else'_coreDirect')
    frozen={**STAGE_BINDINGS,str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-beta-series-normalizer-questions/1',complete=True,
        inputs=xs,graphs={str(k):v for k,v in graphs.items()},program_order=[label(p)for p in specifications],
        source_bindings=frozen,resident_log=resident_log,strict_final=strict_final,selection='Fixed I_x(2,2) bpser windows below .3. '
            'The series is exactly one term: 1-2*x/3, evaluated under declared stores. '
            'No unknown normalizer value is fitted; exact intervals test for one shared resident value.'))
    logs=service.ask('LN',[[client.number(x)]for x in xs])['rows']
    answers=service.ask('BETA.DIST',[[client.number(x),client.number(2.),client.number(2.),client.boolean(True)]for x in xs])['rows']
    models={label(p):dict(common=None,contradiction=None,applied=0)for p in specifications}
    records=[]
    for index,(x,log,answer)in enumerate(zip(xs,logs,answers)):
        if log['type']!='number'or answer['type']!='number':raise ValueError('finite beta series controls required')
        logarithm=native.evaluate(primitive_stage_graph('LN'),[int(x,16)],0x133f)[:10]
        if f'{native.publish(logarithm):016x}'!=log['bits']:
            raise ValueError('public LN primitive differs from native model')
        term=2*(decode(logarithm)if resident_log else binary64(log['bits']))
        center=Q(-math.log(6.))-term
        records.append(dict(x=x,log=log,answer=answer))
        for p in specifications:
            model=models[label(p)]
            if model['contradiction']is not None:continue
            ws,ps,fs,bs,reciprocal=p
            w=operation(-binary64(x),Q(3),'div','stored_pc64'if ws else'resident_pc64')
            product=operation(Q(2),w,'mul','stored_pc64'if ps else'resident_pc64')
            factor=operation(Q(1),product,'add','stored_pc64'if fs else'resident_pc64')
            y=binary64(answer['bits'])
            domain=inverse(Span(y,y),factor,'mul',True,'strict'if strict_final else'stored_pc64')
            if bs:domain=BINARY64.pullback(domain)
            domain=domain.linear(Q(2))if domain else None
            cell=positive_extended_cell(native,graphs[reciprocal],None,center,decreasing=True,
                                        max_radius_bits=24,output_span=domain)if domain else None
            normalizer=None
            if cell:
                exponent=Span(decode(bytes.fromhex(cell['low'])),decode(bytes.fromhex(cell['high'])))
                normalizer=inverse(exponent,-term,'add',True,'resident_pc64')
            model['applied']+=1
            if model['common']is None and model['applied']==1:model['common']=normalizer
            else:model['common']=model['common'].meet(normalizer)if model['common']else None
            if EXTENDED64.coordinates(model['common'])is None:
                model['contradiction']=dict(index=index,x=x,answer=answer,
                    current_normalizer_preimage=normalizer.json()if normalizer else None)
    # Independent published reduction candidates, without moving any digit.
    one=load64(int(bits(1.),16));three=load64(int(bits(3.),16));half=load64(int(bits(.5),16))
    ln2=native.leaf('fldln2',one,one,0x133f)
    w=native.leaf('div',one,three,0x133f)
    lw=native.leaf('fyl2x',ln2,w,0x133f);lz=native.leaf('fyl2x',ln2,half,0x133f)
    reduction=native.leaf('add',lw,lz,0x133f)
    ideal=native.leaf('neg',native.leaf('fyl2x',ln2,load64(int(bits(6.),16)),0x133f),one,0x133f)
    candidates=dict(resident_reduction=decode(reduction),resident_negative_ln6=decode(ideal),stored_negative_ln6=Q(-math.log(6.)))
    summaries={}
    for name,m in models.items():
        span=m['common'];e=EXTENDED64.coordinates(span);d=BINARY64.coordinates(span)
        summaries[name]=dict(applied_constraints=m['applied'],contradiction=m['contradiction'],
            shared_real_interval=span.json()if span else None,resident_count=e[1]-e[0]+1 if e else 0,
            binary64_count=d[1]-d[0]+1 if d else 0,
            contained_candidates={k:bool(span and span.contains(v))for k,v in candidates.items()})
    result=dict(format='dnacalc-beta-series-normalizer/1',complete=True,records=records,models=summaries,
        candidate_values={k:str(v)for k,v in candidates.items()},receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),resident_log=resident_log,strict_final=strict_final,
        limits='Necessary exact constraints under the declared resident EXP monotonicity and selected stored/resident-LN '
               'premises. Empty intersections refute whole programs. A nonempty interval does not identify '
               'a producer. Refuted models stop at their first contradiction; all original target rows remain '
               'retained. Prefactor storage and pre-division EXP storage are indistinguishable at shape2.')
    write(output,result)
    print(json.dumps(dict(rows=len(records),feasible={k:v for k,v in summaries.items()if v['resident_count']})),flush=True)
    return result


def gammaln_stirling_shape(checkout):
    """Expose the existing stored Stirling tree to the plateau selector.

    Input 1 is a separately measured LN leaf. Coefficients are read from
    the application table, never copied or selected by answer counts.
    """
    import re
    from fractions import Fraction as Q
    directory=Path(checkout)/'crates/calc-func/src/ported/excel_numeric'
    chains=directory/'x87_chains.rs';kernel=directory/'gammaln.rs'
    block=chains.read_text().split('const STIRLING_TERMS:',1)[1].split('];',1)[0].split('= [',1)[1]
    coefficients=[]
    for text in block.split(','):
        if not text.strip():continue
        numerator,denominator=text.strip().split('/')
        coefficients.append(bits(float(Q(numerator.strip())/Q(denominator.strip()))))
    if len(coefficients)!=10:raise ValueError('Stirling table changed')
    ls=re.search(r'const LS2PI:.*?0x([0-9A-Fa-f_]+)',kernel.read_text())[1].replace('_','').lower()
    g=Builder(2);constants={}
    def constant(name,raw):
        node=g.constant(raw);constants[name]=node;return node
    q=g.op('mul',g.op('sub',0,g.constant(bits(.5))),1)
    q=g.op('sub',q,0);q=g.op('add',q,constant('LS2PI',ls))
    z=g.op('div',g.one,0);t=g.op('mul',z,z)
    w=constant('C9',coefficients[-1])
    for index in range(8,-1,-1):
        w=g.op('add',g.op('mul',w,t),constant('C'+str(index),coefficients[index]))
    out=g.op('add',q,g.op('mul',z,w))
    return dict(format='dnacalc-plateau-shape/1',function='GAMMALN',argument_inputs=[0],
        graph=g.graph(out,out),constants=constants,cw=0x133f,mask=(1<<g.sites)-1,
        domain=dict(low=bits(8.),high=bits(128.)),
        source_bindings={str(p.resolve()):sha(p)for p in(chains,kernel)},
        leaf='Input1 is independently observed LN(x). Finite no-separator results leave constant bits undetermined.')


def gammaln_stirling_proof(checkout, output, cache, library, count=2048, seed=44005, extra_inputs=None, offline=False):
    """Freeze a fresh stored-program prediction with audited service leaves."""
    import random
    import numpy as np
    shape=gammaln_stirling_shape(checkout);native=Native(library)
    rng=random.Random(seed);xs=[bits(2.**rng.uniform(3.,7.))for _ in range(count)]
    if extra_inputs:
        additional=json.loads(Path(extra_inputs).read_text())
        if any(not isinstance(x,str)or not 8.<=value(x)<=128. for x in additional):
            raise ValueError('extra Stirling controls must be in [8,128]')
        xs=list(dict.fromkeys(xs+additional))
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=90),cache,offline=offline,priority='interactive',allow_unread_inputs=True)
    logs=service.ask('LN',[[client.number(x)]for x in xs])['rows']
    if any(row.get('type')!='number'for row in logs):raise ValueError('numeric LN leaves required')
    inputs=np.array([[int(x,16),int(row['bits'],16)]for x,row in zip(xs,logs)],dtype=np.uint64)
    predictions=[f'{int(v):016x}'for v in native.grid(shape['graph'],inputs,shape['cw'],shape['mask'])[:,0]]
    frozen=Path(str(output)+'.questions.json')
    if frozen.exists():raise FileExistsError('fresh proof requires a new output path')
    write(frozen,dict(xs=xs,shape=shape,predictions=predictions,seed=seed,random_count=count,
        extra_inputs_sha256=sha(extra_inputs)if extra_inputs else None,
        leaf_receipts=service.receipts.copy(),source_bindings=STAGE_BINDINGS,
        library_sha256=sha(library)))
    answers=service.ask('GAMMALN',[[client.number(x)]for x in xs])['rows']
    matched=sum(row==dict(type='number',bits=p)for row,p in zip(answers,predictions))
    result=dict(complete=True,offline=offline,rows=len(xs),matched=matched,failed=len(xs)-matched,
        records=[dict(x=x,prediction=p,answer=a)for x,p,a in zip(xs,predictions,answers)],
        receipts=service.receipts,questions_sha256=sha(frozen),
        limits='A fresh proof of the fixed stored program with independently measured public LN leaves, not unique identification of every store site. Every unread or nonnumeric row fails.')
    write(output,result);print(json.dumps(dict(rows=len(xs),matched=matched,failed=len(xs)-matched)),flush=True)
    return result


def gammaln_table_frequencies(windows, tables_path, output):
    """Compare frozen-table stage frequencies with independent output spectra."""
    from fractions import Fraction as Q
    from exact import exponent_of,two
    observed=json.loads(Path(windows).read_text())
    tables=json.loads(Path(tables_path).read_text())['tables']
    records=[]
    for fn in observed['functions']:
        for window in fn['windows']:
            x=value(f'{int(window["start"],16)+window["count"]//2:016x}')
            xx=Q(x);dx=Q(math.ulp(x));band='B1'if x<1.5 else'B2'
            variable=xx if band=='B1'else xx-2
            stages=[];final={}
            def emit(name,v,d):
                if not v:return
                frequency=d*dx/two(exponent_of(abs(v))-52)
                folded=frequency%1
                if folded>Q(1,2):folded=1-folded
                stages.append(dict(stage=name,value=str(v),derivative=str(d),cycles_per_input=float(folded)))
            for name in ('P','Q'):
                coefficients=[binary64(c)for c in tables[name+band[-1]]]
                if band=='B1':v=coefficients[-1];d=Q(0);coefficients=list(reversed(coefficients[:-1]))
                else:v=Q(0 if name=='P'else 1);d=Q(0)
                for i,c in enumerate(coefficients):
                    derivative=v+variable*d
                    v*=variable;d=derivative;emit(name+'_product_'+str(i),v,d)
                    v+=c;emit(name+'_sum_'+str(i),v,d)
                final[name]=(v,d)
            p,pd=final['P'];q,qd=final['Q'];r=p/q;rd=(pd*q-p*qd)/(q*q)
            emit('ratio',r,rd)
            if band=='B1':emit('result',(xx-1)*r,r+(xx-1)*rd)
            else:
                anchor=binary64('3fdb0ee6072093ce')
                inner=anchor+variable*r
                emit('inner_product',variable*r,r+variable*rd)
                emit('inner_sum',inner,r+variable*rd)
                emit('result',variable*inner,inner+variable*(r+variable*rd))
            peaks=[]
            for line in window['spectrum']['lines']:
                candidates=sorted((abs(s['cycles_per_input']-line['cycles_per_input'])*window['count'],s['stage'])for s in stages)
                peaks.append(dict(line,nearest_stage=candidates[0][1],distance_in_dft_bins=candidates[0][0],
                                  compatible_stages=[name for distance,name in candidates if distance<=1.5]))
            records.append(dict(function=fn['function'],center=bits(x),band=band,stages=stages,peaks=peaks,
                                unmapped_peaks=sum(not p['compatible_stages']for p in peaks)))
    result=dict(format='dnacalc-gammaln-table-frequencies/1',complete=True,records=records,
        source_bindings={str(Path(windows).resolve()):sha(windows),str(Path(tables_path).resolve()):sha(tables_path),**STAGE_BINDINGS},
        limits='Ideal local derivatives of the frozen existing tables, not private Excel stage values. '
               'The fixed 1.5-bin proximity is spectral resolution, not a verdict tolerance. Missing '
               'matches challenge this table/coordinate model; coincident lines do not prove a store, '
               'coefficient or phase. Harmonics and interactions may supply additional output lines.')
    write(output,result)
    print(json.dumps([dict(function=r['function'],band=r['band'],unmapped_peaks=r['unmapped_peaks'],
                           peak_count=len(r['peaks']))for r in records]),flush=True)
    return result


def audit_gammaln_log_change(controls, registered, output):
    """Bind the logarithm change and every new loss to actual invocation."""
    doc=json.loads(Path(controls).read_text())
    actual=json.loads(Path(registered).read_text())
    if not doc['complete']or not actual['complete']:raise ValueError('complete comparisons required')
    observed=doc['records'];computed=actual['records']
    if len(observed)!=len(computed):raise ValueError('registered cohort changed')
    counts=Counter();losses=[];records=[]
    for r,calculated in zip(observed,computed):
        if calculated['original']['x']!=r['x']:raise ValueError('registered operands reordered')
        old=r['replay']['registered'].removeprefix('0x')
        new=calculated['registered'].removeprefix('0x');expected=r['answer']['bits']
        if old!=r['replay']['runtime_subtract']:raise ValueError('old registered replay differs')
        if new!=calculated['worksheet_subtract']:raise ValueError('new registered replay differs')
        if calculated['upper']!=r['replay']['upper']:raise ValueError('upper approximation changed')
        for name,v in (('old_matched',old==expected),('new_matched',new==expected),
                       ('gained',old!=expected and new==expected),('lost',old==expected and new!=expected),
                       ('changed',old!=new)):
            counts[name]+=v
        counts['rows']+=1
        records.append(dict(x=r['x'],expected=expected,old_registered=old,new_registered=new))
        if old==expected and new!=expected:
            y=binary64(expected);ln=binary64(r['log']['bits'])
            needed=inverse(Span(y,y),-ln,'add',True,'strict')
            current=binary64(calculated['upper']);upper=binary64(r['upper']['bits'])
            if needed is None or needed.contains(current)or not needed.contains(upper):
                raise ValueError('loss is not classified by the independent upper-band constraint')
            losses.append(dict(original=r,calculated=calculated,
                required_upper_preimage=needed.json(),classification='unresolved upper-band table',
                proof='The actual upper value is outside the exact final-subtraction preimage; '
                      'the independently observed public upper value is inside it and reproduces the target.'))
    result=dict(format='dnacalc-gammaln-log-change-audit/1',complete=True,counts=counts,records=records,
        classified_losses=losses,source_bindings={**STAGE_BINDINGS,str(Path(controls).resolve()):sha(controls),
            str(Path(registered).resolve()):sha(registered)},
        limits='All before/after values are registered application returns. Every new-cohort loss is '
               'retained and independently localized to the unresolved upper table under the measured '
               'public composition. The new cohort is separate from the canonical bank and corpus counts.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def density_normalizer_questions(checkout, bank, output, cache, library, offline=False):
    """Pair the remaining low density failures with public GAMMALN readings."""
    native=Native(library);client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    cases=[];omitted=Counter()
    for line in Path(bank).read_text().splitlines():
        r=json.loads(line)
        if r['expected']==r['actual']:continue
        a=r['args']
        flag=a[3]if isinstance(a[3],bool)else value(a[3].removeprefix('0x'))!=0
        if flag:omitted['cumulative']+=1;continue
        if not all(isinstance(v,str)and v.startswith('0x')for v in a[:3]):
            omitted['nonnumeric_arguments']+=1;continue
        x,alpha,beta=[v.removeprefix('0x')for v in a[:3]]
        xx,aa,bb=map(value,(x,alpha,beta))
        if not(xx>0 and 0<aa<7 and aa!=1 and bb>0):omitted['outside_low_shape_scope']+=1;continue
        y=native.leaf('div',load64(int(x,16)),load64(int(beta,16)),0x133f)
        yy=value(f'{native.publish(y):016x}')
        if not(aa<=1 or yy<1.1):omitted['stirling_mass_arm']+=1;continue
        cases.append(dict(x=x,alpha=alpha,beta=beta,bank_original=r))
    for alpha in (.125,.375,.5,.75,1.25,1.5,2.5,4.5,6.5):
        cases.append(dict(x=bits(.75),alpha=bits(alpha),beta=bits(1.),fixed_control=True))
    normalizers=sorted({bits(value(r['alpha'])+1.)for r in cases})
    frozen={**STAGE_BINDINGS,str(Path(bank).resolve()):sha(bank),str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-density-normalizer-questions/1',complete=True,
        records=cases,normalizer_arguments=normalizers,omitted=omitted,source_bindings=frozen,
        selection='All original density mismatches with nonwhole/low shapes below7 in the existing logarithmic-mass arm, then nine fixed ordinary controls. No output chooses a kernel.'))
    upper=service.ask('GAMMALN',[[client.number(x)]for x in normalizers])['rows']
    require_numeric=lambda rows:all(r.get('type')=='number'for r in rows)
    if not require_numeric(upper):raise ValueError('positive normalizers must be numeric')
    lookup=dict(zip(normalizers,upper))
    answers=service.ask('GAMMA.DIST',[[client.number(r['x']),client.number(r['alpha']),client.number(r['beta']),client.boolean(False)]for r in cases])['rows']
    records=[dict(r,normalizer_argument=bits(value(r['alpha'])+1.),
                  public_normalizer=lookup[bits(value(r['alpha'])+1.)],live_answer=a)for r,a in zip(cases,answers)]
    result=dict(format='dnacalc-gamma-density-normalizer-observations/1',complete=True,records=records,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Original bank failures remain unchanged. These independent current-build scalar readings '
               'and public normalizers are diagnostic observations, not a new bank capture or engine pass. '
               'The paired replacement is evaluated separately through the current registered baseline.')
    write(output,result);print(json.dumps(dict(rows=len(records),bank_cases=len(cases)-9,ordinary_controls=9)),flush=True)
    return result


def paired(checkout, output, cache, library, count=32, offline=False):
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    cw = 0x133f
    one = load64(int(bits(1), 16))
    exp_graph = primitive_stage_graph('EXP')

    def leaf(op, a, b):
        return native.leaf(op, a, b, cw)

    def store(a):
        return load64(native.publish(a))

    def scalar(a):
        return load64(int(bits(a), 16))

    def exp(a, stored=True):
        negative = bool(a[9] & 128)
        if negative:
            a = leaf('neg', a, a)
        result = native.evaluate(exp_graph, [], cw, extended=a)[:10]
        result = leaf('div', one, result) if negative else result
        return store(result) if stored else result

    # Adjacent windows and power-of-two edges; no current match status is read.
    shapes = (1.25, 1.5, 1.75, 2., 2.5, 3., 4., 6., 10., 19.)
    centers = (.125, .5, .7, 1., 1.1, 1.9, 2.5, 4., 8., 16.)
    cases = []
    for a in shapes:
        for center in centers:
            if center >= a:
                continue
            start = int(bits(center), 16)-count//2
            cases.extend((a, f'{start+k:016x}') for k in range(count))
    xs = sorted({x for _, x in cases})
    logs = dict(zip(xs, service.ask('LN', [[client.number(x)] for x in xs])['rows']))
    lg = dict(zip(shapes, service.ask('GAMMALN', [[client.number(a+1)] for a in shapes])['rows']))
    lg_density = dict(zip(shapes, service.ask('GAMMALN', [[client.number(a)] for a in shapes])['rows']))
    norm_gamma = dict(zip(shapes, service.ask('GAMMA', [[client.number(a+1)] for a in shapes])['rows']))
    questions = [[client.number(x), client.number(a), client.number(1.), client.boolean(flag)]
                 for a, x in cases for flag in (False, True)]
    answers = service.ask('GAMMA.DIST', questions)['rows']

    @lru_cache(maxsize=8192)
    def series(a, x, model):
        rr, term, summ = scalar(a), one, one
        terms = []
        if model in ('sum_zero', 'reverse'):
            summ = scalar(0.)
        xx = load64(int(x, 16))
        for iteration in range(1, 1000):
            rr = store(leaf('add', rr, one))
            if model == 'strict':
                t = value(f'{native.publish(term):016x}')
                z = value(f'{native.publish(rr):016x}')
                term = scalar(t*(value(x)/z))
                summ = scalar(value(f'{native.publish(summ):016x}') + value(f'{native.publish(term):016x}'))
            elif model in ('stored_quotient', 'resident_quotient', 'resident_sum', 'resident_term_sum', 'sum_zero', 'reverse'):
                quotient = leaf('div', xx, rr)
                if model in ('stored_quotient', 'sum_zero', 'reverse'):
                    quotient = store(quotient)
                term = leaf('mul', term, quotient)
                if model != 'resident_term_sum':
                    term = store(term)
                summ = leaf('add', summ, term)
                if model not in ('resident_sum', 'resident_term_sum'):
                    summ = store(summ)
                terms.append(term)
            elif model == 'multiply_then_divide':
                term = store(leaf('div', leaf('mul', term, xx), rr))
                summ = store(leaf('add', summ, term))
            else:
                raise ValueError('unknown isolated series model')
            # Fixed Cephes stopping criterion; no fitted iteration bound.
            test_sum = leaf('add', summ, one) if model in ('sum_zero', 'reverse') else summ
            if decode(term)/decode(test_sum) <= binary64(bits(2.**-53)):
                if model == 'reverse':
                    summ = scalar(0.)
                    for t in reversed(terms):
                        summ = store(leaf('add', summ, t))
                if model in ('sum_zero', 'reverse'):
                    summ = store(leaf('add', summ, one))
                return store(summ), iteration
        raise ValueError('series did not converge')

    totals, by_shape, records = Counter(), {}, []
    models = ('strict', 'stored_quotient', 'resident_quotient', 'multiply_then_divide',
              'resident_sum', 'resident_term_sum', 'sum_zero', 'reverse')
    for i, (a, x) in enumerate(cases):
        density, cdf = answers[2*i:2*i+2]
        if any(r['type'] != 'number' for r in (density, cdf, logs[x], lg[a])):
            raise ValueError('positive finite paired gamma call was not numeric')
        aa, xx = scalar(a), load64(int(x, 16))
        t = store(leaf('mul', aa, load64(int(logs[x]['bits'], 16))))
        argument = leaf('sub', leaf('sub', t, xx), load64(int(lg[a]['bits'], 16)))
        r = exp(argument)
        predicted_pdf = native.publish(leaf('div', store(leaf('mul', r, aa)), xx))
        density_consumers = dict(
            stored_exp_resident_product=native.publish(leaf('div', leaf('mul', r, aa), xx)),
            resident_exp_product=native.publish(leaf('div', leaf('mul', exp(argument, False), aa), xx)),
            ratio_first=native.publish(leaf('mul', r, store(leaf('div', aa, xx)))),
            stored_argument=native.publish(leaf('div', store(leaf('mul', exp(store(argument)), aa)), xx)))
        density_t = store(leaf('mul', scalar(a-1.), load64(int(logs[x]['bits'], 16))))
        density_argument = leaf('sub', leaf('sub', density_t, xx),
                                load64(int(lg_density[a]['bits'], 16)))
        direct_pdf = native.publish(exp(density_argument))
        divided_prefix = store(leaf('div', exp(leaf('sub', t, xx)),
                                    load64(int(norm_gamma[a]['bits'],16))))
        cell = inverse(Span(binary64(density['bits']), binary64(density['bits'])),
                       binary64(x), 'div', True, 'stored_pc64')
        # The intermediate r*a is stored before division by x.
        intermediate = BINARY64.coordinates(cell)
        r_domains = []
        for k in ([] if intermediate is None else range(intermediate[0], intermediate[1]+1)):
            y = BINARY64.value(k)
            domain = inverse(Span(y, y), binary64(bits(a)), 'mul', True, 'stored_pc64')
            coordinates = BINARY64.coordinates(domain)
            if coordinates:
                r_domains.append(coordinates)
        r_coordinates = sorted({k for lo, hi in r_domains for k in range(lo, hi+1)})
        if len(r_coordinates) > 8:
            raise ValueError('unexpectedly broad prefactor preimage')
        prefixes = [BINARY64.value(k) for k in r_coordinates]
        target = Span(binary64(cdf['bits']), binary64(cdf['bits']))
        sum_domains = [inverse(target, p, 'mul', True, 'stored_pc64') for p in prefixes]
        predictions, compatible, iterations = {}, {}, {}
        for model in models:
            summ, n = series(a, x, model)
            predictions[model] = f'{native.publish(leaf("mul", r, summ)):016x}'
            compatible[model] = any(domain.contains(decode(summ)) for domain in sum_domains)
            iterations[model] = n
        strict_sum, _ = series(a, x, 'strict')
        for label, prefactor in (('resident_exp', exp(argument, False)),
                                 ('stored_argument', exp(store(argument)))):
            predictions[label] = f'{native.publish(leaf("mul", prefactor, strict_sum)):016x}'
        predictions['divided_prefix'] = f'{native.publish(leaf("mul", divided_prefix, strict_sum)):016x}'
        row = dict(shape=bits(a), input=x, density=density, cdf=cdf,
                   public_LN=logs[x], public_GAMMALN=lg[a],
                   prefactor=f'{native.publish(r):016x}', pdf_prediction=f'{predicted_pdf:016x}',
                   direct_pdf=f'{direct_pdf:016x}', public_density_GAMMALN=lg_density[a],
                   divided_prefix=f'{native.publish(divided_prefix):016x}', public_GAMMA=norm_gamma[a],
                   density_consumers={k:f'{v:016x}' for k,v in density_consumers.items()},
                   permitted_prefactors=[bits(float(p)) for p in prefixes],
                   cdf_predictions=predictions, series_possible=compatible, iterations=iterations)
        tally = by_shape.setdefault(str(a), Counter())
        for counter in (totals, tally):
            counter['pairs'] += 1
            counter['prefactor_pdf_match'] += predicted_pdf == int(density['bits'], 16)
            counter['direct_pdf_match'] += direct_pdf == int(density['bits'], 16)
            for model, result in density_consumers.items():
                counter[model+'_pdf_match'] += result == int(density['bits'], 16)
            counter['empty_prefactor_preimage'] += not prefixes
            for model in models:
                counter[model+'_possible'] += compatible[model]
                counter[model+'_cdf_match'] += predictions[model] == cdf['bits']
            for model in ('resident_exp', 'stored_argument', 'divided_prefix'):
                counter[model+'_cdf_match'] += predictions[model] == cdf['bits']
        records.append(row)
    result = dict(format='dnacalc-gamma-paired-stages/1', complete=True,
        selection='Fixed shape grid; adjacent windows at dyadic edges and ordinary decimal sites, restricted to x<a<20.',
        counts=totals, by_shape=by_shape, records=records, receipts=service.receipts,
        source_bindings={str(p):sha(p) for p in (Path(__file__), Path(library))},
        limits='The density inverse conditions on stored (r*a)/x and a shared prefactor. Empty constraints refute that family; nonempty constraints do not identify a series. Public GAMMALN and LN are independently observed, never fitted. Four ranked textbook series graphs differ at the recurrence only; no whole-program enumeration or count-selected changes.')
    write(output, result)
    print(json.dumps(dict(counts=totals, by_shape=by_shape), indent=2), flush=True)
    return result


def power_prefactor(checkout, output, cache, library, offline=False):
    """Remove the density's final scaling exactly, then invert one EXP input.

    Shape and x are powers of two, so every proposed a/x consumer is exact.
    The resulting r is an observed binary64, not a fitted prefactor. A common
    c in EXP(stored(a*LN(x))-x-c) must satisfy every x at the same shape.
    """
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    shapes = (2., 4., 8., 16.)
    cases = [(a, 2.**e) for a in shapes for e in range(-12, int(math.log2(a)))]
    questions = [[client.number(x), client.number(a), client.number(1.), client.boolean(c)]
                 for a, x in cases for c in (False, True)]
    answers = service.ask('GAMMA.DIST', questions)['rows']
    xs = sorted({x for _, x in cases})
    logs = dict(zip(xs, service.ask('LN', [[client.number(x)] for x in xs])['rows']))
    lg = dict(zip(shapes, service.ask('GAMMALN', [[client.number(a+1)] for a in shapes])['rows']))
    g = Builder(1)
    t = g.op('mul', 0, g.op('fldl2e', store=False), store=False)
    whole = g.op('frndint', t, store=False)
    small = g.op('f2xm1', g.op('sub', t, whole, store=False), store=False)
    positive = g.op('fscale', g.op('add', g.one, small, store=False), whole, store=False)
    out = g.op('div', g.one, positive, store=False)
    graph = g.graph(out, out)
    records, groups, totals = [], {}, Counter()
    log_groups = {}
    ln2 = native.leaf('fldln2', bytes(10), bytes(10), 0x133f)
    for i, (a, x) in enumerate(cases):
        density, cdf = answers[2*i:2*i+2]
        if density['type'] != 'number' or cdf['type'] != 'number':
            raise ValueError('finite power-of-two density/CDF must be numeric')
        # Exact rescaling, independently checked as a rational.
        prefix = value(density['bits'])*x/a
        if binary64(bits(prefix)) != binary64(density['bits'])*binary64(bits(x))/binary64(bits(a)):
            raise ValueError('density prefactor extraction rounded')
        argument = positive_extended_cell(native, graph, bits(prefix),
                                          binary64(bits(-math.log(prefix))), decreasing=True)
        if argument is None:
            raise ValueError('empty resident EXP cell requires a separate proof')
        m = Span(*(decode(bytes.fromhex(argument[k])) for k in ('low', 'high')))
        t = native.leaf('mul', load64(int(bits(a), 16)), load64(int(logs[x]['bits'], 16)), 0x133f)
        t = load64(native.publish(t))
        base = decode(native.leaf('sub', t, load64(int(bits(x), 16)), 0x133f))
        domain = EXTENDED64.pullback(m).linear(1, base)
        group = groups.setdefault(str(a), dict(domain=None, contradicted=False, rows=[]))
        group['rows'].append(i)
        if not group['contradicted']:
            group['domain'] = domain if group['domain'] is None else group['domain'].meet(domain)
            if group['domain'] is None:
                group['contradicted'] = True
                group['first_contradiction'] = i
        c = binary64(lg[a]['bits'])
        base_domain = EXTENDED64.pullback(m).linear(-1, c)
        t_domain = EXTENDED64.pullback(base_domain).linear(1, binary64(bits(x)))
        log_domain = inverse(t_domain, binary64(bits(a)), 'mul', True, 'stored_pc64')
        log_group = log_groups.setdefault(bits(x), dict(domain=None, contradicted=False, rows=[]))
        log_group['rows'].append(i)
        if not log_group['contradicted']:
            log_group['domain'] = log_domain if log_group['domain'] is None else log_group['domain'].meet(log_domain)
            if log_group['domain'] is None:
                log_group['contradicted'] = True
                log_group['first_contradiction'] = i
        resident_log = native.leaf('fyl2x', ln2, load64(int(bits(x),16)), 0x133f)
        resident_t = native.leaf('mul',load64(int(bits(a),16)),resident_log,0x133f)
        resident_base = decode(native.leaf('sub',resident_t,load64(int(bits(x),16)),0x133f))
        resident_domain = EXTENDED64.pullback(m).linear(1,resident_base)
        rr, term, summ = a, 1., 1.
        while True:
            rr += 1.
            term *= x/rr
            summ += term
            if term/summ <= 2.**-53:
                break
        published = native.publish(native.leaf('mul', load64(int(bits(prefix), 16)),
                                              load64(int(bits(summ), 16)), 0x133f))
        row = dict(shape=bits(a), input=bits(x), density=density, cdf=cdf,
                   extracted_prefactor=bits(prefix), public_LN=logs[x], public_GAMMALN=lg[a],
                   exponential_magnitude_cell=argument, normalizer_domain=domain.json(),
                   public_normalizer_possible=domain.contains(binary64(lg[a]['bits'])),
                   public_normalizer_with_resident_log_product=resident_domain.contains(c),
                   resident_log_product_normalizer_domain=resident_domain.json(),
                   log_domain=log_domain.json() if log_domain else None,
                   cdf_with_extracted_prefactor=f'{published:016x}')
        totals['rows'] += 1
        totals['public_normalizer_possible'] += row['public_normalizer_possible']
        totals['public_normalizer_with_resident_log_product'] += row['public_normalizer_with_resident_log_product']
        totals['cdf_with_extracted_prefactor_match'] += published == int(cdf['bits'], 16)
        records.append(row)
    for group in list(groups.values())+list(log_groups.values()):
        domain = group.pop('domain')
        coords = BINARY64.coordinates(domain) if domain is not None else None
        group.update(domain=domain.json() if domain else None,
                     permitted_binary64_count=0 if coords is None else coords[1]-coords[0]+1,
                     binary64_values=[] if coords is None else [bits(float(BINARY64.value(k)))
                         for k in range(coords[0], min(coords[1]+1, coords[0]+8))])
    result = dict(format='dnacalc-gamma-power-prefactor/1', complete=True, counts=totals,
        groups=groups, log_groups=log_groups, records=records, graph=graph, receipts=service.receipts,
        source_bindings={str(p):sha(p) for p in (Path(__file__), Path(library))},
        limits='Power-of-two density scaling is exact. The extracted r is a conditional observed prefactor; sharing it with the CDF is separately tested. Normalizer domains condition on one stored public LN/product and the declared resident negative EXP chain, with PC64 subtraction. A nonempty common interval is a necessary constraint, never an identified coefficient or a fitted implementation.')
    write(output, result)
    print(json.dumps(dict(counts=totals, groups=groups), indent=2), flush=True)
    return result


def power_primitives(study, output, cache, library, direct_density=False):
    """Invert one exponential primitive with the density scale removed exactly.

    The fixed textbook exponentials vary negative handling and constant
    precision. The unknown normalizer must be common to all x at a
    fixed shape. No constant inferred here is installed in an implementation.
    """
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete power-prefactor observations')
    answers, bindings, ordered = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    native = Native(library)
    graphs = {}
    # Code-shape rank is fixed before the outputs: the ordinary signed
    # textbook sequence, then the public-EXP reciprocal sequence.
    for kind in ('signed', 'reciprocal', 'signed_constant53'):
        g = Builder(1)
        exponent = 0 if kind == 'reciprocal' else g.op('neg', 0, store=False)
        constant = (g.constant('3ff71547652b82fe') if kind == 'signed_constant53'
                    else g.op('fldl2e', store=False))
        t = g.op('mul', exponent, constant, store=False)
        whole = g.op('frndint', t, store=False)
        fraction = g.op('sub', t, whole, store=False)
        small = g.op('f2xm1', fraction, store=False)
        result = g.op('fscale', g.op('add', g.one, small, store=False), whole, store=False)
        if kind == 'reciprocal':
            result = g.op('div', g.one, result, store=False)
        graphs[kind] = g.graph(result, result)
    groups, records = {}, []
    for row in doc['records']:
        a, x = map(binary64, (row['shape'], row['input']))
        density = row['density']
        original_args = [dict(type='number', bits=b) for b in
                         (row['input'], row['shape'], bits(1.))] + [dict(type='boolean', value=False, mode='cell')]
        if answers[('GAMMA.DIST', json.dumps(original_args, sort_keys=True))] != density:
            raise ValueError('power density observation changed')
        target = density['bits'] if direct_density else row['extracted_prefactor']
        prefix = binary64(target)
        if not direct_density and prefix != binary64(density['bits'])*x/a:
            raise ValueError('power-of-two density scaling was not exact')
        coefficient=bits(float(a-1)) if direct_density else row['shape']
        t = native.leaf('mul', load64(int(coefficient, 16)),
                        load64(int(row['public_LN']['bits'], 16)), 0x133f)
        t = load64(native.publish(t))
        base = decode(native.leaf('sub', t, load64(int(row['input'], 16)), 0x133f))
        predictions = {}
        for kind, graph in graphs.items():
            cell = positive_extended_cell(native, graph, target,
                                          binary64(bits(-math.log(float(prefix)))), decreasing=True)
            group = groups.setdefault((kind, row['shape']), dict(domain=None, contradicted=False,
                                                                 rows=0, public_possible=0))
            group['rows'] += 1
            if cell is None:
                domain = None
            else:
                magnitudes = Span(*(decode(bytes.fromhex(cell[k])) for k in ('low', 'high')))
                domain = EXTENDED64.pullback(magnitudes).linear(1, base)
            if not group['contradicted']:
                group['domain'] = domain if group['domain'] is None else (
                    None if domain is None else group['domain'].meet(domain))
                if group['domain'] is None:
                    group['contradicted'] = True
                    group['first_contradiction'] = len(records)
            c = binary64(row['public_GAMMALN']['bits'])
            possible = None if direct_density else domain is not None and domain.contains(c)
            group['public_possible'] += possible is True
            predictions[kind] = dict(exponential_cell=cell,
                                     normalizer_domain=domain.json() if domain else None,
                                     public_normalizer_possible=possible)
        records.append(dict(shape=row['shape'], input=row['input'],
                            density=density, exponential_target=target,
                            extracted_prefactor=row['extracted_prefactor'], primitives=predictions))
    reported = []
    for (kind, shape), group in groups.items():
        domain = group.pop('domain')
        coords = BINARY64.coordinates(domain) if domain else None
        reported.append(dict(primitive=kind, shape=shape, **group,
                             common_normalizer=domain.json() if domain else None,
                             binary64_count=0 if coords is None else coords[1]-coords[0]+1,
                             binary64_values=[] if coords is None else [bits(float(BINARY64.value(k)))
                                 for k in range(coords[0], min(coords[1]+1, coords[0]+8))]))
    result = dict(format=('dnacalc-gamma-density-primitive-inverses/1' if direct_density else
                          'dnacalc-gamma-power-primitive-inverses/1'), complete=True,
                  exponent_coefficient='a-1' if direct_density else 'a',
                  ranked_primitives=list(graphs), graphs=graphs, groups=reported, records=records,
                  original_calls=len(ordered), receipts=doc['receipts'],
                  source_bindings=frozen,
                  original_bindings=bindings,
                  limits=('The observed density directly targets EXP(stored((a-1)*public_LN(x))-x-C(a)). '
                          'The shape-a public GAMMALN is not observed here; no normalizer value is compared. '
                          if direct_density else 'Original power-of-two density scaling removes the consumer exactly. ')+
                  'Each inverse tests one fixed native exponential with a shared but unknown normalizer. '
                  'Empty common constraints refute that primitive/producer composition; nonempty '
                  'constraints do not identify a coefficient or a private logarithm. No kernel is '
                  'chosen from match counts and no inferred constant is installed.')
    write(output, result)
    print(json.dumps(reported, indent=2), flush=True)
    return result


def shared_log_constraints(study, output, cache, library, scaled_log=False):
    """Test a shared unknown log and normalizer without fitting either.

    Powers of two make a*L exact for a binary64 L. For this input grid,
    a*L-x is also exact in PC64 whenever .5 <= abs(L) <= 16 (L(1)=0).
    Therefore each inverse gives one necessary difference constraint on
    L(z)-C(a)/a, with z=x or the explicitly selected exact x/a. A negative cycle refutes the whole declared family even
    after relaxing all interval endpoints to closed bounds.
    """
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete primitive inverses')
    _, bindings, _ = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    from fractions import Fraction as Q
    results = []
    for primitive in doc['ranked_primitives']:
        edges = []
        nodes = set()
        for index, row in enumerate(doc['records']):
            a, x = binary64(row['shape']), binary64(row['input'])
            cell = row['primitives'][primitive]['exponential_cell']
            if cell is None:
                raise ValueError('empty primitive cell needs its own certificate')
            magnitudes = Span(*(decode(bytes.fromhex(cell[k])) for k in ('low', 'high')))
            domain = EXTENDED64.pullback(magnitudes).linear(-1/a, x/a)
            log_input=bits(float(x/a)) if scaled_log else row['input']
            if scaled_log and binary64(log_input)!=x/a:raise ValueError('scaled log argument rounded')
            log, normalizer = 'L:'+log_input, 'C/a:'+row['shape']
            nodes.update((log, normalizer))
            edges.extend([dict(source=normalizer, target=log, bound=domain.high, row=index),
                          dict(source=log, target=normalizer, bound=-domain.low, row=index)])
        distance = {name:Q(0) for name in nodes}
        predecessor = {}
        changed = None
        for _ in range(len(nodes)):
            changed = None
            for i, edge in enumerate(edges):
                candidate = distance[edge['source']]+edge['bound']
                if candidate < distance[edge['target']]:
                    distance[edge['target']] = candidate
                    predecessor[edge['target']] = i
                    changed = edge['target']
            if changed is None:
                break
        cycle = []
        if changed is not None:
            node = changed
            for _ in nodes:
                node = edges[predecessor[node]]['source']
            start = node
            while True:
                edge = edges[predecessor[node]]
                cycle.append(edge)
                node = edge['source']
                if node == start:
                    break
                if len(cycle) > len(nodes):
                    raise ValueError('invalid negative cycle reconstruction')
            total = sum((e['bound'] for e in cycle), Q(0))
            if total >= 0:
                raise ValueError('cycle does not prove a contradiction')
        else:
            total = Q(0)
            if any(distance[e['target']] > distance[e['source']]+e['bound'] for e in edges):
                raise ValueError('feasible witness violates a necessary constraint')
        serial = lambda q:[str(q.numerator), str(q.denominator)]
        results.append(dict(primitive=primitive, continuous_constraints_feasible=changed is None,
                            negative_cycle=[dict(e, bound=serial(e['bound'])) for e in cycle],
                            cycle_bound_sum=serial(total),
                            edges=[dict(e, bound=serial(e['bound'])) for e in edges],
                            feasible_potentials=None if changed is not None else
                                {k:serial(v) for k,v in distance.items()}))
    result = dict(format='dnacalc-gamma-shared-log-constraints/1', complete=True,
                  logarithm_argument='x/a' if scaled_log else 'x',
                  results=results, receipts=doc['receipts'], source_bindings=frozen,
                  limits='Necessary constraints only: one stored binary64 private L at the stated argument, '
                  '.5<=abs(L)<=16 away from x=1, L(1)=0, one normalizer C(a), '
                  'stored a*L, resident differences and one declared exponential, with '
                  'the exact a/x density consumer. Continuous variables and closed '
                  'endpoints relax the model. A negative cycle is a refutation even '
                  'of that relaxation. A feasible witness identifies no private value '
                  'and supplies no implementation or count-selected correction.')
    write(output, result)
    print(json.dumps([dict(primitive=r['primitive'], feasible=r['continuous_constraints_feasible'],
                           cycle_rows=[e['row'] for e in r['negative_cycle']],
                           cycle_bound_sum=r['cycle_bound_sum']) for r in results]), flush=True)
    return result


def split_prefactor(checkout, study, output, cache, library, offline=False):
    """Test the literal power-times-exponential normalizer on dyadic inputs."""
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete power-prefactor study')
    _, bindings, _ = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    xs = sorted({r['input'] for r in doc['records']})
    shapes = sorted({r['shape'] for r in doc['records']})
    # Independent public leaves; original density/CDF answers are reused.
    exps = dict(zip(xs, service.ask('EXP', [[client.number(bits(-value(x)))] for x in xs])['rows']))
    gammas = dict(zip(shapes, service.ask('GAMMA', [[client.number(value(a)+1.)] for a in shapes])['rows']))
    groups = {}
    records = []
    for row in doc['records']:
        a, x = binary64(row['shape']), binary64(row['input'])
        if a.denominator != 1:
            raise ValueError('power experiment needs an integer shape')
        e, gamma = exps[row['input']], gammas[row['shape']]
        if e['type'] != 'number' or gamma['type'] != 'number':
            raise ValueError('positive finite leaf must be numeric')
        product = x**a.numerator*binary64(e['bits'])
        # On this grid the power and its multiplication by EXP(-x) are exact.
        if BINARY64.round(product) != product:
            raise ValueError('dyadic numerator scaling rounded')
        target = binary64(row['extracted_prefactor'])
        domain = inverse(Span(target, target), product, 'div', False, 'stored_pc64')
        group = groups.setdefault(row['shape'], dict(domain=None, contradicted=False, rows=0,
                                                     public_gamma_possible=0))
        group['rows'] += 1
        if not group['contradicted']:
            group['domain'] = domain if group['domain'] is None else (
                None if domain is None else group['domain'].meet(domain))
            if group['domain'] is None:
                group['contradicted'] = True
                group['first_contradiction'] = len(records)
        possible = domain is not None and domain.contains(binary64(gamma['bits']))
        group['public_gamma_possible'] += possible
        records.append(dict(shape=row['shape'], input=row['input'], public_exp=e,
                            public_gamma=gamma, density=row['density'],
                            extracted_prefactor=row['extracted_prefactor'],
                            normalizer_domain=domain.json() if domain else None,
                            public_gamma_possible=possible))
    reported = []
    for shape, group in groups.items():
        domain = group.pop('domain')
        coords = BINARY64.coordinates(domain) if domain else None
        reported.append(dict(shape=shape, **group, common_normalizer=domain.json() if domain else None,
                             binary64_count=0 if coords is None else coords[1]-coords[0]+1,
                             binary64_values=[] if coords is None else [bits(float(BINARY64.value(k)))
                                 for k in range(coords[0], min(coords[1]+1, coords[0]+8))]))
    result = dict(format='dnacalc-gamma-split-prefactor/1', complete=True,
                  groups=reported, records=records, receipts=doc['receipts']+service.receipts,
                  source_bindings=frozen,
                  limits='Exact powers of two remove the power and numerator-product rounding. '
                  'Original density outputs and new independent public EXP/GAMMA leaves constrain '
                  'one unknown shared Gamma(a+1) divisor. A common interval is necessary only; '
                  'no normalizer or algorithm is chosen from match counts or installed.')
    write(output, result)
    print(json.dumps(reported, indent=2), flush=True)
    return result


def poisson_pairs(checkout, study, output, cache, library, offline=False):
    """Check a possible public Poisson dependency before changing gamma helpers."""
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete gamma observations')
    _, bindings, _ = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    questions = [[client.number(value(r['shape'])-1.), client.number(r['input']), client.boolean(flag)]
                 for r in doc['records'] for flag in (False, True)]
    answers = service.ask('POISSON.DIST', questions)['rows']
    counts, records = Counter(), []
    for index, row in enumerate(doc['records']):
        pmf, cdf = answers[2*index:2*index+2]
        same = pmf == row['density']
        complement = cdf['type'] == 'number' and row['cdf'] == dict(type='number', bits=bits(1.-value(cdf['bits'])))
        counts['pairs'] += 1
        counts['pdf_equals_public_poisson'] += same
        counts['cdf_equals_public_poisson_complement'] += complement
        records.append(dict(shape=row['shape'], input=row['input'], gamma_density=row['density'],
                            gamma_cdf=row['cdf'], poisson_pmf=pmf, poisson_cdf=cdf,
                            density_equal=same, cdf_complement_equal=complement))
    result = dict(format='dnacalc-gamma-poisson-public-pairs/1', complete=True, counts=counts,
                  records=records, receipts=doc['receipts']+service.receipts, source_bindings=frozen,
                  limits='Original gamma answers, independently measured Poisson calls at integer '
                  'shape minus one. An equality is a functional observation, not proof of a '
                  'private call or of the noninteger route. Complement roundoff is explicit; '
                  'no correction, selector or implementation is inferred from these counts.')
    write(output, result)
    print(json.dumps(counts), flush=True)
    return result


def gammaln_store_edges(checkout, output, cache, library, offline=False, count=1024):
    """Observe the mean-error jump at each final polynomial's binade.

    A final binary64 sum rounds the constant coefficient on a coarser
    lattice when the polynomial doubles. An unstored final PC64 sum has
    eleven more bits there. Edges come from the existing coefficient table,
    before any new Excel answer. Coefficients are never changed in this study.
    """
    import re
    import struct
    try:
        import mpmath as mp
    except ImportError:
        sys.path.insert(0, '/home/ubuntu/claude-work/gamma/pylib')
        import mpmath as mp
    mp.mp.prec = 200
    path = Path(checkout)/'crates/calc-func/src/ported/excel_numeric/gammaln.rs'
    frozen = {**STAGE_BINDINGS, str(path.resolve()):sha(path),
              str(Path(library).resolve()):sha(library), str(Path(mp.__file__).resolve()):sha(mp.__file__)}
    text = path.read_text()
    tables = {}
    for name in ('P', 'Q'):
        block = text.split('const B3_'+name+':', 1)[1].split('];', 1)[0]
        tables[name] = [v.replace('_', '').lower() for v in re.findall(r'0x([0-9A-Fa-f_]+)', block)]
        if len(tables[name]) != 9:
            raise ValueError('changed B3 coefficient table shape')
    def polynomial(name, x):
        result = mp.mpf(0)
        for b in reversed(tables[name]):
            result = result*x+mp.mpf(value(b))
        return result
    edges = []
    for name, level in (('P', 4096), ('P', 8192), ('Q', 8192), ('Q', 16384)):
        low, high = mp.mpf(2), mp.mpf(3)
        for _ in range(180):
            middle = (low+high)/2
            if polynomial(name, middle) < level:
                low = middle
            else:
                high = middle
        center = float((low+high)/2)
        windows = []
        for delta in (2.**-16, 2.**-20):
            for sign in (-1, 1):
                first = int(bits(center+sign*delta), 16)-count//2
                windows.append(dict(center=bits(center+sign*delta), delta=bits(delta), sign=sign,
                                    inputs=[f'{first+k:016x}' for k in range(count)]))
        edges.append(dict(polynomial=name, level=level, estimated_edge=bits(center), windows=windows))
    graphs = {}
    for numerator_stored, denominator_stored in ((False,False), (True,False), (False,True)):
        g = Builder(1)
        def op(kind, a, b, stored=True):
            index = g.op(kind, a, b, store=False)
            if stored:
                g.nodes[index]['site'] = 0
                g.sites = 1
            return index
        def chain(name, stored):
            v = g.constant(tables[name][-1])
            for i in range(7, -1, -1):
                v = op('mul', v, 0)
                v = op('add', v, g.constant(tables[name][i]), stored=i != 0 or stored)
            return v
        p, q = chain('P', numerator_stored), chain('Q', denominator_stored)
        ratio = op('div', p, q)
        difference = op('sub', 0, g.constant(bits(2.)), stored=False)
        result = op('mul', difference, ratio)
        label = 'P'+('stored' if numerator_stored else 'resident')+'_Q'+('stored' if denominator_stored else 'resident')
        graphs[label] = g.graph(result, result)
    manifest = dict(format='dnacalc-gammaln-store-edge-questions/1', complete=True,
                    tables=tables, edges=edges, graphs=graphs, source_bindings=frozen,
                    selection='Existing polynomial crossings of four fixed powers of two; '
                    'two fixed symmetric distances and consecutive-double windows. No Excel output selects a site.')
    write(str(output)+'.questions.json', manifest)
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    records = []
    for edge in edges:
        for window in edge['windows']:
            xs = window['inputs']
            answers = service.ask('GAMMALN', [[client.number(x)] for x in xs])['rows']
            center = mp.mpf(value(window['center']))
            f0, f1, f2 = mp.loggamma(center), mp.digamma(center), mp.polygamma(1, center)
            sums = dict(Excel=mp.mpf(0), **{k:mp.mpf(0) for k in graphs})
            squares = dict(Excel=mp.mpf(0), **{k:mp.mpf(0) for k in graphs})
            for x, answer in zip(xs, answers):
                if answer['type'] != 'number':
                    raise ValueError('GAMMALN edge answer must be numeric')
                xx = mp.mpf(value(x)); dx = xx-center
                truth = f0+f1*dx+f2*dx*dx/2
                unit = mp.mpf(math.ulp(value(answer['bits'])))
                predictions = {k:native.publish(native.evaluate(graph, [int(x,16)], 0x133f, mask=1)[:10])
                               for k,graph in graphs.items()}
                observed = dict(Excel=int(answer['bits'],16), **predictions)
                for name, raw in observed.items():
                    residual = (mp.mpf(struct.unpack('>d',struct.pack('>Q',raw))[0])-truth)/unit
                    sums[name] += residual
                    squares[name] += residual*residual
            n = len(xs)
            means = {k:float(v/n) for k,v in sums.items()}
            sem = {k:float(mp.sqrt(max(mp.mpf(0), squares[k]/n-(sums[k]/n)**2)/n)) for k in sums}
            records.append(dict(polynomial=edge['polynomial'], level=edge['level'],
                                center=window['center'], delta=window['delta'], sign=window['sign'],
                                first_input=xs[0], last_input=xs[-1], calls=n,
                                mean_error_in_result_ulps=means, descriptive_standard_error=sem))
    result = dict(format='dnacalc-gammaln-final-sum-store-edges/1', complete=True,
                  records=records, receipts=service.receipts, source_bindings=frozen,
                  questions_sha256=sha(str(output)+'.questions.json'),
                  limits='Stage-binade observations only. Means use a 200-bit mathematical '
                  'log-gamma and a quadratic expansion in each consecutive window; '
                  'the third-derivative remainder is below 1e-20 result ULP on this grid. '
                  'Standard errors describe within-window residual variation, not independent '
                  'random sampling. Existing coefficients remain approximate and unchanged; '
                  'no table, store or threshold is chosen by exact-match counts.')
    write(output, result)
    print(json.dumps(records, indent=2), flush=True)
    return result


def unit_shape_exponential(checkout,bank,manifest,output,cache,library,offline=False):
    """Use shape-one density to expose an independent exponential wrapper.

    Existing original density inputs choose centers only. Their answers are
    retained, never re-asked; both successor inputs get paired new calls.
    """
    import gzip
    source=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz'
    old=[]
    for line in gzip.open(source,'rt'):
        raw=json.loads(line);args=raw['args']
        if len(args)!=1 or len(args[0])!=4 or any(a['kind']!='number' for a in args[0]):continue
        x,a,b,c=(v['value'] for v in args[0])
        if a==b==1. and c==0. and 0.<x<700.:old.append(raw)
    if not old:raise ValueError('no original unit-shape densities')
    primitive=Path(manifest);doc=json.loads(primitive.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete exponential graphs')
    graphs=doc['graphs'];known={r['args'][0][0]['bits_hex'].removeprefix('0x') for r in old}
    inputs=sorted({f'{int(x,16)+offset:016x}' for x in known for offset in (1,2)}-known)
    if len(inputs)>1500:raise ValueError('bounded adjacent-input study only')
    native=Native(library);predictions={x:{name:f'{native.publish(native.evaluate(graph,[int(x,16)],0x133f)[:10]):016x}'
                                        for name,graph in graphs.items()} for x in inputs}
    frozen={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(primitive.resolve()):sha(primitive),
            str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-unit-shape-questions/1',complete=True,
        original_centers=old,inputs=inputs,predictions=predictions,graphs=graphs,source_bindings=frozen,
        selection='First two successors of every original positive finite unit-scale shape-one density below700. All original density arguments are excluded before new calls; no current match result selects inputs.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    densities=service.ask('GAMMA.DIST',[[client.number(x),client.number(1.),client.number(1.),client.number(0.)] for x in inputs])['rows']
    exps=service.ask('EXP',[[client.number(bits(-value(x)))] for x in inputs])['rows']
    counts=Counter();records=[]
    for x,d,e in zip(inputs,densities,exps):
        if d['type']!='number' or e['type']!='number':raise ValueError('nonnumeric bounded density/EXP')
        counts['paired_rows']+=1;counts['same_public_EXP']+=d==e
        for name,prediction in predictions[x].items():
            counts['density/'+name]+=d==dict(type='number',bits=prediction)
            counts['EXP/'+name]+=e==dict(type='number',bits=prediction)
        records.append(dict(input=x,density=d,public_EXP=e,predictions=predictions[x]))
    result=dict(format='dnacalc-gamma-unit-shape-exponential/1',complete=True,counts=counts,records=records,
        original_centers=old,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Paired original current-build calls at new adjacent inputs. Shape1/unit scale '
        'mathematically isolates exp(-x); this does not assume other shapes call the same '
        'private exponential or identify a resident-argument ABI. The three fixed native '
        'graphs are comparisons, never selected by their number of matches. Historical '
        'center answers retain their separate environment and are not new coverage.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def unit_shape_prefactor(checkout,study,output,cache,library,offline=False):
    """Test the literal log/prefactor/divide route without cancelling x."""
    from geomean import verified_study_answers
    import itertools
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete unit-shape pairs')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    questions=Path(str(source)+'.questions.json')
    if sha(questions)!=doc['questions_sha256']:raise ValueError('changed original question manifest')
    graphs=json.loads(questions.read_text())['graphs'];inputs=[r['input'] for r in doc['records']]
    profiles=[dict(logarithm=log,argument_stored=arg,prefix_stored=prefix,division=division,exponential=exp)
              for log,arg,prefix,division,exp in itertools.product(
                  ('public_LN','resident_LN'),(True,False),(True,False),('stored_pc64','pc53'),graphs)]
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-unit-prefactor-questions/1',complete=True,
        inputs=inputs,profiles=profiles,source_bindings=STAGE_BINDINGS,
        selection='EXP(LN(x)-x)/x, three fixed exponentials and one-at-a-time publications. '
        'No Gamma1=exp(-x) algebraic cancellation, inferred coefficient or fitted normalizer.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    logs=service.ask('LN',[[client.number(x)] for x in inputs])['rows']
    native=Native(library);ln2=native.leaf('fldln2',bytes(10),bytes(10),0x133f)
    counts=[Counter() for _ in profiles];records=[]
    for row,log in zip(doc['records'],logs):
        if log['type']!='number':raise ValueError('nonnumeric public logarithm')
        xx=load64(int(row['input'],16));ll=load64(int(log['bits'],16))
        resident=native.leaf('fyl2x',ln2,xx,0x133f);predictions=[]
        for index,p in enumerate(profiles):
            magnitude=native.leaf('sub',xx,ll if p['logarithm']=='public_LN' else resident,0x133f)
            if p['argument_stored']:magnitude=load64(native.publish(magnitude))
            prefix=native.evaluate(graphs[p['exponential']],[],0x133f,extended=magnitude)[:10]
            if p['prefix_stored']:prefix=load64(native.publish(prefix))
            raw=native.leaf('div',prefix,xx,0x133f if p['division']=='stored_pc64' else 0x123f)
            predicted=f'{native.publish(raw):016x}';predictions.append(predicted)
            counts[index]['matched' if row['density']==dict(type='number',bits=predicted) else 'mismatched']+=1
        records.append(dict(input=row['input'],density=row['density'],public_EXP=row['public_EXP'],public_LN=log,predictions=predictions))
    frozen={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source),str(questions.resolve()):sha(questions),
            str(Path(library).resolve()):sha(library)}
    result=dict(format='dnacalc-gamma-unit-prefactor-stores/1',complete=True,profiles=profiles,counts=counts,
        records=records,receipts=doc['receipts']+service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='All48 literal profiles are fixed before new independent LN leaves. Original '
        'density and EXP answers remain unchanged. A full agreement establishes only this '
        'sampled functional composition, not uniqueness, other shapes or an unmeasured '
        'normalizer. A failed profile is never replaced by its nearest answer or chosen by counts.')
    write(output,result);print(json.dumps(dict(profiles=len(profiles),complete_survivors=[(i,p) for i,(p,c) in enumerate(zip(profiles,counts)) if not c['mismatched']])));return result


def unit_shape_constants(study,output,cache,library):
    """Constrain one common normalization or exponent shift, without installing it."""
    from geomean import verified_study_answers
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete shape-one pairs')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    questions=Path(str(source)+'.questions.json')
    if sha(questions)!=doc['questions_sha256']:raise ValueError('changed original primitive questions')
    graphs=json.loads(questions.read_text())['graphs'];native=Native(library)
    groups={name:dict(domain=None,contradicted=False,rows=0) for name in
            ('divide_stored_pc64','multiply_stored_pc64',*(('shift/'+k) for k in graphs))}
    records=[]
    for index,row in enumerate(doc['records']):
        target=binary64(row['density']['bits']);e=binary64(row['public_EXP']['bits']);x=binary64(row['input'])
        domains={name:inverse(Span(target,target),e,op,left,'stored_pc64') for name,op,left in
                 [('divide_stored_pc64','div',False),('multiply_stored_pc64','mul',True)]}
        cells={}
        for kind,graph in graphs.items():
            cell=positive_extended_cell(native,graph,row['density']['bits'],binary64(bits(-math.log(float(target)))),decreasing=True,max_radius_bits=56)
            cells[kind]=cell
            if cell is None:domain=None
            else:
                magnitudes=Span(*(decode(bytes.fromhex(cell[k])) for k in ('low','high')))
                domain=EXTENDED64.pullback(magnitudes).linear(1,-x)
            domains['shift/'+kind]=domain
        for name,domain in domains.items():
            group=groups[name];group['rows']+=1
            if not group['contradicted']:
                group['domain']=domain if group['domain'] is None else (None if domain is None else group['domain'].meet(domain))
                if group['domain'] is None:group['contradicted']=True;group['first_contradiction']=index
        records.append(dict(input=row['input'],density=row['density'],public_EXP=row['public_EXP'],
            domains={k:v.json() if v is not None else None for k,v in domains.items()},exponential_cells=cells))
    reported={}
    for name,group in groups.items():
        domain=group.pop('domain');coords=BINARY64.coordinates(domain) if domain is not None else None
        reported[name]=dict(group,common_interval=domain.json() if domain is not None else None,
            binary64_count=0 if coords is None else coords[1]-coords[0]+1)
    result=dict(format='dnacalc-gamma-unit-shape-normalizer-constraints/1',complete=True,groups=reported,records=records,
        source_bindings={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source),str(questions.resolve()):sha(questions),str(Path(library).resolve()):sha(library)},
        receipts=doc['receipts'],
        limits='One common multiplicative/divisive normalization of the measured public EXP '
        'or one common additive shift before a declared native exponential. Intersections '
        'are necessary only; no constant is chosen or installed. An empty cell rules out '
        'that complete declared scalar normalization, not a log/power route or a private '
        'exponential. Shape-one density is mathematically exp(-x), but its implementation '
        'is not assumed to be an exponential-only wrapper.')
    write(output,result);print(json.dumps(reported),flush=True);return result


def unit_shape_consumers(study,output,cache):
    """Invert only the final density consumer, with an arbitrary prefix."""
    from geomean import verified_study_answers
    from exact import operation
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete original unit-shape study')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    profiles=[dict(kind=kind,model=model) for kind in ('divide_by_x','multiply_reciprocal_x')
              for model in ('strict','stored_pc64')]
    counts=[Counter() for _ in profiles];records=[]
    for row in doc['records']:
        x=binary64(row['input']);y=binary64(row['density']['bits']);constraints=[]
        for index,p in enumerate(profiles):
            if p['kind']=='divide_by_x':constant=x;op='div'
            else:constant=operation(1,x,'div',p['model']);op='mul'
            domain=inverse(Span(y,y),constant,op,True,p['model'])
            coords=BINARY64.coordinates(domain) if domain is not None else None
            resident=EXTENDED64.coordinates(domain) if domain is not None else None
            n=0 if coords is None else coords[1]-coords[0]+1
            tally=counts[index];tally['rows']+=1;tally['empty_binary64_prefix']+=n==0
            tally['empty_resident64_prefix']+=resident is None
            neighbors=[]
            if n==0 and domain is not None:
                ks=[BINARY64.floor_coordinate(domain.low),BINARY64.ceiling_coordinate(domain.high)]
                for k in ks:
                    prefix=BINARY64.value(k);result=operation(prefix,constant,op,p['model'])
                    neighbors.append(dict(prefix=bits(float(prefix)),calculated=bits(float(result))))
                if not binary64(neighbors[0]['calculated'])<y<binary64(neighbors[1]['calculated']):
                    raise ValueError('binary64 prefix gap not certified by its actual neighboring results')
            constraints.append(dict(profile=index,required_prefix=domain.json() if domain is not None else None,
                binary64_prefix_count=n,resident64_prefix_count=0 if resident is None else resident[1]-resident[0]+1,
                adjacent_binary64_outputs=neighbors))
        records.append(dict(input=row['input'],density=row['density'],constraints=constraints))
    result=dict(format='dnacalc-gamma-unit-final-consumer-preimages/1',complete=True,profiles=profiles,counts=counts,records=records,
        receipts=doc['receipts'],source_bindings={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source)},
        limits='Any positive prefix, without a log, exponential, coefficient or normalizer producer. '
        'A binary64 gap excludes that one final stored-prefix consumer; its exact neighboring '
        'calculated outputs certify the gap. A resident preimage is necessary only, not an '
        'ABI or private algorithm identification. No prefix is fitted or installed and '
        'these feasibility counts are never product pass counts.')
    write(output,result);print(json.dumps(list(zip(profiles,counts))),flush=True);return result


def scaled_density_consumers(checkout,study,bank,output,cache,offline=False):
    """Hold shape and exact x/beta fixed while changing the final scale."""
    from geomean import verified_study_answers
    from exact import operation
    import gzip
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete power-density anchors')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    profiles=[dict(kind=kind,model=model,prefix_role='post_alpha_product' if kind=='divide' else 'unscaled_prefactor')
              for kind in ('divide','stored_ratio_multiply') for model in ('strict','stored_pc64')]+[
                  dict(kind='resident_ratio_multiply',model='stored_pc64',prefix_role='unscaled_prefactor')]
    rows=[];groups=[]
    for g,anchor in enumerate(doc['records']):
        a,y=binary64(anchor['shape']),binary64(anchor['input'])
        prefix=binary64(anchor['density']['bits'])*y/a
        if BINARY64.round(prefix)!=prefix:raise ValueError('anchor rescaling was not exact')
        groups.append(dict(shape=anchor['shape'],ratio=anchor['input'],density=anchor['density'],prefix_publication=bits(float(prefix))))
        for beta in (3.,5.,7.,10.):
            x=y*int(beta)
            if binary64(bits(float(x)))!=x:raise ValueError('scaled input rounded')
            args=[bits(float(x)),anchor['shape'],bits(beta),bits(0.)]
            rows.append(dict(group=g,args=args))
    original=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz';lookup={};wanted={tuple(r['args']) for r in rows}
    for line in gzip.open(original,'rt'):
        raw=json.loads(line);aa=raw['args']
        if len(aa)==1 and len(aa[0])==4 and all(a['kind']=='number' for a in aa[0]):
            key=tuple(a['bits_hex'].removeprefix('0x').lower() for a in aa[0])
            if key in wanted:lookup[key]=raw
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-scaled-consumer-questions/1',complete=True,
        rows=rows,groups=groups,profiles=profiles,original_counterparts=list(lookup.values()),
        source_bindings={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(original.resolve()):sha(original)},
        selection='Four fixed non-dyadic scales for each prior power-of-two input/shape. x/beta is exact and unchanged; original matching bank questions are reused, never re-asked.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    missing=[r for r in rows if tuple(r['args']) not in lookup]
    fresh=service.ask('GAMMA.DIST',[[client.number(x) for x in r['args']] for r in missing])['rows'] if missing else []
    answers={tuple(r['args']):a for r,a in zip(missing,fresh)}
    for key,raw in lookup.items():
        if raw['outcome']['kind']!='number':raise ValueError('nonnumeric historical scaled density')
        answers[key]=dict(type='number',bits=raw['outcome']['bits_hex'].removeprefix('0x'))
    accum=[]
    for g in groups:
        d,y=binary64(g['density']['bits']),binary64(g['ratio']);r=binary64(g['prefix_publication'])
        accum.append([inverse(Span(d,d),y,'div',True,p['model']) if p['kind']=='divide'
                      else BINARY64.pullback(Span(r,r)) for p in profiles])
    dead=[[v is None for v in group] for group in accum];records=[]
    for row in rows:
        answer=answers[tuple(row['args'])]
        if answer['type']!='number':raise ValueError('nonnumeric scaled-density answer')
        x,a=map(binary64,row['args'][:2]);y=binary64(answer['bits']);domains=[]
        for index,p in enumerate(profiles):
            if p['kind']=='divide':constant=x;op='div'
            else:
                constant=(EXTENDED64.round(a/x) if p['kind']=='resident_ratio_multiply'
                          else operation(a,x,'div',p['model']));op='mul'
            domain=inverse(Span(y,y),constant,op,True,p['model']);domains.append(domain.json() if domain is not None else None)
            gi=row['group']
            if not dead[gi][index]:
                accum[gi][index]=None if domain is None else accum[gi][index].meet(domain)
                if accum[gi][index] is None:dead[gi][index]=True
        records.append(dict(row,answer=answer,original_bank=lookup.get(tuple(row['args'])),required_prefixes=domains))
    counts=[Counter() for _ in profiles]
    for gi,g in enumerate(groups):
        g['prefix_constraints']=[]
        for i,p in enumerate(profiles):
            domain=accum[gi][i];b=BINARY64.coordinates(domain) if domain is not None else None;e=EXTENDED64.coordinates(domain) if domain is not None else None
            counts[i]['groups']+=1;counts[i]['binary64_nonempty']+=b is not None;counts[i]['resident64_nonempty']+=e is not None
            g['prefix_constraints'].append(dict(required_shared_prefix=domain.json() if domain is not None else None,
                binary64_count=0 if b is None else b[1]-b[0]+1,resident64_count=0 if e is None else e[1]-e[0]+1))
    result=dict(format='dnacalc-gamma-scaled-final-consumer-preimages/1',complete=True,profiles=profiles,counts=counts,groups=groups,records=records,
        receipts=doc['receipts']+service.receipts,source_bindings={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source),str(original.resolve()):sha(original)},
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='One arbitrary shared prefix per fixed shape/exact x/beta. For division it is '
        'the post-alpha product S, constrained directly by old D=S/y; for multiplication it '
        'is unscaled R, whose publication is fixed by the exact alpha/y power-of-two anchor. '
        'This allows every possible rounded alpha product even for fractional shapes. '
        'Both binary64 and resident64 possibilities '
        'are counted without choosing a value. Original bank counterparts retain their own '
        'environment. A contradiction refutes only that shared-prefix final-consumer model; '
        'a nonempty intersection identifies neither a prefix producer nor a complete gamma '
        'algorithm. No store, coefficient, threshold or per-row helper is installed.')
    write(output,result);print(json.dumps(dict(new_calls=len(fresh),original_bank_calls=len(lookup),profiles=list(zip(profiles,counts)))),flush=True);return result


def density_scale_holdout(checkout,bank,output,cache,offline=False):
    """Fresh fractional shapes with alpha/x an exact power of two at the anchor."""
    import gzip
    shapes=(.375,.625,.875,1.125,1.375,1.625,2.25,2.75,3.5,5.5,7.5,12.5,24.,48.)
    rows=[[bits(a*2.**e),bits(a),bits(1.),bits(0.)] for a in shapes for e in (-4,-1,0,1,2)]
    original=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz';wanted={tuple(r) for r in rows};known={}
    for line in gzip.open(original,'rt'):
        raw=json.loads(line);aa=raw['args']
        if len(aa)==1 and len(aa[0])==4 and all(a['kind']=='number' for a in aa[0]):
            key=tuple(a['bits_hex'].removeprefix('0x').lower() for a in aa[0])
            if key in wanted:known[key]=raw
    anchors=Path(str(output)+'.anchors.json')
    write(str(anchors)+'.questions.json',dict(format='dnacalc-gamma-scale-holdout-questions/1',complete=True,
        rows=rows,shapes=shapes,original_bank_counterparts=list(known.values()),source_bindings=STAGE_BINDINGS,
        selection='Fourteen new fixed fractional/integer shapes, five exact alpha/x powers of two. All historical matches are reused before asking new anchors; every anchor then gets all four predeclared scales.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    missing=[r for r in rows if tuple(r) not in known]
    fresh=service.ask('GAMMA.DIST',[[client.number(x) for x in r] for r in missing])['rows'] if missing else []
    answers={tuple(r):a for r,a in zip(missing,fresh)}
    for key,raw in known.items():
        if raw['outcome']['kind']!='number':raise ValueError('nonnumeric original anchor')
        answers[key]=dict(type='number',bits=raw['outcome']['bits_hex'].removeprefix('0x'))
    records=[]
    for r in rows:
        a=answers[tuple(r)]
        if a['type']!='number' or value(a['bits'])<=0:raise ValueError('positive normal holdout anchor required')
        records.append(dict(shape=r[1],input=r[0],density=a,original_bank=known.get(tuple(r))))
    write(anchors,dict(format='dnacalc-gamma-scale-holdout-anchors/1',complete=True,records=records,
        receipts=service.receipts,source_bindings={**STAGE_BINDINGS,str(original.resolve()):sha(original)},
        questions_sha256=sha(str(anchors)+'.questions.json'),
        limits='Original positive density anchors only; exact power-of-two scaling constrains '
        'the prefix publication without a logarithm, exponential or gamma approximation. '
        'No response selects a shape, scale, coefficient or store programme.'))
    return scaled_density_consumers(checkout,anchors,bank,output,cache,offline)


def density_exponents(checkout,bank,output,cache,offline=False):
    """Fixed normal-input, ratio-overflow and zero controls before a consumer change."""
    import gzip
    shapes=(.375,.5,1.,1.25,2.,4.,16.,64.)
    xs=[0.,2.**-1022,2.**-1000,2.**-500,2.**-100,2.**-10,1.,2.**100,2.**500,value('7fefffffffffffff')]
    rows=[[bits(x),bits(a),bits(b),bits(0.)] for a in shapes for x in xs for b in (2.**-1022,1.,value('7fefffffffffffff'))]
    for a in (4.,16.,64.):
        center=int(bits(a/value('7fefffffffffffff')),16)
        for k in range(-4,5):
            x=f'{center+k:016x}'
            if value(x)>=2.**-1022:
                rows.extend([[x,bits(a),bits(b),bits(0.)] for b in (2.**-1022,1.,value('7fefffffffffffff'))])
    rows=list(dict.fromkeys(map(tuple,rows)));wanted=set(rows);known={}
    original=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz'
    for line in gzip.open(original,'rt'):
        raw=json.loads(line);aa=raw['args']
        if len(aa)==1 and len(aa[0])==4 and all(a['kind']=='number' for a in aa[0]):
            key=tuple(a['bits_hex'].removeprefix('0x').lower() for a in aa[0])
            if key in wanted:known[key]=raw
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-density-exponent-questions/1',complete=True,
        rows=rows,original_bank_counterparts=list(known.values()),source_bindings=STAGE_BINDINGS,
        selection='Fixed positive shapes, zero and normal exponent edges, extreme scales and adjacent normal doubles around alpha/MAX. No subnormal input or output-dependent cutoff.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    missing=[r for r in rows if r not in known]
    fresh=service.ask('GAMMA.DIST',[[client.number(x) for x in r] for r in missing])['rows'] if missing else []
    records=[dict(args=list(r),answer=a,original_bank=None) for r,a in zip(missing,fresh)]
    for r,raw in known.items():
        o=raw['outcome'];a=(dict(type='number',bits=o['bits_hex'].removeprefix('0x')) if o['kind']=='number'
                          else dict(type='error',code=o['code']))
        records.append(dict(args=list(r),answer=a,original_bank=raw))
    result=dict(format='dnacalc-gamma-density-exponents/1',complete=True,records=records,receipts=service.receipts,
        source_bindings={**STAGE_BINDINGS,str(original.resolve()):sha(original)},questions_sha256=sha(str(output)+'.questions.json'),
        limits='Original numeric/error readings at predeclared admitted input extremes. No prediction '
        'is an Excel answer; exponent/ratio observations justify no unmeasured numerical cutoff.')
    write(output,result);print(json.dumps(dict(new_calls=len(fresh),original_calls=len(known),kinds=Counter(r['answer']['type'] for r in records))),flush=True);return result


def classify_consumer_losses(prefixes,output):
    """Classify every verified old match without inventing a private prefix."""
    from exact import operation
    source=Path(prefixes);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete registered prefix replay')
    records=[];counts=Counter()
    for row in doc['records']:
        if row['classification']!='prefix_replay_verified':raise ValueError('unresolved typed loss requires a separate route audit')
        old=row['original'];want=old['expected'];r=value(row['prefix']);factor=value(row['stored_factor'])
        if not want.startswith('0x') or not math.isfinite(r) or not math.isfinite(factor) or factor==0:
            classification='exception_or_scale_domain';constraints=[];alternatives=[]
        else:
            target=binary64(want);rr=binary64(row['prefix'])
            parse=lambda a:binary64(bits(float(a))) if isinstance(a,bool) else binary64(a)
            x,a=map(parse,old['args'][:2]);constraints=[];alternatives=[];any_binary64=False
            for factor_model in ('strict','stored_pc64'):
                f=operation(a,x,'div',factor_model)
                if factor_model=='stored_pc64' and f!=binary64(row['stored_factor']):raise ValueError('independent factor and registered replay differ')
                for multiply_model in ('strict','stored_pc64'):
                    predicted=operation(rr,f,'mul',multiply_model)
                    domain=inverse(Span(target,target),f,'mul',True,multiply_model)
                    coords=BINARY64.coordinates(domain) if domain is not None else None
                    any_binary64|=coords is not None
                    possible=domain is not None and domain.contains(rr)
                    if possible!=(predicted==target):raise ValueError('prefix inverse and actual arithmetic disagree')
                    if possible:alternatives.append(dict(factor_word=factor_model,multiply_word=multiply_model))
                    constraints.append(dict(factor_word=factor_model,multiply_word=multiply_model,
                        factor=bits(float(f)),prediction=bits(float(predicted)),
                        required_prefix=domain.json() if domain is not None else None,
                        binary64_prefix_count=0 if coords is None else coords[1]-coords[0]+1,
                        current_prefix_possible=possible))
            if any(p['factor_word']==p['multiply_word']=='stored_pc64' for p in alternatives):
                raise ValueError('listed loss still matches through the current consumer')
            classification=('arithmetic_words_remain_open' if alternatives else
                'current_prefix_outside_all_declared_word_preimages' if any_binary64 else
                'prefix_return_storage_remains_open')
        counts[classification]+=1;records.append(dict(row,loss_class=classification,word_alternatives=alternatives,prefix_constraints=constraints))
    result=dict(format='dnacalc-gamma-consumer-loss-classification/1',complete=True,counts=counts,records=records,
        source_bindings={**STAGE_BINDINGS,str(source.resolve()):sha(source)},
        limits='Each old matched return and new registered return is verified by the parent '
        'cleared-prefix replay. Necessary preimages allow both quotient/multiply words. '
        'Word alternatives are explicitly unresolved, never silently attributed to a prefix. '
        'Otherwise the retained prefix lies outside every declared preimage; its logarithm, '
        'exponential, normalizer or return-storage cause is not uniquely identified. Original '
        'answers stay unchanged and every loss stays a failure; no correction or selector is installed.')
    write(output,result);print(json.dumps(dict(rows=len(records),counts=counts)),flush=True);return result


def gammaln_near_stores(checkout,manifest,output,cache,library,offline=False,count=4096):
    """Resolve the existing fixed polynomial crossings with adjacent inputs."""
    import struct
    try:
        import mpmath as mp
    except ImportError:
        sys.path.insert(0,'/home/ubuntu/claude-work/gamma/pylib')
        import mpmath as mp
    mp.mp.prec=200
    source=Path(manifest);prior=json.loads(source.read_text())
    if prior.get('complete') is not True:raise ValueError('incomplete earlier edge manifest')
    frozen={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(Path(library).resolve()):sha(library),
            str(Path(mp.__file__).resolve()):sha(mp.__file__)}
    windows=[];questions=[]
    for edge in prior['edges']:
        center=int(edge['estimated_edge'],16);xs=[f'{center+i-count//2:016x}' for i in range(count)]
        windows.append(dict(polynomial=edge['polynomial'],level=edge['level'],center=edge['estimated_edge'],inputs=xs))
        questions.extend([[dict(type='number',bits=x)] for x in xs])
    write(str(output)+'.questions.json',dict(format='dnacalc-gammaln-near-store-questions/1',complete=True,
        windows=windows,source_bindings=frozen,selection='Four earlier polynomial crossings, each centered on4096 consecutive binary64 inputs. No new answer selects a center or a model.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('GAMMALN',questions)['rows'];native=Native(library);records=[]
    graphs=prior['graphs'];remainder=mp.mpf(0)
    for wi,window in enumerate(windows):
        center=mp.mpf(value(window['center']));f0=mp.loggamma(center);f1=mp.digamma(center);f2=mp.polygamma(1,center)
        samples=[]
        for index,(x,a) in enumerate(zip(window['inputs'],answers[wi*count:(wi+1)*count])):
            if a['type']!='number':raise ValueError('nonnumeric positive gamma observation')
            dx=mp.mpf(value(x))-center;truth=f0+f1*dx+f2*dx*dx/2
            unit=mp.mpf(math.ulp(value(a['bits'])));remainder=max(remainder,abs(dx)**3/(6*unit))
            predictions={name:f'{native.publish(native.evaluate(graph,[int(x,16)],0x133f,mask=1)[:10]):016x}' for name,graph in graphs.items()}
            ys=dict(Excel=a['bits'],**predictions)
            residuals={name:float((mp.mpf(value(raw))-truth)/unit) for name,raw in ys.items()}
            samples.append(dict(offset=index-count//2,input=x,excel=a,predictions=predictions,error_in_result_ulps=residuals))
        blocks=[]
        for start in range(0,count,128):
            part=samples[start:start+128]
            blocks.append(dict(first_offset=part[0]['offset'],last_offset=part[-1]['offset'],
                means={name:sum(p['error_in_result_ulps'][name] for p in part)/len(part) for name in ys}))
        records.append(dict(polynomial=window['polynomial'],level=window['level'],center=window['center'],samples=samples,blocks=blocks))
    result=dict(format='dnacalc-gammaln-near-store-windows/1',complete=True,records=records,
        source_bindings=frozen,receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        quadratic_remainder_bound_in_ulps=str(remainder),
        limits='Exact adjacent-input observations at previously chosen centers; block means are descriptive, '
        'not independent random samples or a fitted change point. The old coefficient table and '
        'all three predeclared store graphs remain frozen. On2<=x<=3, abs(psi2(x))<1 bounds '
        'the stated Taylor remainder. Neither a table nor a store change follows from a pattern alone.')
    write(output,result);print(json.dumps(dict(windows=len(records),calls=len(answers),remainder=str(remainder))),flush=True);return result


def scaled_quotient_words(checkout, output, cache, library, offline=False, count=64):
    """Independent public anchors at strict versus PC64 quotient ties.

    Select input divisions before any density answer. Shape is twice the
    PC64 quotient, making its unit-scale prefactor extraction exact. The
    competing strict anchor retains its whole exact consumer preimage.
    """
    import random
    from exact import operation
    native = Native(library)
    randomizer = random.Random(0x31a07)
    cases=[];attempts=0
    while len(cases)<count:
        attempts+=1
        if attempts>1000000:
            raise ValueError('bounded primitive tie mining exhausted')
        raw=[f'{(randomizer.randrange(1022,1026)<<52)|randomizer.getrandbits(52):016x}' for _ in range(2)]
        x,beta=raw
        strict=bits(value(x)/value(beta))
        extended=native.leaf('div',load64(int(x,16)),load64(int(beta,16)),0x133f)
        pc64=f'{native.publish(extended):016x}'
        if strict==pc64:
            continue
        xx,bb=binary64(x),binary64(beta)
        if bits(float(operation(xx,bb,'div','strict')))!=strict or bits(float(operation(xx,bb,'div','stored_pc64')))!=pc64:
            raise ValueError('native and rational quotient discriminants differ')
        shape=bits(2*value(pc64))
        cases.append(dict(x=x,beta=beta,shape=shape,strict_quotient=strict,pc64_quotient=pc64))
    frozen={**STAGE_BINDINGS,str(Path(library).resolve()):sha(library)}
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    anchors=[[client.number(y),client.number(c['shape']),client.number(1.),client.number(0.)]
             for c in cases for y in (c['pc64_quotient'],c['strict_quotient'])]
    questions=[[client.number(c['x']),client.number(c['shape']),client.number(c['beta']),client.number(0.)]for c in cases]
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-scaled-quotient-questions/1',complete=True,
        cases=cases,primitive_attempts=attempts,anchors=anchors,scaled_questions=questions,source_bindings=frozen,
        selection='Fixed RNG primitive division ties only, before any density answer. '
                  'Shape twice the PC64 quotient makes that anchor scale exact. '
                  'Every strict-anchor prefix preimage is retained; no producer is assumed.'))
    anchor_answers=service.ask('GAMMA.DIST',anchors)['rows']
    scaled_answers=service.ask('GAMMA.DIST',questions)['rows']
    records=[];counts=Counter()
    for index,c in enumerate(cases):
        pc64_anchor,strict_anchor=anchor_answers[2*index:2*index+2]
        target=scaled_answers[index]
        if any(a['type']!='number' for a in (pc64_anchor,strict_anchor,target)):
            raise ValueError('positive finite density controls must be numeric')
        a,x,y=map(binary64,(c['shape'],c['x'],c['strict_quotient']))
        prefix=binary64(pc64_anchor['bits'])/2
        if BINARY64.round(prefix)!=prefix:
            raise ValueError('anchor rescaling was not exact')
        scaled_factor=operation(a,x,'div','stored_pc64')
        predicted=bits(float(operation(prefix,scaled_factor,'mul','stored_pc64')))
        anchor_factor=operation(a,y,'div','stored_pc64')
        possible=inverse(Span(binary64(strict_anchor['bits']),binary64(strict_anchor['bits'])),anchor_factor,'mul',True,'stored_pc64')
        coordinates=BINARY64.coordinates(possible) if possible else None
        if coordinates and coordinates[1]-coordinates[0]>4:
            raise ValueError('unexpected unbounded strict prefix publication')
        strict_prefixes=[] if coordinates is None else [BINARY64.value(k)for k in range(coordinates[0],coordinates[1]+1)]
        strict_outputs=sorted({bits(float(operation(r,scaled_factor,'mul','stored_pc64')))for r in strict_prefixes})
        target_domain=inverse(Span(binary64(target['bits']),binary64(target['bits'])),scaled_factor,'mul',True,'stored_pc64')
        shared_domain=possible.meet(target_domain) if possible and target_domain else None
        resident_coordinates=EXTENDED64.coordinates(shared_domain) if shared_domain else None
        counts['rows']+=1;counts['pc64_anchor_prediction_matches']+=predicted==target['bits']
        counts['strict_anchor_has_binary64_prefix']+=coordinates is not None
        counts['strict_anchor_scaled_possible']+=target['bits']in strict_outputs
        counts['pc64_only']+=predicted==target['bits']and target['bits']not in strict_outputs
        counts['strict_only']+=predicted!=target['bits']and target['bits']in strict_outputs
        counts['neither']+=predicted!=target['bits']and target['bits']not in strict_outputs
        counts['strict_anchor_shared_resident_possible']+=resident_coordinates is not None
        counts['strict_anchor_shared_real_possible']+=shared_domain is not None
        counts['pc64_only_even_resident_prefix']+=predicted==target['bits']and resident_coordinates is None
        counts['pc64_only_even_real_prefix']+=predicted==target['bits']and shared_domain is None
        records.append(dict(c,pc64_anchor=pc64_anchor,strict_anchor=strict_anchor,scaled_answer=target,
            extracted_pc64_prefix=bits(float(prefix)),scaled_factor=bits(float(scaled_factor)),
            pc64_prediction=predicted,strict_prefix_domain=possible.json()if possible else None,
            strict_prefixes=[bits(float(r))for r in strict_prefixes],strict_predictions=strict_outputs,
            strict_shared_prefix_domain=shared_domain.json()if shared_domain else None,
            strict_shared_resident_count=0 if resident_coordinates is None else resident_coordinates[1]-resident_coordinates[0]+1))
    result=dict(format='dnacalc-gamma-scaled-quotient-words/1',complete=True,records=records,
        counts=counts,receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Conditional division-word discrimination through the independently measured stored alpha/x consumer. '
               'Each competing quotient uses its own original public density anchor, with every '
               'permitted binary64 prefix retained. Separate shared-prefix intersections also allow '
               'an arbitrary resident PC64 or real prefix, so discriminants conditional on prefix '
               'storage are distinguished from robust exclusions. A complete match identifies '
               'this sampled quotient publication, not a private prefix producer or a resident-quotient route. '
               'Primitive ties alone select inputs; no answer or match count selects a case.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def scaled_prefix_storage(checkout,study,prior_cache,output,cache,offline=False):
    """Use primitive multiplication ties to resolve one shared prefix's lattice."""
    import random
    from fractions import Fraction as Q
    from geomean import verified_study_answers
    from exact import operation
    source=Path(study);doc=json.loads(source.read_text())
    _,bindings,_=verified_study_answers(doc,Path(prior_cache))
    frozen={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source)}
    rng=random.Random(3111);groups=[];rows=[];excluded=[]
    for anchor in doc['records']:
        a,y=binary64(anchor['shape']),binary64(anchor['input'])
        if y not in (Q(1,4),Q(1),Q(2),Q(4)):continue
        r=binary64(anchor['density']['bits'])*y/a
        if BINARY64.round(r)!=r:raise ValueError('prefix publication extraction rounded')
        gi=len(groups);groups.append(dict(shape=anchor['shape'],ratio=anchor['input'],
            anchor=anchor['density'],prefix_publication=bits(float(r))))
        selected={-1:0,1:0}
        for trial in range(500000):
            xbits=((1023+rng.randint(-8,8))<<52)|rng.getrandbits(52)
            x=binary64(f'{xbits:016x}');beta=x/y
            if binary64(bits(float(beta)))!=beta:raise ValueError('holding scaled input rounded')
            factor=operation(a,x,'div','stored_pc64')
            strict=operation(r,factor,'mul','strict');stored=operation(r,factor,'mul','stored_pc64')
            if strict==stored:continue
            sign=1 if stored>strict else -1
            if selected[sign]==4:continue
            selected[sign]+=1
            rows.append(dict(group=gi,args=[f'{xbits:016x}',anchor['shape'],bits(float(beta)),bits(0.)],
                stored_ratio=bits(float(factor)),strict_stored_prefix=bits(float(strict)),
                pc64_stored_prefix=bits(float(stored)),direction=sign,selection_trial=trial))
            if all(n==4 for n in selected.values()):break
        if not all(n==4 for n in selected.values()):
            excluded.append(dict(group=gi,selected=selected,trials=500000,
                reason='Fixed bounded primitive search incomplete; no impossibility claim.'))
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-prefix-storage-questions/1',complete=True,
        groups=groups,rows=rows,excluded=excluded,source_bindings=frozen,
        selection='Thirteen inherited dyadic shape/input anchors. Seed3111; first four final-product strict/PC64 ties in each direction per anchor, before any new answer. x/beta stays exactly the original y.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('GAMMA.DIST',[[client.number(v) for v in row['args']] for row in rows])['rows']
    accum=[[BINARY64.pullback(Span(binary64(g['prefix_publication']),binary64(g['prefix_publication'])))
            for word in ('strict','stored_pc64')] for g in groups]
    dead=[[False,False] for g in groups];counts=[Counter(),Counter()]
    for row,answer in zip(rows,answers):
        row['answer']=answer;row['required_prefixes']=[]
        gi=row['group'];factor=binary64(row['stored_ratio'])
        for i,word in enumerate(('strict','stored_pc64')):
            domain=(inverse(Span(binary64(answer['bits']),binary64(answer['bits'])),factor,'mul',True,word)
                    if answer['type']=='number' else None)
            row['required_prefixes'].append(domain.json() if domain else None)
            counts[i]['calls']+=1
            prediction=row['strict_stored_prefix' if word=='strict' else 'pc64_stored_prefix']
            counts[i]['stored_anchor_matches']+=answer==dict(type='number',bits=prediction)
            if not dead[gi][i]:
                accum[gi][i]=accum[gi][i].meet(domain) if domain else None
                dead[gi][i]=accum[gi][i] is None
    for gi,g in enumerate(groups):
        g['constraints']=[]
        for i,word in enumerate(('strict','stored_pc64')):
            domain=None if dead[gi][i] else accum[gi][i]
            b=BINARY64.coordinates(domain) if domain else None
            e=EXTENDED64.coordinates(domain) if domain else None
            counts[i]['groups']+=1;counts[i]['binary64_possible']+=b is not None
            counts[i]['resident64_possible']+=e is not None
            counts[i]['resident64_unique']+=e is not None and e[0]==e[1]
            g['constraints'].append(dict(word=word,domain=domain.json() if domain else None,
                binary64_count=0 if b is None else b[1]-b[0]+1,resident64_count=0 if e is None else e[1]-e[0]+1,
                unique_resident_is_anchor=e is not None and e[0]==e[1] and EXTENDED64.value(e[0])==binary64(g['prefix_publication'])))
    result=dict(format='dnacalc-gamma-shared-prefix-storage/1',complete=True,
        groups=groups,records=rows,counts=counts,excluded=excluded,receipts=service.receipts,
        source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Conditional shared prefix for fixed shape and exactly unchanged x/beta; the '
               'final alpha/x is stored from PC64. The old dyadic anchor fixes its binary64 '
               'publication without selecting hidden low bits. A surviving lattice value is '
               'a necessary constraint, never a replacement prefix. A unique resident value '
               'at the anchor would identify the sampled value, not a literal machine store. '
               'Empty common cells refute this common-prefix/consumer composition. No '
               'coefficient or model is chosen by counts; no application code is changed.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def final_sum_quantization(study, edges, output):
    """Coefficient-residue consequences of the two numerator binade jumps.

    For a stored product on lattice D, adding one fixed constant and then
    storing on D/2 or D changes its rounding offset by zero or half D.
    This algebra does not use the unknown higher polynomial coefficients.
    The empirical jump remains a descriptive observation, not a fitted table.
    """
    from fractions import Fraction as Q
    from exact import operation
    source, edge_path = Path(study), Path(edges)
    doc, manifest = json.loads(source.read_text()),json.loads(edge_path.read_text())
    rows=[];allowed=set(range(4096));ambiguous=set()
    p0=manifest['tables']['P'][0]
    for window in doc['records']:
        if window['polynomial']!='P':
            continue
        samples=window['samples'];n=len(samples)//2
        jump=sum(r['error_in_result_ulps']['Excel'] for r in samples[n:])/n-sum(r['error_in_result_ulps']['Excel'] for r in samples[:n])/n
        x=binary64(window['center']);den=Q(0)
        for raw in reversed(manifest['tables']['Q']):den=den*x+binary64(raw)
        level=window['level'];coarse=Q(level,1<<52);fine=coarse/2
        unit=Q(math.ulp(value(samples[n]['excel']['bits'])))
        quantum=fine*(x-2)/den/unit
        candidates=set();ties=set();residue_shifts={}
        for residue in range(4096):
            constant=-2-Q(residue,1<<51)
            errors={False:set(),True:set()}
            # Both preceding-product parities on each side. The integer
            # part of c0 is a whole number of coarse lattice units, so it
            # only translates the arbitrary preceding product.
            for k in range(-3,6):
                product=level+2+k*coarse;exact_sum=product+constant
                stored=operation(product,constant,'add','stored_pc64')
                errors[exact_sum>=level].add(stored-exact_sum)
            shifts={(after-before)/fine for before in errors[False] for after in errors[True]}
            residue_shifts[residue]=sorted(str(q)for q in shifts)
            if len(shifts)!=1:
                ties.add(residue);continue
            shift=next(iter(shifts))
            if shift and (shift>0)==(jump>0):candidates.add(residue)
        allowed &= candidates;ambiguous |= ties
        rows.append(dict(level=level,center=window['center'],observed_half_mean_jump_ulps=jump,
            one_quantum_in_result_ulps=float(quantum),observed_quanta=jump/float(quantum),
            current_table_constant_shift_quanta=residue_shifts[int(p0,16)&4095],
            permitted_non_tie_residues=sorted(candidates),ambiguous_residues=sorted(ties),
            ambiguous_shift_quanta={str(r):residue_shifts[r]for r in sorted(ties)}))
    result=dict(format='dnacalc-gammaln-final-sum-quantization/1',complete=True,records=rows,
        joint_non_tie_residues=sorted(allowed),current_constant=p0,current_low12=int(p0,16)&4095,
        source_bindings={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(edge_path.resolve()):sha(edge_path)},
        limits='Conditional on a stored numerator last sum, an already stored preceding product '
               'and the measured negative/positive jump signs. Only c0 modulo4096 of its binary64 '
               'units is constrained; no coefficient is selected. Both product parities are '
               'evaluated through exact PC64 addition and binary64 publication. Exact ties and '
               'PC64 double-rounding midpoint neighbors remain ambiguous. Quantum sizes use the frozen '
               'approximate Q only as a scale diagnostic; no Q coefficient is inferred. '
               'Descriptive residual means establish neither independent random samples nor a '
               'unique private implementation. The current table cannot exhibit either observed '
               'nonzero jump merely by adding the final store; a store-only numerical change '
               'still needs an independent classification of its losses.')
    write(output,result)
    print(json.dumps(dict(joint_residue_count=len(allowed),current_low12=result['current_low12'],
        observed_quanta=[r['observed_quanta'] for r in rows])),flush=True)
    return result


def base2_prefactor_constraints(study, output, cache, library, signed=False):
    """Shared-normalizer inverses for literal base-two gamma prefixes.

    Exact dyadic y and power-of-two shapes remove the power logarithm. The
    remaining literal programs use the x87 log2(e) ROM constant and one
    arbitrary normalizer per shape. This tests a different arithmetic family
    from a stored natural logarithm; no normalizer is fitted or installed.
    """
    from geomean import verified_study_answers
    from exact import operation
    source = Path(study)
    doc = json.loads(source.read_text())
    _, bindings, _ = verified_study_answers(doc,Path(cache))
    native = Native(library)
    frozen = {**STAGE_BINDINGS, **bindings, str(source.resolve()):sha(source),
              str(Path(library).resolve()):sha(library)}
    g = Builder(1)
    exponent = g.op('neg',0,store=False) if signed else 0
    whole = g.op('frndint',exponent,store=False)
    fraction = g.op('sub',exponent,whole,store=False)
    mantissa = g.op('add',g.one,g.op('f2xm1',fraction,store=False),store=False)
    positive = g.op('fscale',mantissa,whole,store=False)
    result = positive if signed else g.op('div',g.one,positive,store=False)
    graph = g.graph(result,result)
    l2e = decode(native.leaf('fldl2e',bytes(10),bytes(10),0x133f))
    profiles = [dict(form='sum_then_scale',sum_stored=s) for s in (False,True)]
    profiles += [dict(form='separate_scales',x_stored=x,c_stored=c) for x,c in
                 ((False,False),(True,False),(False,True),(True,True))]
    profiles += [dict(form='subtract_power_then_normalizer',difference_stored=s) for s in (False,True)]
    groups = {}
    records = []
    for row in doc['records']:
        a,x = binary64(row['shape']),binary64(row['input'])
        if x.numerator & (x.numerator-1) or x.denominator & (x.denominator-1):
            raise ValueError('exact dyadic powers required')
        e = x.numerator.bit_length()-x.denominator.bit_length()
        ae = a*e
        if BINARY64.round(ae) != ae:
            raise ValueError('power exponent rounded')
        cell = positive_extended_cell(native,graph,row['extracted_prefactor'],
            binary64(bits(-math.log2(value(row['extracted_prefactor'])))),decreasing=True)
        if cell is None:
            raise ValueError('empty base-two primitive preimage')
        interval = Span(*(decode(bytes.fromhex(cell[k])) for k in ('low','high')))
        group = groups.setdefault(row['shape'],[dict(domain=None,dead=False,rows=0) for _ in profiles])
        domains = []
        for index,p in enumerate(profiles):
            if p['form'] == 'sum_then_scale':
                # u = (x+c)*log2(e) - a*log2(x).
                domain = inverse(interval,-ae,'add',True,'resident_pc64')
                domain = inverse(domain,l2e,'mul',True,'resident_pc64') if domain else None
                domain = inverse(domain,x,'add',True,'stored_pc64' if p['sum_stored'] else 'resident_pc64') if domain else None
            elif p['form'] == 'separate_scales':
                # u = (x*log2(e) + c*log2(e)) - a*log2(x).
                yy = operation(x,l2e,'mul','stored_pc64' if p['x_stored'] else 'resident_pc64')
                domain = inverse(interval,-ae,'add',True,'resident_pc64')
                domain = inverse(domain,yy,'add',True,'resident_pc64') if domain else None
                domain = inverse(domain,l2e,'mul',True,'stored_pc64' if p['c_stored'] else 'resident_pc64') if domain else None
            else:
                # u = (x*log2(e)-a*log2(x)) + c*log2(e).
                yy = operation(x,l2e,'mul','resident_pc64')
                base = operation(yy,-ae,'add','stored_pc64' if p['difference_stored'] else 'resident_pc64')
                domain = inverse(interval,base,'add',True,'resident_pc64')
                domain = inverse(domain,l2e,'mul',True,'resident_pc64') if domain else None
            domains.append(domain.json() if domain else None)
            item = group[index]
            item['rows'] += 1
            if not item['dead']:
                item['domain'] = domain if item['domain'] is None else item['domain'].meet(domain) if domain else None
                if item['domain'] is None:
                    item['dead'] = True
                    item['first_contradiction'] = len(records)
        records.append(dict(shape=row['shape'],input=row['input'],density=row['density'],
            extracted_prefactor=row['extracted_prefactor'],base2_magnitude_cell=cell,
            normalizer_domains=domains))
    for shape, group in groups.items():
        for item in group:
            domain = item.pop('domain')
            coords = BINARY64.coordinates(domain) if domain else None
            item.update(domain=domain.json() if domain else None,
                binary64_count=0 if coords is None else coords[1]-coords[0]+1)
    result = dict(format='dnacalc-gamma-base2-prefix-constraints/1',complete=True,
        profiles=profiles,groups=groups,records=records,graph=graph,receipts=doc['receipts'],signed=signed,
        source_bindings=frozen,
        limits='Necessary shared-normalizer intervals only, from original dyadic density anchors. '
               'The declared signed or reciprocal base-two exponential and fixed PC64 arithmetic are assumptions. '
               'An empty intersection rejects only its declared program. '
               'A feasible interval identifies no coefficient, private routine or implementation.')
    write(output,result)
    print(json.dumps({shape:[item['binary64_count'] for item in group] for shape,group in groups.items()}),flush=True)
    return result


def stirling_stage_constraints(checkout, output, cache, library, offline=False, count=64, consumers=False):
    """Invert the last GAMMA product without choosing Stirling coefficients.

    Independent public LN and EXP leaves remove the corresponding producer
    assumptions. Each declared stage graph admits an arbitrary binary64 w(x).
    Common-window and monotone constraints are additional, explicitly
    conditional tests, not fitted values or a replacement for the application.
    """
    if not 2 <= count <= 256:
        raise ValueError('small adjacent windows required')
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    centers = (10., 12., 16., 20., 24., 30., 32., 36., 48., 64., 80., 100., 128., 160.)
    prefix_profiles = ['resident_product_resident_sqrt', 'stored_product_resident_sqrt',
                       'resident_product_stored_sqrt', 'stored_product_stored_sqrt']
    forms = ('constant_then_correction', 'correction_then_constant', 'stored_constant_correction',
             'resident_constant_correction', 'resident_constant_body', 'resident_square_body') if consumers else ('constant_then_correction',)
    profiles = [dict(prefix=p,finish=f) for p in prefix_profiles for f in forms]
    groups = [dict(center=bits(c), inputs=[f'{int(bits(c),16)+i:016x}' for i in range(-count//2,count//2+1) if i])
              for c in centers]
    # Stay in the high band even at its first exact boundary.
    groups[0]['inputs'] = [f'{int(bits(10.),16)+i:016x}' for i in range(1,count+1)]
    xs = [x for group in groups for x in group['inputs']]
    frozen = {**STAGE_BINDINGS, str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json', dict(format='dnacalc-gamma-stirling-stage-questions/1',
        complete=True, groups=groups, profiles=profiles, source_bindings=frozen,
        selection='Fourteen fixed integer and exponent boundaries, consecutive inputs in increasing order. '
                  'Four literal product/square-root store graphs precede every new answer. '
                  'All logarithms and exponentials are independently measured public calls.'))
    logs = dict(zip(xs,service.ask('LN',[[client.number(x)] for x in xs])['rows']))
    ems = dict(zip(xs,service.ask('EXP',[[client.number(-value(x))] for x in xs])['rows']))
    answers = dict(zip(xs,service.ask('GAMMA',[[client.number(x)] for x in xs])['rows']))
    def op(kind, a, b):
        return native.leaf(kind, a, b, 0x133f)
    def scalar(x):
        return load64(int(bits(x),16))
    def stored(a):
        return load64(native.publish(a))
    root_arguments = {}
    for x in xs:
        half = .5*value(x)
        root_arguments[x] = f'{native.publish(op("mul",scalar(half-math.floor(half)),load64(int(logs[x]["bits"],16)))):016x}'
    roots_unique = sorted(set(root_arguments.values()))
    roots = dict(zip(roots_unique,service.ask('EXP',[[client.number(x)] for x in roots_unique])['rows']))
    def power(x, n):
        product, square = None, load64(int(x,16))
        while n:
            if n & 1:
                product = square if product is None else stored(op('mul',product,square))
            n >>= 1
            if n:
                square = stored(op('mul',square,square))
        return product if product is not None else scalar(1.)
    records = []
    counts = [Counter() for _ in profiles]
    for group in groups:
        common = [None for _ in profiles]
        common_dead = [False for _ in profiles]
        monotone_upper = [None for _ in profiles]
        monotone_dead = [False for _ in profiles]
        group_records = []
        for x in group['inputs']:
            target = answers[x]
            if any(a['type'] != 'number' for a in (target,logs[x],ems[x],roots[root_arguments[x]])):
                raise ValueError('high-band leaves must be finite numbers')
            v = stored(op('mul',power(x,math.floor(.5*value(x))),load64(int(roots[root_arguments[x]]['bits'],16))))
            square_root = op('sqrt',load64(int(x,16)),scalar(0.))
            product = op('mul',v,load64(int(ems[x]['bits'],16)))
            domains = []
            for index, profile in enumerate(profiles):
                label, finish = profile['prefix'], profile['finish']
                pp = stored(product) if label.startswith('stored_product') else product
                qq = stored(square_root) if label.endswith('stored_sqrt') else square_root
                u = stored(op('div',pp,qq))
                square = stored(op('mul',u,v))
                constant = load64(0x40040d931ff62706)
                body = stored(op('mul',square,constant))
                interval = Span(binary64(target['bits']),binary64(target['bits']))
                if finish in ('constant_then_correction','resident_constant_body','resident_square_body'):
                    if finish == 'resident_constant_body':
                        body = op('mul',square,constant)
                    elif finish == 'resident_square_body':
                        body = stored(op('mul',op('mul',u,v),constant))
                    domain = inverse(interval,decode(body),'mul',True,'stored_pc64')
                elif finish == 'correction_then_constant':
                    interval = inverse(interval,decode(constant),'mul',True,'stored_pc64')
                    domain = inverse(interval,decode(square),'mul',True,'stored_pc64') if interval else None
                elif finish in ('stored_constant_correction','resident_constant_correction'):
                    interval = inverse(interval,decode(square),'mul',True,'stored_pc64')
                    model = 'stored_pc64' if finish.startswith('stored_') else 'resident_pc64'
                    domain = inverse(interval,decode(constant),'mul',True,model) if interval else None
                else:
                    raise ValueError('unknown fixed final association')
                coordinates = BINARY64.coordinates(domain) if domain else None
                counts[index]['rows'] += 1
                counts[index]['empty_binary64_w'] += coordinates is None
                if not common_dead[index]:
                    common[index] = domain if common[index] is None else common[index].meet(domain) if domain else None
                    common_dead[index] = common[index] is None
                if not monotone_dead[index]:
                    if coordinates is None:
                        monotone_dead[index] = True
                    else:
                        low, high = coordinates
                        upper = high if monotone_upper[index] is None else min(monotone_upper[index],high)
                        monotone_dead[index] = upper < low
                        monotone_upper[index] = upper
                domains.append(dict(body=f'{native.publish(body):016x}',resident_body=body.hex(),square=f'{native.publish(square):016x}',
                    required_w=domain.json() if domain else None,
                    binary64_count=0 if coordinates is None else coordinates[1]-coordinates[0]+1,
                    binary64_values=[] if coordinates is None else [bits(float(BINARY64.value(k)))
                        for k in range(coordinates[0],min(coordinates[1]+1,coordinates[0]+4))]))
            group_records.append(dict(input=x,answer=target,public_log=logs[x],public_negative_exp=ems[x],
                root_argument=root_arguments[x],public_root=roots[root_arguments[x]],
                powered=f'{native.publish(v):016x}',resident_sqrt=square_root.hex(),profiles=domains))
        results = []
        for index in range(len(profiles)):
            coordinates = BINARY64.coordinates(common[index]) if common[index] else None
            counts[index]['groups'] += 1
            counts[index]['common_w_possible'] += coordinates is not None
            counts[index]['monotone_w_possible'] += not monotone_dead[index]
            results.append(dict(common_w=common[index].json() if common[index] else None,
                common_binary64_count=0 if coordinates is None else coordinates[1]-coordinates[0]+1,
                monotone_w_possible=not monotone_dead[index]))
        records.append(dict(center=group['center'],samples=group_records,profiles=results))
    result = dict(format='dnacalc-gamma-stirling-stage-constraints/1',complete=True,
        profiles=profiles,counts=counts,records=records,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Exact final-product preimages with an arbitrary binary64 correction per input. '
               'Empty cells refute the declared preceding body and final product. '
               'Common or nonincreasing corrections across a window are stronger conditional hypotheses; '
               'a failed one does not refute an arbitrary nonmonotonic computed correction. '
               'The source integer-power route and stored public LN/EXP leaves stay fixed. '
               'No correction, coefficient or stage is selected by counts or installed.')
    write(output,result)
    print(json.dumps(dict(profiles=profiles,counts=counts)),flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout')
    p.add_argument('output')
    p.add_argument('--cache', required=True)
    p.add_argument('--library')
    p.add_argument('--count', type=int, default=32)
    p.add_argument('--offline', action='store_true')
    p.add_argument('--power-prefactor', action='store_true')
    p.add_argument('--power-primitives')
    p.add_argument('--density-primitives')
    p.add_argument('--shared-log-constraints')
    p.add_argument('--scaled-log-constraints')
    p.add_argument('--split-prefactor')
    p.add_argument('--poisson-pairs')
    p.add_argument('--gammaln-store-edges', action='store_true')
    p.add_argument('--gammaln-near-stores',metavar='MANIFEST')
    p.add_argument('--unit-shape-exp',nargs=2,metavar=('BANK','GRAPHS'))
    p.add_argument('--unit-prefactor',metavar='STUDY')
    p.add_argument('--unit-constants',metavar='STUDY')
    p.add_argument('--unit-consumers',metavar='STUDY')
    p.add_argument('--scaled-consumers',nargs=2,metavar=('STUDY','BANK'))
    p.add_argument('--scale-holdout',metavar='BANK')
    p.add_argument('--density-exponents',metavar='BANK')
    p.add_argument('--classify-losses',metavar='PREFIX_REPLAY')
    p.add_argument('--stirling-stages',action='store_true')
    p.add_argument('--stirling-consumers',action='store_true')
    p.add_argument('--base2-prefactor',metavar='STUDY')
    p.add_argument('--base2-signed',action='store_true')
    p.add_argument('--final-sum-quantization',nargs=2,metavar=('NEAR_STUDY','EDGE_MANIFEST'))
    p.add_argument('--prefix-storage',nargs=2,metavar=('STUDY','PRIOR_CACHE'))
    p.add_argument('--scaled-quotient-words',action='store_true')
    p.add_argument('--seed-phase-tracks',metavar='RETAINED_JSONL')
    p.add_argument('--radius',type=float,default=.001)
    p.add_argument('--seed-phase-validation',metavar='STUDY')
    p.add_argument('--retained-phase-dir')
    p.add_argument('--seed-phase-absolute',metavar='VALIDATED_STUDY')
    p.add_argument('--low-recurrence',action='store_true')
    p.add_argument('--low-shared-seed',metavar='STUDY')
    p.add_argument('--gammaln-low-composition',action='store_true')
    p.add_argument('--gammaln-low-tables',action='store_true')
    p.add_argument('--beta-shape-two',action='store_true')
    p.add_argument('--gammaln-log-controls',metavar='REPLAY')
    p.add_argument('--beta-log-replay',nargs=2,metavar=('STUDY','REPLAY'))
    p.add_argument('--beta-integer-unit-shape',action='store_true')
    p.add_argument('--beta-unit-prior',metavar='STUDY')
    p.add_argument('--beta-series-normalizer',action='store_true')
    p.add_argument('--beta-public-f-tails',metavar='BANK_TRACE')
    p.add_argument('--beta-f-unit-prefix',metavar='PRIOR_STUDY')
    p.add_argument('--beta-f-power-pairs',action='store_true')
    p.add_argument('--beta-resident-log',action='store_true')
    p.add_argument('--beta-strict-final',action='store_true')
    p.add_argument('--gammaln-table-frequencies',nargs=2,metavar=('WINDOWS','TABLES'))
    p.add_argument('--gammaln-stirling-shape',action='store_true')
    p.add_argument('--gammaln-stirling-proof',action='store_true')
    p.add_argument('--stirling-extra-inputs')
    p.add_argument('--seed',type=int,default=44005)
    p.add_argument('--audit-gammaln-log-change',nargs=2,metavar=('CONTROLS','REGISTERED'))
    p.add_argument('--density-normalizer-questions',metavar='BANK_TRACE')
    args = p.parse_args()
    # Keep the exact reader epoch before a later study extends this tool.
    snapshots=Path(args.output).resolve().parent/'producers'
    snapshots.mkdir(parents=True,exist_ok=True)
    for filename,digest in STAGE_BINDINGS.items():
        data=Path(filename).read_bytes()
        if hashlib.sha256(data).hexdigest()!=digest:raise ValueError('study source changed')
        destination=snapshots/(digest+'-'+Path(filename).name)
        if not destination.exists():
            with destination.open('xb')as stream:stream.write(data)
    if not 2 <= args.count <= (8192 if args.gammaln_stirling_proof else 4096 if args.gammaln_store_edges or args.gammaln_near_stores or args.seed_phase_tracks else 256):
        raise ValueError('small adjacent windows only')
    if not (args.gammaln_stirling_shape or args.seed_phase_tracks or args.seed_phase_validation or args.seed_phase_absolute or args.low_shared_seed or args.gammaln_log_controls or args.gammaln_table_frequencies or args.audit_gammaln_log_change or args.beta_public_f_tails) and not args.library:
        p.error('--library is required for instruction studies')
    if args.gammaln_stirling_shape:
        write(args.output,gammaln_stirling_shape(args.checkout))
    elif args.gammaln_stirling_proof:
        gammaln_stirling_proof(args.checkout,args.output,args.cache,args.library,args.count,args.seed,args.stirling_extra_inputs,args.offline)
    elif args.beta_f_unit_prefix:
        beta_f_unit_prefix(args.checkout,args.beta_f_unit_prefix,args.output,args.cache,args.library,args.offline,args.beta_f_power_pairs)
    elif args.beta_public_f_tails:
        beta_public_f_tails(args.checkout,args.beta_public_f_tails,args.output,args.cache,args.offline)
    elif args.density_normalizer_questions:
        density_normalizer_questions(args.checkout,args.density_normalizer_questions,args.output,args.cache,args.library,args.offline)
    elif args.audit_gammaln_log_change:
        audit_gammaln_log_change(*args.audit_gammaln_log_change,args.output)
    elif args.gammaln_table_frequencies:
        gammaln_table_frequencies(*args.gammaln_table_frequencies,args.output)
    elif args.beta_series_normalizer:
        beta_series_normalizer(args.checkout,args.output,args.cache,args.library,args.offline,args.count,args.beta_resident_log,args.beta_strict_final)
    elif args.beta_integer_unit_shape:
        beta_integer_unit_shape(args.checkout,args.output,args.cache,args.library,args.offline,args.count,args.beta_unit_prior)
    elif args.beta_log_replay:
        beta_log_replay(*args.beta_log_replay,args.output,args.library)
    elif args.gammaln_log_controls:
        gammaln_log_controls(args.checkout,args.gammaln_log_controls,args.output,args.cache,args.offline)
    elif args.beta_shape_two:
        beta_shape_two(args.checkout,args.output,args.cache,args.library,args.offline,args.count)
    elif args.gammaln_low_tables:
        gammaln_low_tables(args.checkout,args.output,args.cache,args.library,args.offline,args.count)
    elif args.gammaln_low_composition:
        gammaln_low_composition(args.checkout,args.output,args.cache,args.library,args.offline,args.count)
    elif args.low_shared_seed:
        low_gamma_shared_seed(args.low_shared_seed,args.output)
    elif args.seed_phase_absolute:
        seed_phase_absolute(args.seed_phase_absolute,args.output)
    elif args.low_recurrence:
        low_gamma_recurrence(args.checkout,args.output,args.cache,args.library,args.offline,args.count)
    elif args.seed_phase_validation:
        seed_phase_validation(args.seed_phase_validation,args.output,args.cache,args.retained_phase_dir)
    elif args.seed_phase_tracks:
        seed_phase_tracks(args.checkout,args.output,args.cache,args.seed_phase_tracks,args.offline,args.radius,args.count)
    elif args.scaled_quotient_words:
        scaled_quotient_words(args.checkout,args.output,args.cache,args.library,args.offline,args.count)
    elif args.prefix_storage:
        scaled_prefix_storage(args.checkout,args.prefix_storage[0],args.prefix_storage[1],args.output,args.cache,args.offline)
    elif args.final_sum_quantization:
        final_sum_quantization(args.final_sum_quantization[0],args.final_sum_quantization[1],args.output)
    elif args.base2_prefactor:
        base2_prefactor_constraints(args.base2_prefactor,args.output,args.cache,args.library,args.base2_signed)
    elif args.stirling_stages or args.stirling_consumers:
        stirling_stage_constraints(args.checkout,args.output,args.cache,args.library,args.offline,args.count,args.stirling_consumers)
    elif args.classify_losses:
        classify_consumer_losses(args.classify_losses,args.output)
    elif args.density_exponents:
        density_exponents(args.checkout,args.density_exponents,args.output,args.cache,args.offline)
    elif args.scale_holdout:
        density_scale_holdout(args.checkout,args.scale_holdout,args.output,args.cache,args.offline)
    elif args.scaled_consumers:
        scaled_density_consumers(args.checkout,args.scaled_consumers[0],args.scaled_consumers[1],args.output,args.cache,args.offline)
    elif args.unit_consumers:
        unit_shape_consumers(args.unit_consumers,args.output,args.cache)
    elif args.unit_constants:
        unit_shape_constants(args.unit_constants,args.output,args.cache,args.library)
    elif args.unit_prefactor:
        unit_shape_prefactor(args.checkout,args.unit_prefactor,args.output,args.cache,args.library,args.offline)
    elif args.unit_shape_exp:
        unit_shape_exponential(args.checkout,args.unit_shape_exp[0],args.unit_shape_exp[1],args.output,args.cache,args.library,args.offline)
    elif args.gammaln_near_stores:
        gammaln_near_stores(args.checkout,args.gammaln_near_stores,args.output,args.cache,args.library,args.offline,args.count)
    elif args.gammaln_store_edges:
        gammaln_store_edges(args.checkout, args.output, args.cache, args.library, args.offline, args.count)
    elif args.poisson_pairs:
        poisson_pairs(args.checkout, args.poisson_pairs, args.output, args.cache, args.library, args.offline)
    elif args.split_prefactor:
        split_prefactor(args.checkout, args.split_prefactor, args.output, args.cache, args.library, args.offline)
    elif args.scaled_log_constraints:
        shared_log_constraints(args.scaled_log_constraints,args.output,args.cache,args.library,scaled_log=True)
    elif args.shared_log_constraints:
        shared_log_constraints(args.shared_log_constraints, args.output, args.cache, args.library)
    elif args.density_primitives:
        power_primitives(args.density_primitives,args.output,args.cache,args.library,direct_density=True)
    elif args.power_primitives:
        power_primitives(args.power_primitives, args.output, args.cache, args.library)
    elif args.power_prefactor:
        power_prefactor(args.checkout, args.output, args.cache, args.library, args.offline)
    else:
        paired(args.checkout, args.output, args.cache, args.library, args.count, args.offline)
