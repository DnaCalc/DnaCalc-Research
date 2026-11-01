"""Exact unknown-normalizer constraints. Never select a denominator value."""
from pathlib import Path
import collections
import json
import sys
from fractions import Fraction as Q
from intervals import B64,E64,E53,Interval,inverse_op,inverse_add,number
from certify import exact_root


def trace_records(path):
    with Path(path).open() as stream:
        for line in stream:
            yield json.loads(line)


def normal_cdf_linear(trace, output):
    """Shared arbitrary coefficient feasibility, never choose or fit its bits."""
    import hashlib,itertools,time
    start=time.monotonic();output=Path(output);output.mkdir(exist_ok=False)
    groups=collections.defaultdict(list);counts=collections.Counter()
    for index,row in enumerate(trace_records(trace)):
        args=row['args']
        flag=args[1] if isinstance(args[1],bool) else number(args[1])
        if not flag:counts['density_control']+=1;continue
        x=Q(number(args[0]))
        if x==0 or abs(x)>Q(1,1<<27):counts['outside_fixed_tiny_domain']+=1;continue
        counts['scoped']+=1
        key=(json.dumps(row['reading'],sort_keys=True),row['source'],json.dumps(args[1]))
        groups[key].append((index,row,x))
    profiles=[]
    for form,sites in [(0,2),(1,3),(2,3)]:
        for models in itertools.product(range(5),repeat=sites):
            profiles.append(dict(id=len(profiles),form=form,models=list(models)))
    def span_for(row,x,p):
        y=Q(number(row['expected']));span=B64.preimage(Interval(y,y));ms=p['models']
        if p['form']==0:
            span=inverse_add(span,ms[1],Q(1,2))
        elif p['form']==1:
            span=inverse_op(span,ms[2],'mul',Q(1,2),True)
            span=inverse_add(span,ms[1],Q(1)) if span else None
        else:
            span=inverse_add(span,ms[2],Q(1,2))
            span=inverse_op(span,ms[1],'mul',Q(1,2),True) if span else None
        return inverse_op(span,ms[0],'mul',x,True) if span else None
    alive={p['id']:p for p in profiles};certificates=[]
    for key,rows in groups.items():
        if not alive:break
        # Independent input order: largest |z| first; ties retain source order.
        rows=sorted(rows,key=lambda r:(-abs(r[2]),r[0]))
        for cid,p in list(alive.items()):
            current=None;first=True;bounds=[]
            for index,row,x in rows:
                required=span_for(row,x,p)
                new=required if first else current.intersect(required) if current else None
                first=False
                if new is None or E64.preimage(new) is None:
                    ids=sorted(set(bounds+[index]))
                    cert=dict(profile=p,group=dict(reading=row['reading'],source=row['source'],flag=row['args'][1]),
                        witnesses=[dict(trace_index=i,witness=r)for i,r,_ in rows if i in ids],
                        reason='empty_real_intersection'if new is None else'empty_ext80_lattice',
                        remaining_interval=new.json()if new else None)
                    certificates.append(cert);del alive[cid];break
                if current is None:bounds=[index,index]
                else:
                    if new.lo!=current.lo or new.lc!=current.lc:bounds[0]=index
                    if new.hi!=current.hi or new.hc!=current.hc:bounds[1]=index
                current=new
    (output/'profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
    (output/'certificates.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in certificates))
    (output/'survivors.json').write_text(json.dumps(list(alive.values()),indent=2)+'\n')
    result=dict(complete=True,profiles=len(profiles),excluded=len(certificates),surviving_constraint_profiles=len(alive),
        counts=counts,exact_environment_source_flag_groups=len(groups),seconds=time.monotonic()-start,
        input_sha256=hashlib.sha256(Path(trace).read_bytes()).hexdigest(),generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        forms=['RNpub(RNa(.5+RNm(c*z)))','RNpub(RNs(.5*RNa(1+RNm(c*z))))','RNpub(RNa(.5+RNs(.5*RNm(c*z))))'],
        models=['strict','stored_pc64','resident_pc64','stored_pc53','resident_pc53'],
        scope='Fixed analytical tiny region 0<|z|<=2^-27; coefficient may be any common finite Ext80 value per exact reading/source/flag group. No coefficient value is selected. Not an assertion that Excel uses a linear branch.',
        claim='Exact necessary shared-coefficient intersections; empty real sets exclude even arbitrary real coefficients under the specified operation/store assumptions. Nonempty sets are unresolved constraints, never accepted programs.',
        constant_fitted=False,new_excel_answers=0)
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def normal_cdf_symmetry(trace, output):
    """Bounded sign/reconstruction relations with arbitrary shared latent values."""
    import hashlib,itertools,time
    start=time.monotonic();output=Path(output);output.mkdir(exist_ok=False)
    grouped=collections.defaultdict(dict);counts=collections.Counter()
    for index,row in enumerate(trace_records(trace)):
        a=row['args'];flag=a[1] if isinstance(a[1],bool)else number(a[1])
        if not flag:counts['density_control']+=1;continue
        bits=int(a[0],16);magnitude=bits&((1<<63)-1)
        if magnitude==0 or not row['expected'].startswith('0x'):counts['unsupported_control']+=1;continue
        key=(json.dumps(row['reading'],sort_keys=True),row['source'],str(a[1]),magnitude)
        grouped[key][bits>>63]=dict(trace_index=index,witness=row)
    pairs=[]
    for group in grouped.values():
        if len(group)!=2:counts['unpaired_inputs']+=len(group);continue
        positive,negative=group[0],group[1]
        yp,yn=Q(number(positive['witness']['expected'])),Q(number(negative['witness']['expected']))
        if not (Q(0)<yn<Q(1,2)<yp<Q(1)):
            counts['endpoint_or_nonordinary_pair_controls']+=1;continue
        pairs.append(dict(positive=positive,negative=negative))
    pairs.sort(key=lambda p:p['positive']['trace_index'])
    profiles=[]
    for form,n in [(0,1),(1,1),(2,2),(3,4)]:
        for models in itertools.product(range(5),repeat=n):profiles.append(dict(id=len(profiles),form=form,models=list(models)))
    assert len(profiles)==660
    def allowed(entry,p,negative):
        y=Q(number(entry['witness']['expected']));span=B64.preimage(Interval(y,y));ms=p['models']
        if p['form']==0:
            if not negative:
                span=inverse_add(span,ms[0],Q(1))
                if span:span=span.affine(Q(-1))
        elif p['form']==1:
            if negative:
                span=inverse_add(span,ms[0],Q(1))
                if span:span=span.affine(Q(-1))
        elif p['form']==2:
            span=inverse_add(span,ms[0 if negative else 1],Q(1,2))
            if span and negative:span=span.affine(Q(-1))
        else:
            base=0 if negative else 2
            span=inverse_op(span,ms[base+1],'mul',Q(1,2),True)
            span=inverse_add(span,ms[base],Q(1))if span else None
            if span and negative:span=span.affine(Q(-1))
        return span
    alive={p['id']:p for p in profiles};certificates=[];survivor_constraints=[]
    for order,pair in enumerate(pairs):
        for cid,p in list(alive.items()):
            negative,positive=allowed(pair['negative'],p,True),allowed(pair['positive'],p,False)
            result=negative.intersect(positive)if negative else None
            if result is not None:
                result=result.intersect(Interval(Q(0),Q(1,2)if p['form']==2 else Q(1)))
            if result is None or E64.preimage(result)is None:
                certificates.append(dict(profile=p,order=order,pair=pair,reason='empty_real_intersection'if result is None else'empty_ext80_lattice',remaining_interval=result.json()if result else None))
                del alive[cid]
    for pair in pairs:
        for p in alive.values():
            a,b=allowed(pair['negative'],p,True),allowed(pair['positive'],p,False)
            i=a.intersect(b).intersect(Interval(Q(0),Q(1,2)if p['form']==2 else Q(1)))
            survivor_constraints.append(dict(profile_id=p['id'],positive_trace_index=pair['positive']['trace_index'],negative_trace_index=pair['negative']['trace_index'],interval=i.json()))
    (output/'profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
    (output/'pairs.json').write_text(json.dumps(pairs,indent=2)+'\n')
    (output/'certificates.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in certificates))
    (output/'survivors.json').write_text(json.dumps(list(alive.values()),indent=2)+'\n')
    (output/'survivor-constraints.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in survivor_constraints))
    result=dict(profiles=660,excluded=len(certificates),surviving_constraint_profiles=len(alive),pairs=len(pairs),counts=counts,complete=True,
        seconds=time.monotonic()-start,input_sha256=hashlib.sha256(Path(trace).read_bytes()).hexdigest(),generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        forms=['negative=s,positive=1-s','negative=1-s,positive=s','negative=.5-s,positive=.5+s','negative=.5*(1-s),positive=.5*(1+s)'],
        scope='One arbitrary common finiteExt80 latent per exact |z|,reading,source,flag pair; each indicated operation ranges over five arithmetic models. Ordinary interior CDF targets only. Pair sharing is a hypothesis, not an Excel fact.',
        claim='Every excluded sign/reconstruction relation has an exact original paired contradiction. Surviving intervals are necessary shared-latent constraints, not complete programs, fitted values or identified helpers.',constant_fitted=False,new_excel_answers=0)
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def complement_prediction_sets(target_span):
    """Exact finite outputs over every Ext80 latent publishing into target_span."""
    import struct
    def bits(x):return '0x'+struct.pack('>d',float(x)).hex()
    span=B64.preimage(target_span)
    lo,hi=E64.ceil(span.lo),E64.floor(span.hi)
    if lo==span.lo and not span.lc:lo=E64.next(lo)
    if hi==span.hi and not span.hc:hi=E64.prev(hi)
    def complement(s,m):
        if m==0:return B64.rn(1-B64.rn(s))
        return B64.rn((E64 if m in(1,2)else E53).rn(1-s))
    sets=[]
    for model in range(5):
        low,high=complement(hi,model),complement(lo,model)
        a,b=int(bits(low),16),int(bits(high),16)
        assert 0<=b-a<=8 # exact monotone endpoint range, never an approximate search window
        allowed=[]
        for raw in range(a,b+1):
            t=Q(number(f'{raw:016x}'));p=B64.preimage(Interval(t,t))
            p=inverse_add(p,model,Q(1));p=p.affine(Q(-1))if p else None
            possible=span.intersect(p)
            if possible and E64.preimage(possible)is not None:allowed.append(f'0x{raw:016x}')
        assert allowed
        sets.append(tuple(allowed))
    return sets


def normal_cdf_questions(trace, symmetry, packet):
    """Jointly distinguish surviving relations; predictions are conditional sets."""
    import hashlib,datetime,itertools,struct,math
    packet=Path(packet);packet.mkdir(exist_ok=True)
    assert not (packet/'requests-NORM.S.DIST.json').exists()
    profiles=json.loads((Path(symmetry)/'survivors.json').read_text())
    assert profiles==[dict(id=i,form=0,models=[i])for i in range(5)]
    def bits(x):return '0x'+struct.pack('>d',float(x)).hex()
    def key(args):return tuple(f'0x{int(a,16):016x}'if isinstance(a,str)and a.startswith('0x')else a for a in args)
    def predictions(target):
        y=Q(number(target))
        return complement_prediction_sets(Interval(y,y))
    rows=list(trace_records(trace))
    known={key(r['args'])for r in rows};pool=[]
    cache={}
    for index,r in enumerate(rows):
        a=r['args']
        if len(a)!=2 or a[1]!='0x3ff0000000000000' or int(a[0],16)>>63!=1:continue
        if not r['expected'].startswith('0x')or not Q(0)<Q(number(r['expected']))<Q(1,2):continue
        if r['reading']['excel_build']not in ('20131','20430'):continue
        positive=[f'0x{int(a[0],16)&((1<<63)-1):016x}',a[1]]
        if r['expected']not in cache:cache[r['expected']]=predictions(r['expected'])
        sets=cache[r['expected']]
        pool.append(dict(trace_index=index,raw=r,positive_args=positive,predictions=sets,
                         positive_novel=key(positive)not in known))
    separated=set();selected=[];reasons=[]
    pairs=list(itertools.combinations(range(5),2))
    def score(r):return {pair for pair in pairs if set(r['predictions'][pair[0]]).isdisjoint(r['predictions'][pair[1]])}-separated
    # Prefer novel positive counterparts; maximize additional disjoint-model pairs.
    for _ in range(4):
        options=[r for r in pool if r not in selected and r['positive_novel']]
        best=min(options,key=lambda r:(-len(score(r)),r['trace_index']))
        gain=score(best)
        if not gain:break
        selected.append(best);reasons.append('additional disjoint prediction-set separation');separated|=gain
    # Published coordinates supply controls; no boundary is inferred from residuals.
    for target in (0.66291,math.sqrt(32),7.07106781186547,1.0):
        options=[r for r in pool if r not in selected]
        best=min(options,key=lambda r:(abs(abs(number(r['raw']['args'][0]))-target),r['trace_index']))
        selected.append(best);reasons.append(f'published/mathematical coordinate control near |z|={target!r}; not a claimed Excel branch')
    for build in ('20430','20131'):
        if len(selected)>=8:break
        best=min((r for r in pool if r not in selected and r['raw']['reading']['excel_build']==build),key=lambda r:r['trace_index'])
        selected.append(best);reasons.append('build/source transport and stability control')
    while len(selected)<8:
        best=next(r for r in pool if r not in selected);selected.append(best);reasons.append('fixed input-order control')
    requests=[];predicted=[];selection=[]
    for item,reason in zip(selected,reasons):
        old=item['raw'];negative_id=f'r10-norm-{len(requests):02d}';positive_id=f'r10-norm-{len(requests)+1:02d}'
        requests.extend([dict(id=negative_id,args=old['args'],input_transport='cell_value2_bulk'),
                         dict(id=positive_id,args=item['positive_args'],input_transport='cell_value2_bulk')])
        predicted.append(dict(id=negative_id,role='historical anchor recheck',predicted=old['expected'],historical_observation=old,
            condition='same admitted typed inputs and build-invariant old negative-tail result; not yet confirmed on20527'))
        predicted.append(dict(id=positive_id,role='conditional sign/reconstruction discriminator',
            condition=dict(negative_request=negative_id,negative_outcome=old['expected'],input_kept=True,
                           same_target_profile=True,latent_sharing='one Ext80 value for this exact |z| pair'),
            alternatives=[dict(profile_id=p['id'],possible_outcomes=list(values))for p,values in zip(profiles,item['predictions'])],
            not_a_complete_program_prediction=True))
        selection.append(dict(negative_id=negative_id,positive_id=positive_id,trace_index=item['trace_index'],
            historical_id=old['id'],source=old['source'],reading=old['reading'],reason=reason,
            positive_typed_input_novel=item['positive_novel'],argument_types=['binary64','binary64']))
    assert len(requests)==16 and len({key(r['args'])for r in requests})==16
    (packet/'requests-NORM.S.DIST.json').write_text(json.dumps(dict(function='NORM.S.DIST',witnesses=requests),indent=2)+'\n')
    (packet/'NORM-predictions.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in predicted))
    (packet/'NORM-selection.json').write_text(json.dumps(dict(pool=len(pool),selected_pairs=8,questions=16,
        additional_disjoint_model_pairs=sorted(separated),control_and_discriminator_selection=selection,
        surviving_profiles=profiles,method='greedy additional disjoint set separation, then independent published-coordinate and build/source controls; no numerical closeness score'),indent=2)+'\n')
    (packet/'README-NORM.S.DIST.txt').write_text('''R10 normal-CDF questions:16 typed Value2 calls in8 sign pairs, zero answers.
 cargo xtask witness read --func NORM.S.DIST --inputs <packet>/requests-NORM.S.DIST.json
Do not use --new-only: the8negative calls intentionally re-read old observations.
NORM-predictions.jsonl is separate and is never submitted as Excel evidence.
The5survivors are sign/reconstruction RELATIONS with an arbitrary latent tail,
not complete CDF programs: all declared concrete Cody/Hart kernels were refuted.
Positive predictions are exact finite SETS conditional on the paired negative
reproducing its old outcome on the current build and on admitted inputs/profile.
If that prerequisite fails, mark the positive prediction inapplicable and keep
the build/transport failure; do not silently substitute a returned anchor value.
Read both signs in one session and retain requested/submitted/read-back inputs,
raw outcomes and environment. Actual build/CPU/settings must be recorded.
These are fresh discovery questions, not validation of an identified helper.
''')
    files=['requests-NORM.S.DIST.json','NORM-predictions.jsonl','NORM-selection.json','README-NORM.S.DIST.txt']
    manifest=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),function='NORM.S.DIST',questions=16,
        historical_rechecks=8,positive_companions=8,novel_positive_inputs=sum(r['positive_novel']for r in selected),new_excel_answers=0,
        source_trace_sha256=hashlib.sha256(Path(trace).read_bytes()).hexdigest(),symmetry_survivors_sha256=hashlib.sha256((Path(symmetry)/'survivors.json').read_bytes()).hexdigest(),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),files_sha256={n:hashlib.sha256((packet/n).read_bytes()).hexdigest()for n in files},
        claim='Frozen conditional set predictions for5surviving reconstruction relations, not complete numerical programs or independently validated Excel internals.')
    (packet/'NORM-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))


def normal_cdf_lattice_questions(packet):
    """Construct probability-lattice probes from mathematics, not Excel fits."""
    import datetime,hashlib,struct
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'checker'))
    from reference import normal_bounds,LIB
    packet=Path(packet)
    assert not (packet/'requests-NORM.S.DIST-lattice.json').exists()
    p0,w=Q(1,1<<54),Q(1,1<<66)
    intervals=[Interval(p0-w,p0,True,False),Interval(p0,p0),Interval(p0,p0+w,False,True)]
    cases=[dict(negative_outcome_interval=i.json(),alternatives=[dict(profile_id=j,possible_outcomes=list(v))for j,v in enumerate(complement_prediction_sets(i))])for i in intervals]
    def raw(x):return int.from_bytes(struct.pack('>d',x),'big')
    def propose(target):
        low,high=raw(8.0),raw(9.0)
        while high-low>1:
            middle=(low+high)//2
            ref=normal_bounds((1<<63)|middle)
            a,b=Q(*map(int,ref['lower'])),Q(*map(int,ref['upper']))
            if a>target:low=middle
            elif b<target:high=middle
            else:raise ValueError('mathematical bisection precision insufficient')
        def distance(bits):
            ref=normal_bounds((1<<63)|bits)
            return abs((Q(*map(int,ref['lower']))+Q(*map(int,ref['upper'])))/2-target)
        return min((low,high),key=lambda b:(distance(b),b)),[low,high]
    root,bracket=propose(p0)
    robust,robust_bracket=propose(p0+Q(1,1<<67))
    requests=[];links=[]
    for bits,reason in [(root-1,'adjacent input below mathematical tail midpoint'),(root,'nearest input to mathematical tail midpoint'),
                        (root+1,'adjacent input above mathematical tail midpoint'),(robust,'interior of the exact double-rounding discriminator window')]:
        negative,positive=f'r10-norm-lattice-{len(requests):02d}',f'r10-norm-lattice-{len(requests)+1:02d}'
        a,b=f'0x{(1<<63)|bits:016x}',f'0x{bits:016x}'
        requests.extend([dict(id=negative,args=[a,'0x3ff0000000000000'],input_transport='cell_value2_bulk'),
                         dict(id=positive,args=[b,'0x3ff0000000000000'],input_transport='cell_value2_bulk')])
        links.append(dict(negative_request=negative,positive_request=positive,reason=reason,
            mathematical_reference_only=normal_bounds(a),
            prerequisite='admitted inputs; returned negative result lies in one of the3frozen intervals; same-profile shared latent for this sign pair',
            positive_prediction='use only the predeclared case table; outside its domain is a failed prerequisite, never a revised prediction'))
    assert len({tuple(q['args'])for q in requests})==8
    (packet/'requests-NORM.S.DIST-lattice.json').write_text(json.dumps(dict(function='NORM.S.DIST',witnesses=requests),indent=2)+'\n')
    predictions=dict(cases=cases,pairs=links,conditional_relations_only=True,complete_function_programs=0,
        probability_midpoint=['1',str(1<<54)],window=['1',str(1<<66)],
        mathematical_proposal=dict(mpfr_version=LIB.mpfr_get_version().decode(),root_bracket=bracket,robust_bracket=robust_bracket,
            targets=[['1',str(1<<54)],[str((p0+Q(1,1<<67)).numerator),str((p0+Q(1,1<<67)).denominator)]],
            source='DLMF7.20.1 normal CDF with directed MPFR enclosures; choose inputs, not coefficients or target answers'),
        explanation='For negative publication p just above2^-54, strict and PC53 complement predict prev(1), while PC64 then binary64 predicts1. The returned negative determines a frozen case, not an adapted helper value.')
    (packet/'NORM-lattice-predictions.json').write_text(json.dumps(predictions,indent=2)+'\n')
    base=json.loads((packet/'requests-NORM.S.DIST.json').read_text())['witnesses']
    combined=base+requests
    assert len(combined)==len({tuple(q['args'])for q in combined})==24
    (packet/'requests-NORM.S.DIST-all.json').write_text(json.dumps(dict(function='NORM.S.DIST',witnesses=combined),indent=2)+'\n')
    (packet/'RUN-ROUND10.txt').write_text('''Run exactly these two batches with the updated DnaCalc master on Windows:
 cargo xtask witness read --func PMT --inputs <packet>/requests-PMT.json
 cargo xtask witness read --func NORM.S.DIST --inputs <packet>/requests-NORM.S.DIST-all.json
64PMT calls plus24normal-CDF calls =88 total, zero supplied answers.
Do not also run the two smaller NORM request files: the all file includes them.
Do not use --new-only: there are deliberate historical negative-tail rechecks.
Retain raw worker request/result files, frozen model predictions, readings and audits.
The8lattice calls are4new sign pairs from mathematical probability-rounding cells.
NORM-lattice-predictions.json contains a frozen3-case relation table: classify only
by the returned negative probability in its predeclared exact intervals. A return
outside all cases is an unmet prerequisite; do not revise the table after answers.
No mathematical-reference value is an Excel answer or an identified private helper.
Other NORM predictions keep their exact historical-anchor prerequisites unchanged.
Literal PMT readings remain counted transport audits, never ordinary semantic rows.
Record actual build/CPU/settings and requested/submitted/read-back inputs. No current
build20527 reproduction or complete numerical-program validation is claimed yet.
''')
    files=['requests-NORM.S.DIST-lattice.json','NORM-lattice-predictions.json','requests-NORM.S.DIST-all.json','RUN-ROUND10.txt']
    manifest=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),new_lattice_calls=8,total_normal_calls=24,total_PMT_calls=64,
        new_excel_answers=0,coefficient_fitted=False,profile_hypotheses=5,
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        files_sha256={n:hashlib.sha256((packet/n).read_bytes()).hexdigest()for n in files},
        parent_normal_request_sha256=hashlib.sha256((packet/'requests-NORM.S.DIST.json').read_bytes()).hexdigest())
    (packet/'NORM-lattice-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest))


if len(sys.argv)>1 and sys.argv[1]=='normal-cdf-lattice-questions':
    normal_cdf_lattice_questions(*sys.argv[2:])
    sys.exit(0)

if len(sys.argv)>1 and sys.argv[1]=='normal-cdf-questions':
    normal_cdf_questions(*sys.argv[2:])
    sys.exit(0)

if len(sys.argv)>1 and sys.argv[1]=='normal-cdf-symmetry':
    normal_cdf_symmetry(*sys.argv[2:])
    sys.exit(0)

if len(sys.argv)>1 and sys.argv[1]=='normal-cdf-linear':
    normal_cdf_linear(*sys.argv[2:])
    sys.exit(0)

root=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(exist_ok=False)
meta=json.loads((root/'summary.json').read_text());programs=meta['numerator_programs']*3
alive={i:None for i in range(programs)};bound_rows={};killed=set();stats=[];counts=collections.Counter()
with (out/'first-kills.jsonl').open('x')as f:
    for line in (root/'matrix.jsonl').open():
        row=json.loads(line);order=row['order'];before=len(alive)
        if alive:
            expected=row['witness']['expected']
            assert expected.startswith('0x')and number(expected)>0
            target=Q(number(expected));target=B64.preimage(Interval(target,target))
            for group in row['values']:
                n=exact_root(group['value']);assert n is not None and n>0
                for model in range(3):
                    required=inverse_op(target,model,'div',n,False)
                    for numerator in group['numerator_ids']:
                        cid=numerator*3+model
                        if cid not in alive:continue
                        old=alive[cid]
                        new=required if old is None else old.intersect(required)
                        if new is None or E64.preimage(new)is None:
                            assert cid not in killed;killed.add(cid)
                            reason='empty_real_intersection'if new is None else'empty_ext80_lattice'
                            counts[reason]+=1
                            record=dict(candidate_id=cid,numerator_id=numerator,division_model=model,order=order,
                                trace_index=row['trace_index'],expected=expected,numerator=group['value'],
                                reason=reason,previous_interval=old.json()if old else None,
                                required_interval=required.json()if required else None,
                                remaining_interval=new.json()if new else None,
                                previous_bound_rows=bound_rows.get(cid),new_row=order)
                            f.write(json.dumps(record,separators=(',',':'))+'\n');del alive[cid];bound_rows.pop(cid,None)
                        else:
                            if old is None:bound_rows[cid]=[order,order]
                            else:
                                if new.lo!=old.lo or new.lc!=old.lc:bound_rows[cid][0]=order
                                if new.hi!=old.hi or new.hc!=old.hc:bound_rows[cid][1]=order
                            alive[cid]=new
        stat=dict(order=order,trace_index=row['trace_index'],id=row['witness']['id'],killed=before-len(alive),remaining=len(alive));stats.append(stat)
        if stat['killed']:f.flush();print(json.dumps(stat),flush=True)
with(out/'survivor-constraints.jsonl').open('x')as f:
    for cid,i in alive.items():f.write(json.dumps(dict(candidate_id=cid,interval=i.json()if i else None,bound_rows=bound_rows.get(cid)))+'\n')
with(out/'witness-kills.jsonl').open('x')as f:
    for s in stats:f.write(json.dumps(s)+'\n')
summary=dict(program_profiles=programs,profiles_killed_once=len(killed),surviving_constraint_profiles=len(alive),
    rows=len(stats),available_rows=meta['available_build20131_rows'],counts=counts,
    all_profiles_eliminated=not alive,complete_for_exported_rows=True,
    claim='Each killed numerator/division profile has no single count-only denominator shared across its recorded means. A surviving interval is a necessary constraint, not a selected constant or identified program.',
    scope='Original three arithmetic models, normalization exp(-(S+d))/G, build20131, count1; first64 original scoped rows. Other normalizations, five-model refinements, and different numerator methods are outside this projection.',
    constant_fitted=False,new_excel_answers=0)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
