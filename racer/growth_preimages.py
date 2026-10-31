"""Exact growth-factor constraints isolated by zero-payment PV/FV observations.

Tier 1 annuity algebra and IEEE 754. No rounded inversion or fitted helper.
The result constrains each explicit final arithmetic model separately.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import sys
from intervals import B64, E64, E53, Interval, inverse_op, number, stored
from certify import exact_root


def key(raw):
    return tuple(raw['args'][:2])+(json.dumps(raw['reading'],sort_keys=True),)


def preimage(case,model):
    raw=case['witness']
    a=list(map(number,raw['args']))
    y=number(raw['expected'])
    if not all(math.isfinite(x)for x in a+[y])or not y or not a[3]:
        return None
    target=Interval(Q(y),Q(y))
    # Resident arithmetic still publishes through a final binary64 store.
    # Model1 includes that store inside inverse_op; model2 needs it here.
    if model==2:
        target=B64.preimage(target)
    # FV(r,n,0,pv)=-pv*g; PV(r,n,0,fv)=-fv/g.
    i=inverse_op(target,model,'mul' if raw['function']=='FV'else'div',Q(-a[3]),False)
    return i.intersect(Interval(Q(0),E64.maximum,False,True))if i else None


def run(input_path,out):
    out.mkdir(parents=True,exist_ok=False)
    cases=[json.loads(l)for l in open(input_path)]
    groups=defaultdict(list)
    excluded=[]
    records=[]
    for idx,c in enumerate(cases):
        raw=c['witness'];a=list(map(number,raw['args']))
        if a[0]<=-1 or not math.isfinite(a[0])or not math.isfinite(a[1])or not number(raw['expected']):
            excluded.append(dict(case=idx,reason='outside positive finite nonzero growth isolation'));continue
        models=[]
        for model in range(3):
            i=preimage(c,model)
            models.append(dict(model=model,interval=i.json()if i else None))
        rec=dict(case=idx,witness=raw,models=models)
        records.append(rec);groups[key(raw)].append((idx,c))
    with (out/'row-preimages.jsonl').open('x')as f:
        for rec in records:f.write(json.dumps(rec,separators=(',',':'))+'\n')
    intersections=[]
    for k,cs in sorted(groups.items()):
        if len({c['witness']['function']for _,c in cs})<2:
            continue
        for pm in range(3):
            for fm in range(3):
                current=None;first=True;bounds=[]
                for idx,c in cs:
                    i=preimage(c,pm if c['witness']['function']=='PV'else fm)
                    nxt=i if first else current.intersect(i)if current is not None and i else None
                    if first or nxt!=current:bounds.append(idx)
                    first=False;current=nxt
                    if current is None:break
                intersections.append(dict(rate=k[0],nper=k[1],environment=json.loads(k[2]),pv_model=pm,fv_model=fm,
                    row_count=len(cs),real_empty=current is None,ext80_value_exists=E64.preimage(current)is not None if current else False,
                    interval=current.json()if current else None,bound_cases=bounds))
    (out/'joint-intersections.jsonl').write_text(''.join(json.dumps(r,separators=(',',':'))+'\n'for r in intersections))
    summary=dict(cases=len(cases),processed=len(records),excluded=excluded,
        shared_pv_fv_groups=len({(x['rate'],x['nper'],json.dumps(x['environment'],sort_keys=True))for x in intersections}),
        model_pair_intersections=len(intersections),empty_intersections=sum(x['real_empty']for x in intersections),
        input_sha256=hashlib.sha256(input_path.read_bytes()).hexdigest(),
        precision='Exact Fraction intervals, nearest-even ties, open endpoints, strict input stores and PC64 double rounding.',
        limitation='Growth-only observations isolate a factor conditional on one final multiply/divide. Other cashflows, domain branches and a more complex final program are not constrained by this projection.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


def compare(preimages,diagnostics,out):
    out.mkdir(parents=True,exist_ok=False)
    constraints={x['case']:x for x in map(json.loads,open(preimages/'row-preimages.jsonl'))}
    counts=defaultdict(Counter);kills={};group_sets=defaultdict(set)
    def from_json(x):
        return Interval(Q(*map(int,x['lower'])),Q(*map(int,x['upper'])),x['lower_closed'],x['upper_closed'])if x else None
    with (out/'diagnostic-constraints.jsonl').open('x')as f:
        for row in map(json.loads,open(diagnostics)):
            if row['case']not in constraints:continue
            c=constraints[row['case']];raw=c['witness']
            memberships=[]
            for candidate in row['candidates']:
                name=candidate['name'];q=exact_root(candidate['value'])
                for m in c['models']:
                    i=from_json(m['interval']);hit=q is not None and i is not None and i.contains(q)
                    k=(name,raw['function'],m['model'])
                    counts[k]['observations']+=1;counts[k]['within_preimage']+=hit
                    if not hit and k not in kills:kills[k]=dict(case=row['case'],witness=raw,value=candidate['value'],interval=m['interval'])
                    if hit:group_sets[(raw['function'],m['model'])].add(name)
                    memberships.append(dict(helper=name,model=m['model'],contains=hit))
            f.write(json.dumps(dict(case=row['case'],memberships=memberships),separators=(',',':'))+'\n')
    result=dict(conditional_profiles=[dict(helper=n,function=fn,model=m,**c,first_contradiction=kills.get((n,fn,m)))for(n,fn,m),c in sorted(counts.items())],
        claim='Named growth diagnostics tested against exact conditional preimages. Applicability limits are explicit; a compatible finite sample is not algorithm identification.')
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(profiles=len(counts),compatible_on_all_applicable=sum(c['observations']==c['within_preimage']for c in counts.values()))))


def loan_prepare(traces,out):
    """Preserve original observations; derived math queries contain no answers."""
    import struct
    out.mkdir(parents=True,exist_ok=False)
    pp=list(map(json.loads,open(traces/'PPMT.jsonl')))
    ip=list(map(json.loads,open(traces/'IPMT.jsonl')))
    loans=[];queries=[]
    def env(r):return json.dumps(r['reading'],sort_keys=True)
    wanted={(tuple([r['args'][0]]+r['args'][2:]),env(r))for r in pp}
    pm={k:[]for k in wanted}
    for row in map(json.loads,open(traces/'PMT.jsonl')):
        k=(tuple(row['args']),env(row))
        if k in pm:pm[k].append(row)
    for r in pp:
        a=list(map(number,r['args']))
        assert len(a)==6 and a[0]>0 and a[4]==a[5]==0
        assert a[1].is_integer()and a[2].is_integer()and 1<=a[1]<=a[2]
        k=(tuple([r['args'][0]]+r['args'][2:]),env(r))
        p=pm[k];i=[s for s in ip if s['args']==r['args']and env(s)==env(r)]
        assert p and i
        remaining=a[2]-a[1]+1
        assert remaining.is_integer()and abs(remaining)<2**53
        queries.append(dict(witness=dict(args=[r['args'][0],'0x'+struct.pack('>d',-remaining).hex()]),
            kind='derived mathematical growth query; not an Excel observation'))
        loans.append(dict(ppmt=r,ipmt=i,pmt=p,remaining_periods=remaining))
    (out/'observations.json').write_text(json.dumps(loans,indent=2)+'\n')
    (out/'math-queries.jsonl').write_text(''.join(json.dumps(x,separators=(',',':'))+'\n'for x in queries))
    print(json.dumps(dict(ppmt_rows=len(loans),pmt_loans=len(wanted),queries=len(queries))))


def loan_compare(prepared,diagnostics,out):
    """Constrain one latent PMT with all observed publications, never fit PMT."""
    out.mkdir(parents=True,exist_ok=False)
    rows=json.loads((prepared/'observations.json').read_text())
    ds=list(map(json.loads,open(diagnostics)))
    assert len(rows)==len(ds)
    names=set.intersection(*[{c['name']for c in d['candidates']}for d in ds])
    assert len({tuple(r['pmt'][0]['args'])for r in rows})==1
    initial=None
    for row in rows:
        for p in row['pmt']:
            y=Q(number(p['expected']));i=B64.preimage(Interval(y,y))
            initial=i if initial is None else initial.intersect(i)
    assert initial is not None
    def multiply_preimage(raw,m,c):
        y=Q(number(raw['expected']))
        return inverse_op(B64.preimage(Interval(y,y)),m,'mul',c,True)
    def lattice(i,fmt):
        x=fmt.ceil(i.lo)if i else None
        if x is not None and x==i.lo and not i.lc:x=fmt.next(x)
        return x is not None and i.contains(x)
    def rnd(x,m):
        if m==0:return B64.rn(x)
        if m==1:return stored(x)
        if m==2:return E64.rn(x)
        y=E53.rn(x)
        return B64.rn(y)if m==3 and y is not None else y
    def scalar(d,m):
        if m==0:d=B64.rn(d)
        return rnd(1-d,m)
    profiles=[]
    for name in sorted(names):
        values=[exact_root(next(c['value']for c in d['candidates']if c['name']==name))for d in ds]
        assert all(x is not None and x>0 for x in values)
        for ppmodel in range(5):
            current=initial;bound=[];firstkill=None
            for idx,(row,d)in enumerate(zip(rows,values)):
                i=multiply_preimage(row['ppmt'],ppmodel,d)
                nxt=current.intersect(i)if current is not None and i else None
                if nxt!=current:bound.append(idx)
                current=nxt
                if current is None:firstkill=idx;break
            record=dict(helper=name,ppmt_model=ppmodel,ppmt_interval=current.json()if current else None,
                ppmt_first_contradiction=firstkill,bound_rows=bound,
                ppmt_ext80_pmt_exists=lattice(current,E64),ppmt_binary64_pmt_exists=lattice(current,B64))
            ips=[]
            if current:
                for sm in range(5):
                    for im in range(5):
                        both=current;kill=None
                        for idx,(row,d)in enumerate(zip(rows,values)):
                            c=scalar(d,sm)
                            assert c is not None and c!=0
                            for ip in row['ipmt']:
                                i=multiply_preimage(ip,im,c)
                                both=both.intersect(i)if both is not None and i else None
                            if both is None:kill=idx;break
                        ips.append(dict(subtract_model=sm,ipmt_model=im,first_contradiction=kill,
                            interval=both.json()if both else None,ext80_pmt_exists=lattice(both,E64)))
            record['joint_ipmt_profiles']=ips;profiles.append(record)
    (out/'profiles.json').write_text(json.dumps(profiles,indent=2)+'\n')
    summary=dict(loans=1,ppmt_rows=len(rows),ipmt_rows=sum(len(x['ipmt'])for x in rows),
        growth_providers=len(names),ppmt_profiles=len(profiles),
        ppmt_real_compatible=sum(x['ppmt_interval']is not None for x in profiles),
        ppmt_ext80_compatible=sum(x['ppmt_ext80_pmt_exists']for x in profiles),
        joint_profiles_tested=sum(len(x['joint_ipmt_profiles'])for x in profiles),
        joint_ext80_compatible=sum(y['ext80_pmt_exists']for x in profiles for y in x['joint_ipmt_profiles']),
        grammar='D=growth(rate,-(nper-period+1)); PPMT=P*D; IPMT=P*(1-D); each arithmetic node has a fixed one-of-five nearest/even model; one latent P shared across every period and constrained by original PMT publications.',
        prepared_sha256=hashlib.sha256((prepared/'observations.json').read_bytes()).hexdigest(),
        diagnostics_sha256=hashlib.sha256(diagnostics.read_bytes()).hexdigest(),
        limitation='Conditional helper/root projection for one loan, not a PMT implementation or whole-family identification. No observed PMT value is installed as a program constant.')
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))


if __name__=='__main__':
    if sys.argv[1]=='loan-prepare':loan_prepare(Path(sys.argv[2]),Path(sys.argv[3]))
    elif sys.argv[1]=='loan-compare':loan_compare(Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]))
    elif sys.argv[1]=='compare':compare(Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]))
    else:run(Path(sys.argv[1]),Path(sys.argv[2]))
