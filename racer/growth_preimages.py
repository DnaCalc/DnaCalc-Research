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
from intervals import B64, E64, Interval, inverse_op, number
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


if __name__=='__main__':
    if sys.argv[1]=='compare':compare(Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]))
    else:run(Path(sys.argv[1]),Path(sys.argv[2]))
