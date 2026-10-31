"""Complete first-kill elimination of fixed combines with an unknown helper.

One arbitrary binary64/Ext80 expm1 value is shared by every row of a fixed
rate/nper/environment group. Once a combine has an exact contradiction it is
never rescued. Later groups kill zero. This is a stronger constraint than any
particular expm1 implementation, not a fitted or row-specific helper model.
"""
import collections
import hashlib
import json
from pathlib import Path
import sys
import time
import math
import intersect as engine
from intervals import E64, Interval, number


def verify(record, profile, rows):
    killed=record['first_empty']
    if killed:
        ids=sorted(set(killed['bound_rows']+[killed['new_row']]))
    else:
        ids=record['bound_rows']
    assert ids
    cur=None
    first=True
    for idx in ids:
        i=engine.row_preimage(rows[idx],profile)
        cur=i if first else (cur.intersect(i) if cur is not None else None)
        first=False
    assert cur is None or E64.preimage(cur) is None,record
    return dict(reason='empty_real_intersection' if cur is None else 'empty_ext80_lattice',
        sufficient_row_indices=ids, remaining_interval=cur.json() if cur else None,
        witnesses=[dict(trace_index=rows[idx]['trace_index'],witness=rows[idx]['raw'])for idx in ids])


def run(input_path,catalogue,output,interest=False):
    start=time.monotonic()
    rows=[]
    groups=collections.defaultdict(list)
    selection_counts=collections.Counter()
    for line in open(input_path):
        r=json.loads(line)
        if number(f'{r["leaves"][2]:016x}')!=0:
            selection_counts['FV_nonzero_control']+=1;continue
        if interest:
            rate=number(f'{r["leaves"][0]:016x}')
            tau=number(r['raw']['args'][1])*number(f'{r["log1p"]:016x}')
            if not(0<rate<=0.25 and 0<tau<1):
                selection_counts['outside_interest_research_scope']+=1;continue
        idx=len(rows);rows.append(r);groups[engine.key(r)].append(idx)
        selection_counts['scoped']+=1
    engine.ROWS=rows
    if interest:
        precondition=json.loads(catalogue.read_text())
        allowed=precondition['allowed_model_indices']
        descriptor={k:precondition.get(k,default)for k,default in [('sum_kind','add'),('amortization_numerator','interest'),('amortization_operation','div')]}
        es=[dict(family='shared_interest',**descriptor,programme='I=PV*rate; amort=(I or PV) op helper; output=-(I+amort)/tf or -fma(PV,rate,amort)/tf; descriptor and all five site models are explicit')]
        names=precondition.get('arithmetic_models',['strict','stored_pc64','resident_pc64'])
        radix=len(names)
        cs=[]
        for code in range(radix**5):
            digits=[];n=code
            for _ in range(5):digits.append(n%radix);n//=radix
            if digits[0]not in allowed:continue
            if descriptor['sum_kind']=='fma'and digits[3]!=0:continue
            if descriptor['sum_kind']=='fma'and descriptor['amortization_numerator']=='pv'and digits[0]!=0:continue
            cs.append(dict(id=code,family='shared_interest',expression=0,models=[names[i]for i in digits],**descriptor))
    else:
        es=json.loads((catalogue/'expressions.json').read_text())
        cs=json.loads((catalogue/'combines.json').read_text())
    profiles=[]
    for c in cs:
        ms=[engine.MODELS[m] for m in c['models']]
        if interest:expr=es[0]
        else:
            expr,end=engine.compile_expr(es[c['expression']],ms,1)
            assert end==4
        profiles.append((c,expr,ms))
    output.mkdir(parents=True,exist_ok=True)
    result={}
    for cohort in ['all','nonunit','timing','nonunit_timing']:
        if 'timing' in cohort:
            selected_groups=collections.defaultdict(list)
            for idx,row in enumerate(rows):selected_groups[engine.key(row)+(row['leaves'][3],)].append(idx)
        else:
            selected_groups=groups
        alive={p[0]['id']:p for p in profiles}
        counts=collections.Counter()
        killers=[]
        real_kills=0
        lattice_kills=0
        seen_rows=0
        seen_groups=0
        with (output/f'{cohort}-first-kills.jsonl').open('x') as kills, (output/f'{cohort}-group-kills.jsonl').open('x') as stats:
            for gid,ids in enumerate(selected_groups.values()):
                row=rows[ids[0]]
                if cohort.startswith('nonunit') and number(row['raw']['args'][1])==1:
                    continue
                seen_groups+=1
                seen_rows+=len(ids)
                before=len(alive)
                if alive:
                    engine.PROFILES=list(alive.values())
                    _,records,subcounts=engine.process((gid,ids))
                    counts.update(subcounts)
                    for record in records:
                        if not record['ext80_value_exists'] and record['status']!='unsupported':
                            cid=record['combine']
                            proof=verify(record,alive[cid],rows)
                            real_kills+=proof['reason']=='empty_real_intersection'
                            lattice_kills+=proof['reason']=='empty_ext80_lattice'
                            certificate=dict(candidate_id=cid,group=gid,rate=f'0x{row["leaves"][0]:016x}',
                                nper=row['raw']['args'][1],environment=row['raw']['reading'],canonical_timing=row['leaves'][3] if 'timing' in cohort else None,
                                combine=alive[cid][0],expression=es[alive[cid][0]['expression']],
                                full_group_rows=len(ids),first_empty=record['first_empty'],proof=proof)
                            kills.write(json.dumps(certificate,separators=(',',':'))+'\n')
                            del alive[cid]
                stat=dict(group=gid,rate=f'0x{row["leaves"][0]:016x}',nper=row['raw']['args'][1],
                    environment=row['raw']['reading'],rows=len(ids),killed=before-len(alive),remaining=len(alive))
                stats.write(json.dumps(stat,separators=(',',':'))+'\n')
                if before!=len(alive):killers.append(stat);kills.flush();stats.flush();print(cohort,json.dumps(stat),flush=True)
        with (output/f'{cohort}-survivors.json').open('x') as f:json.dump(list(alive),f)
        result[cohort]=dict(profiles=len(profiles),survivors=len(alive),groups=seen_groups,rows=seen_rows,
            real_intersection_kills=real_kills,ext80_lattice_kills=lattice_kills,
            first_kill_statistics=killers,computed_group_profile_counts=counts,complete=True)
    summary=dict(cohorts=result,seconds=time.monotonic()-start,
        family='shared_interest'if interest else'original_four_leaf_combine',selection_counts=selection_counts,
        programme_descriptor=descriptor if interest else None,
        input_sha256=hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
        claim='Each killed fixed combine profile has a rate/nper/environment group with no common finite binary64 or Ext80 expm1 value. Real-empty certificates also exclude arbitrary real em. Nonunit cohorts exclude a one-period-identity rescue; timing cohorts permit the helper to depend on canonical payment timing. Nonfinite/zero helper values cannot produce the finite nonzero killing targets.',
        limitation='This excludes only the killed specified combines/precision programs with a helper shared independently of PV (and independently of timing in cohorts without timing in their name). Surviving profiles would be constraints, not identified expm1 programs. It does not exclude PV-dependent arithmetic, additional operations/branches or other programs.',
        proof='Exact rational nearest/even preimages; each killed profile has a verified sufficient subset of at most three original rows. Later groups are zero-kill groups, never passes.')
    with (output/'summary.json').open('x') as f:json.dump(summary,f,indent=2)
    print(json.dumps(summary))


if __name__=='__main__':
    run(sys.argv[1],Path(sys.argv[2]),Path(sys.argv[3]),interest='--interest'in sys.argv[4:])
