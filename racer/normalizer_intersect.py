"""Exact unknown-normalizer constraints. Never select a denominator value."""
from pathlib import Path
import collections
import json
import sys
from fractions import Fraction as Q
from intervals import B64,E64,Interval,inverse_op,number
from certify import exact_root

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
