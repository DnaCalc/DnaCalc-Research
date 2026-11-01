"""Compare raw native instruction measurements; never ingest Excel evidence."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path


def table(path, result):
    with Path(path).open(newline='')as file:
        reader=csv.DictReader(file,delimiter='\t')
        expected=['id','kind','cw','a','b']+(['result']if result else[])
        if reader.fieldnames!=expected:raise ValueError('changed native probe header')
        rows=list(reader)
    if len({r['id']for r in rows})!=len(rows):raise ValueError('duplicate native leaf ID')
    for row in rows:
        if row['kind']not in ('fyl2x','fyl2xp1','f2xm1','fldln2','fldl2e'):raise ValueError('unknown instruction')
        if row['cw']not in ('123f','133f'):raise ValueError('unknown control word')
        for key in ('a','b')+ (('result',)if result else()):
            if len(row[key])!=20 or len(bytes.fromhex(row[key]))!=10:raise ValueError('not an 80-bit encoding')
    return rows


def compare(requests,reference,actual,output):
    paths=[Path(p)for p in (requests,reference,actual)]
    frozen={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}
    requested=table(requests,False);baseline=table(reference,True);measured=table(actual,True)
    strip=lambda rows:[{k:v for k,v in r.items()if k!='result'}for r in rows]
    if requested!=strip(baseline)or requested!=strip(measured):
        raise ValueError('missing, reordered, renamed or changed native instruction inputs')
    differences=[dict(id=a['id'],kind=a['kind'],cw=a['cw'],a=a['a'],b=a['b'],
                     reference=r['result'],actual=a['result'])for r,a in zip(baseline,measured)if r['result']!=a['result']]
    result=dict(complete=True,instructions=len(requested),matches=len(requested)-len(differences),
        mismatches=len(differences),mismatches_by_kind=dict(Counter(r['kind']for r in differences)),
        differences=differences,all_certificate_leaves_reproduced=not differences,input_sha256=frozen,
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limits='Native CPU instruction comparison only. Equal leaves support transfer of this finite certificate on these operands; no general CPU equivalence or Excel-function match. Changed leaves require a new run and new certificates on that CPU, never replacement of Excel answers.')
    if frozen!={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths}:raise ValueError('probe files changed')
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('requests','reference','actual','output'):parser.add_argument(key)
    a=parser.parse_args();r=compare(a.requests,a.reference,a.actual,a.output)
    print(json.dumps({k:v for k,v in r.items()if k!='differences'}))
    raise SystemExit(1 if r['mismatches'] else 0)
