"""Consolidate the 180 unanswered calls and freeze instruction-boundary audits.

Requests use the existing witness read schema. No result is put in a request.
Original packet bytes, IDs, prediction files and manifests remain unchanged.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

from x87_programs import Builder, Native, ROOT, WORDS, digest, write

ORIGINALS = (('round11',64),('structural-r11-v2',52),('round12-pmt-groups',64))
ANCHORS = ('r10-pmt-22','r10-pmt-20','r10-pmt-02','r10-pmt-06','r10-pmt-14','r10-pmt-17')


def standalone(function, args, reciprocal=False):
    g=Builder(len(args))
    if function=='LN':
        coefficient=g.op('fldln2',store=False)
        result=g.op('fyl2x',coefficient,0)
        return g.graph(result,result)
    if function=='EXP':
        coefficient=g.op('fldl2e',store=False)
        t=g.op('mul',0,coefficient)
    else:
        t=g.op('fyl2x',1,0)
    # Fixed signed or positive-then-reciprocal reconstruction; branch is solely
    # on the requested exponent sign, never on a returned Excel answer.
    exponent_negative = bool(int(args[0 if function=='EXP' else 1],16)>>63)
    if reciprocal and exponent_negative:
        t=g.op('neg',t,store=False)
    k=g.op('frndint',t,store=False)
    f=g.op('sub',t,k)
    z=g.op('f2xm1',f)
    m=g.op('add',g.one,z)
    result=g.op('fscale',m,k)
    if reciprocal and exponent_negative:
        result=g.op('div',g.one,result)
    return g.graph(result,result)


def create(trace,gaps,library,output,leaves=None,environment=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    native=Native(library)
    batches=[];reading=[];origin_counts=Counter();sources={}
    for name,calls in ORIGINALS:
        src=ROOT/'packets'/name;dst=output/'prior'/name
        manifest=json.loads((src/'manifest.json').read_text())
        for path,expected_digest in manifest['files_sha256'].items():
            if digest(src/path)!=expected_digest:raise ValueError('original packet changed: '+name+'/'+path)  # ubs:ignore - Public packet integrity hashes, never secrets or authentication.
        shutil.copytree(src,dst)
        found=0
        for path in sorted(dst.glob('requests-*.json')):
            doc=json.loads(path.read_text());f=doc['function'];found+=len(doc['witnesses'])
            batches.append(dict(function=f,path=path.relative_to(output).as_posix(),
                                calls=len(doc['witnesses']),origin=name))
            for row in doc['witnesses']:
                reading.append(dict(function=f,origin=name,**row))
                origin_counts[name]+=1
        if found!=calls:raise ValueError('original request coverage differs')
        sources[name]=dict(manifest_sha256=digest(src/'manifest.json'),calls=calls)
    rows=[json.loads(line)for line in Path(trace).read_text().splitlines()]
    byid={r['id']:r for r in rows}
    gap_ids=[json.loads(line)['id']for line in Path(gaps).read_text().splitlines()]
    requests={f:[]for f in ('PMT','POWER','LN','EXP')};selection=[];seen={f:{}for f in requests}
    def add(f,args,reason,source_id,derivation=None):
        key=tuple(args)
        if key not in seen[f]:
            row=dict(id=f'r13-{f.lower()}-{len(requests[f]):03d}',args=args,input_transport='cell_value2_bulk')
            requests[f].append(row);seen[f][key]=row['id']
        selection.append(dict(id=seen[f][key],function=f,source_id=source_id,
                              reason=reason,derivation=derivation))
    for id in gap_ids:
        add('PMT',byid[id]['args'],'Recheck every original final-product gap with retained Value2 read-back',id)
    for id in ANCHORS:
        row=byid[id];r,n=map(lambda x:int(x,16),row['args'][:2])
        base=Builder(1);b=base.op('add',base.one,0)
        base_bits=native.publish(native.evaluate(base.graph(b,b),[r],0x133f)[:10])
        base_arg=f'0x{base_bits:016x}'
        add('LN',[base_arg],'Observe the worksheet log boundary at the explicitly stored 1+r',id,
            dict(rule='FADD PC64 (1,r), FSTP binary64',rate=row['args'][0]))
        for exponent in (n,n^(1<<63)):
            add('POWER',[base_arg,f'0x{exponent:016x}'],
                'Compare worksheet POWER with fixed native reconstruction hypotheses; not an internal PMT trace',id,
                dict(rule='stored base and original signed periods',rate=row['args'][0]))
        for cw in WORDS:
            exponent=Builder(2);rom=exponent.op('fldln2',store=False)
            log=exponent.op('fyl2xp1',rom,0)
            t=exponent.op('mul',1,log)
            for sign in (n,n^(1<<63)):
                value=native.evaluate(exponent.graph(t,t),[r,sign],cw)[:10]
                bits=native.publish(value,cw)
                add('EXP',[f'0x{bits:016x}'],
                    'Observe EXP at a predetermined native exponent boundary; hypothetical helper input',id,
                    dict(rule='FYL2XP1(FLDLN2,r), FMUL signed n, FSTP binary64',
                         cw=cw,rate=row['args'][0],periods=f'0x{sign:016x}',resident_t=value.hex()))
    for f,rs in requests.items():
        write(output/f'requests-{f}.json',dict(function=f,
            selection='R13 instruction-boundary and original-gap audit; intentional rechecks, no --new-only',witnesses=rs))
        batches.append(dict(function=f,path=f'requests-{f}.json',calls=len(rs),origin='r13'))
        for row in rs:
            reading.append(dict(function=f,origin='r13',**row));origin_counts['r13']+=1
    (output/'predictions').mkdir()
    # Freeze every elementary template mask; none is a surviving PMT program.
    with (output/'predictions'/'native-elementary.jsonl').open('x')as out:
        for f in ('LN','EXP','POWER'):
            for row in requests[f]:
                for reciprocal in ((False,)if f=='LN'else(False,True)):
                    graph=standalone(f,row['args'],reciprocal)
                    for cw in WORDS:
                        values=native.masks(graph,[int(x,16)for x in row['args']],cw)
                        predictions=[f'0x{native.publish(values[i:i+10],cw):016x}'for i in range(0,len(values),20)]
                        out.write(json.dumps(dict(id=row['id'],function=f,reciprocal=reciprocal,
                            cw=cw,graph=graph,predictions_by_store_mask=predictions))+'\n')
    write(output/'selection-r13.json',selection)
    write(output/'batches.json',batches)
    with (output/'READING-LIST.tsv').open('x')as out:
        out.write('ordinal\torigin\tfunction\tid\ttransport\targuments_json\n')
        for i,r in enumerate(reading,1):
            out.write(f"{i}\t{r['origin']}\t{r['function']}\t{r['id']}\t{r.get('input_transport','cell_value2_bulk')}\t{json.dumps(r['args'],separators=(',',':'))}\n")
    (output/'RUN-R13.ps1').write_text('''param([Parameter(Mandatory=$true)][string]$DnaCalcCheckout)
$ErrorActionPreference = 'Stop'
$packet = $PSScriptRoot
Push-Location $DnaCalcCheckout
try {
    $batches = Get-Content -LiteralPath (Join-Path $packet 'batches.json') -Raw | ConvertFrom-Json
    foreach ($batch in $batches) {
        $inputFile = Join-Path $packet $batch.path
        Write-Host ("R13: {0} / {1}: {2} calls" -f $batch.origin, $batch.function, $batch.calls)
        & cargo xtask witness read --func $batch.function --inputs $inputFile
        if ($LASTEXITCODE -ne 0) { throw "Reading failed: $($batch.path). Preserve all raw worker output before continuing." }
    }
} finally { Pop-Location }
'''.replace('\n+','\n'))
    (output/'RUN-R13.txt').write_text('''Round 13 consolidated Excel reading list. No answers have been read.
Run in Windows PowerShell from a copy of this complete directory:
  .\\RUN-R13.ps1 -DnaCalcCheckout C:\\Work\\DnaCalc\\DnaCalc
The script only invokes the existing cargo xtask witness read --func F --inputs FILE
interface. It performs no git operation. Use the Windows lane's approved product
checkout; respect the current product history-rewrite freeze. This packet does
not authorize a product commit/push or any DnaCalc-engine access.

READING-LIST.tsv is the complete list; batches.json supplies the reading order.
Read all original 64 R11 cumulative, 52 R11 structural-v2 and 64 R12 group calls,
then the R13 calls. Original batches and prediction files are byte-identical in
prior/. Separate batches preserve repeated controls: merging by function would
violate the reader's duplicate-call rule. Do not use --new-only. The superseded
original structural-r11 packet is deliberately excluded.

Keep each entire raw worker directory, request.json, result.json, input_audits,
requested/submitted/read-back typed arguments, formulas, raw outcomes, errors,
and the reader's frozen application predictions. Record actual Excel build,
CPU/bitness, Compatibility Version and calculation settings. Retain every failed
or changed-input obligation; never normalize it into a same-input match. A failed
read is a failure and must remain separately visible even if later retried.

Score the immutable prior/ prediction files before subsequent elimination.
The R13 PMT calls repeat all 24 original gaps. POWER/LN/EXP calls observe explicit
worksheet boundaries, using original inputs or fixed native-ISA calculations.
selection-r13.json records every derivation, including duplicate derivations.
predictions/native-elementary.jsonl is model output only: every store mask of
the declared standalone hypotheses, never an Excel answer or surviving PMT code.
Worksheet elementary outputs do not reveal PMT's internal extended operands.
No request contains an expected answer. Do not ingest a prediction as evidence.

Native CPU follow-up: if raw-x87-probe.rs and native-leaf-requests.tsv are present,
run the separate probe commands in NATIVE-PROBE.txt. Its results are CPU instruction
measurements, not Excel readings. Keep those two result sets separate.
''')
    if leaves:
        leaves=Path(leaves)
        shutil.copyfile(leaves/'native-leaf-requests.tsv',output/'native-leaf-requests.tsv')
        shutil.copyfile(leaves/'native-leaf-reference.tsv',output/'predictions'/'native-leaf-reference.tsv')
        shutil.copyfile(ROOT/'racer/x87_probe.rs',output/'raw-x87-probe.rs')
        shutil.copyfile(ROOT/'checker/compare_x87_probe.py',output/'compare-native-probe.py')
        if environment:shutil.copyfile(environment,output/'native-reference-environment.json')
        (output/'RUN-NATIVE-PROBE.ps1').write_text('''param([string]$OutputDirectory = (Join-Path (Split-Path -Parent $PSScriptRoot) 'excel-results-r13-native'))
$ErrorActionPreference = 'Stop'
$packet = $PSScriptRoot
$null = New-Item -ItemType Directory -Force -Path $OutputDirectory
$OutputDirectory = (Resolve-Path -LiteralPath $OutputDirectory).Path
$stamp = Get-Date -Format 'yyyyMMdd-HHmmss'
$exe = Join-Path $OutputDirectory ("raw-x87-probe-{0}.exe" -f $stamp)
& rustc --edition=2024 -O (Join-Path $packet 'raw-x87-probe.rs') -o $exe
if ($LASTEXITCODE -ne 0) { throw 'Native probe compilation failed' }
$result = Join-Path $OutputDirectory ("native-results-{0}.tsv" -f $stamp)
& $exe (Join-Path $packet 'native-leaf-requests.tsv') $result
if ($LASTEXITCODE -ne 0) { throw 'Native probe failed: preserve the partial result as a failure' }
Write-Host ("Native CPU results: {0}" -f $result)
Write-Host 'Return this file with the Excel readings. It contains no Excel answers.'
'''.replace('\n+','\n'))
        (output/'NATIVE-PROBE.txt').write_text('''Run separately from the Excel readings, on the SAME physical Windows CPU:
  .\\RUN-NATIVE-PROBE.ps1
Requires the x86_64 Rust toolchain (the product already uses Rust). This compiles
one standalone file with no product dependency and executes only the requested
FYL2X, FYL2XP1, F2XM1 and ROM-load instructions. Every word is PC53 or PC64, RN.
It restores the caller's control word and refuses to overwrite an output file.
The executable and results go in sibling excel-results-r13-native/ by default;
the frozen packet is unchanged. Override with -OutputDirectory if needed.
No source in either product repository is changed. Never use DnaCalc-engine.

Return native-results-<timestamp>.tsv and its console CPU/OS line. The request
file contains only IDs, instruction names, words and exact input operands. The
reference file is a Linux EPYC instruction measurement, not an Excel answer.
On any machine with Python, compare complete inputs and results with:
  python compare-native-probe.py native-leaf-requests.tsv predictions/native-leaf-reference.tsv native-results-<timestamp>.tsv native-comparison.json
Changed inputs, duplicate/missing rows or reordered IDs reject the comparison.
Any changed leaf keeps the corresponding Linux exclusion conditional on EPYC;
re-run the programs on the Windows CPU before transferring that exclusion.
Equal leaves support only this finite certificate, never general CPU parity.
Do not put either native output or a model prediction into an Excel bank.
''')
    # Schema/coverage validation mirrors the existing reader's admission checks,
    # including duplicate typed calls inside each batch, but performs no read.
    for batch in batches:
        doc=json.loads((output/batch['path']).read_text());ids=set();calls=set()
        if doc['function']!=batch['function']:raise ValueError('batch function differs')
        for row in doc['witnesses']:
            key=(json.dumps(row['args'],sort_keys=True),json.dumps(row.get('result_index')),
                 row.get('input_transport','cell_value2_bulk'))
            if row['id']in ids or key in calls:raise ValueError('duplicate call in batch')
            if set(row)-{'id','args','result_index','input_transport'}:raise ValueError('unexpected request fields')
            ids.add(row['id']);calls.add(key)
    files={str(p.relative_to(output)):digest(p)for p in sorted(output.rglob('*'))if p.is_file()}
    manifest=dict(frozen_utc=datetime.now(timezone.utc).isoformat(),new_excel_answers=0,
        calls=len(reading),origin_calls=dict(origin_counts),batches=len(batches),
        prior_packets=sources,generator_sha256=digest(__file__),trace_sha256=digest(trace),
        original_gap_ids_sha256=digest(gaps),native_library_sha256=digest(library),files_sha256=files,
        purpose='One reading list preserving all180 unanswered obligations plus fixed instruction-boundary audits. No fitted helper or surviving complete PMT program.',
        evidence_policy='Requests only enter Excel. Immutable predictions are scored first and never ingested as answers. CPU microcode outputs remain separate from Excel observations.')
    write(output/'manifest.json',manifest)
    print(json.dumps({k:v for k,v in manifest.items()if k not in ('files_sha256','prior_packets')}))
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('trace','gaps','library','output'):p.add_argument(name)
    p.add_argument('--leaves');p.add_argument('--environment')
    a=p.parse_args();create(a.trace,a.gaps,a.library,a.output,a.leaves,a.environment)
