"""Adaptive binary64 domain questions through the existing Excel witness reader.

Original MIT research tooling. No Excel implementation or second capture path.
The transition search assumes one numeric-to-NUM transition in each supplied
bracket; all queried points and local neighbor checks are retained to test it.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with Path(path).open('x', encoding='utf-8') as out:
        json.dump(value, out, indent=2)
        out.write('\n')


def bits(value):
    return '0x'+struct.pack('>d', float(value)).hex()


def number(value):
    return struct.unpack('>d', bytes.fromhex(value[2:]))[0]


def midpoint(low, high):
    require(0 < number(low) < number(high) and int(high, 16)-int(low, 16) > 1,
            'midpoint needs a nonadjacent positive finite bracket')
    return f'0x{(int(low,16)+int(high,16))//2:016x}'


def create(trace, output):
    """Use R16 readings as discovery; every endpoint is re-read on Windows."""
    rows = [json.loads(line) for line in Path(trace).read_text().splitlines()]
    groups = defaultdict(list)
    for row in rows:
        if row['id'].startswith('r16-') and row['function'] in ('BESSELJ','BESSELY'):
            x,n = map(number, row['args'])
            require(abs(x) >= 2**26 and n == int(n), 'unexpected discovery scope')
            groups[row['function'], int(n), x < 0].append(row)
    plan = []
    for (function,order,negative), group in sorted(groups.items()):
        numeric = [r for r in group if r['expected'].startswith('0x')]
        errors = [r for r in group if r['expected'] == 'error:Num']
        require(numeric, 'missing numeric discovery endpoint')
        low_row = max(numeric,key=lambda r:abs(number(r['args'][0])))
        low = abs(number(low_row['args'][0]))
        candidates = [r for r in errors if abs(number(r['args'][0])) > low]
        high_row = min(candidates,key=lambda r:abs(number(r['args'][0]))) if candidates else None
        high = abs(number(high_row['args'][0])) if high_row else 2**27+4.0
        plan.append(dict(id=f'{function.lower()}-{order}-'+('negative' if negative else 'positive'),
            function=function, order=bits(order), negative=negative, low=bits(low), high=bits(high),
            discovery_low_id=low_row['id'], discovery_high_id=high_row['id'] if high_row else None))
    require(len(plan) == 11, 'R16 J/Y discovery bracket inventory differs')
    # Recheck the settled I boundary, including unmeasured higher orders and
    # negative high orders. These are controls, not assumed endpoints.
    controls = []
    for order in (0,1,2,3,10,100,1570,10000):
        for x in (709.7, math.nextafter(709.7, math.inf)):
            for sign in (1,-1):
                controls.append(dict(function='BESSELI',args=[bits(sign*x),bits(order)]))
    # Separately labelled, unverified phase-guard hypothesis: a private trig
    # limit2^27-1 after subtracting a nine-decimal-place phase. The generic
    # bisection does not depend on it being correct; observed values alone
    # may tighten the starting brackets.
    for state in plan:
        phase = math.ceil((3 if number(state['order']) == 1 else 1)*math.pi/4*10**9)/10**9
        center = int(bits((2**27-1)+phase),16)
        for offset in range(-4,5):
            x = number(f'0x{center+offset:016x}')
            controls.append(dict(function=state['function'],args=[bits(-x if state['negative']else x),state['order']],
                reason='Unverified private phase guard2^27-1; fixed nine-decimal phase, neighboring input lattice.'))
    output = Path(output)
    output.mkdir(parents=True,exist_ok=False)
    write(output/'plan.json',dict(format='dnacalc-binary64-domain-bisection/1',brackets=plan,controls=controls,
        hypothesis='Private Bessel trig may reject a stored phase beyond2^27-1. Phase0/1 are the ceilings of pi/4 and3pi/4 at9decimal places. '
                   'This is an unanswered question, not an adopted cutoff or observed result.',
        transition_assumption='One numeric-to-NUM transition per bracket; observed reversals abort. '
            'Adjacent endpoints prove the measured local transition, not monotonicity at unqueried inputs.',
        maximum_rounds=64,maximum_endpoint_expansions=8,neighbor_radius=2,
        source_trace_sha256=sha(trace),sources='Original R16 Excel Value2 observations; IEEE754 ordered positive binary64 lattice.'))
    shutil.copyfile(__file__,output/'domain_bisection.py')
    (output/'README.txt').write_text('''R19 adaptive Bessel domain packet. Python 3 and Cargo are required on Windows.
Run from any directory, after pulling the published DnaCalc master:
  python domain_bisection.py run --checkout C:\\Work\\DnaCalc\\DnaCalc-win --evidence C:\\Work\\DnaCalc\\OxFunc-Evidence --output C:\\Work\\DnaCalc\\codex-runs\\excel-results-r19-domains

Use a new output directory outside both checkouts and this packet. The script
uses only cargo xtask witness read --inputs; no --new-only and no new oracle.
It rechecks all endpoints, bisects positive binary64 magnitude coordinates
independently per function/order/sign, then asks both endpoints and two exterior
neighbors. I's adjacent 709.7 boundary is separately rechecked through order10000.
Initial controls also surround an unverified2^27-1 private-phase guard; only
their observed answers can shorten a bracket. The full bisection remains valid
if that hypothesis is wrong. No predicted answer is submitted to the reader.
Every requested/submitted/Value2 operand and raw result must agree before a
reading changes a bracket. Predictions never count as answers. Missing results,
unexpected errors, operand changes, conflicting repeated results and observed
nonmonotonicity abort. A batch failure preserves the partial attempt.

The bisection assumes one numeric-to-NUM transition inside each starting bracket;
the packet cannot prove global monotonicity or identify Excel's hidden program.
Initial high endpoints lacking an R16 NUM answer expand by powers of two only
if Excel returns a number. A completed report identifies adjacent local endpoints.
Return the whole output directory: packet, requests, canonical WitnessSets,
worker requests/results, frozen predictions, step journal and file hashes.
The original canonical readings also stay in OxFunc-Evidence/readings.
Do not submit plan.json, summaries or predictions to the witness reader.
''',encoding='utf-8')
    write(output/'manifest.json',dict(format='dnacalc-r19-domain-packet/1',
        files_sha256={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()},
        brackets=len(plan),fixed_controls=len(controls),new_excel_answers=0))


def verify_reading(request, doc, frozen, submitted, worker):
    """Require exact operands and canonical/raw audit agreement before branching."""
    require(doc['function'] == request['function'] and not doc['unread'], 'missing or wrong-function readings')
    provenance = doc['capture_provenance']
    require(provenance['frozen_model_predictions']['inputs_sha256'] == request['_sha256'], 'request hash differs')
    require(provenance['environment']['excel_input_plumbing'] == 'cell_value2_bulk', 'unexpected input surface')
    original = {k:v for k,v in request.items() if k != '_sha256'}
    require(frozen['function'] == submitted['function'] == request['function']
            and frozen['input_document'] == original and frozen['inputs_sha256'] == request['_sha256'],
            'original worker/frozen request differs')
    require(worker['fatal'] is None and not worker['failures'] and not doc['transport_audit_readings'],
            'worker failure or transport obligation')
    require(worker['automation_security_before'] == worker['automation_security_after'], 'worker changed security settings')
    for field in ('automation_security_before','automation_security_after','calculation_settings'):
        require(provenance[field] == worker[field], 'worker provenance differs')
    for canonical,raw in (('excel_version','excel_version'),('excel_build','excel_build'),
                          ('excel_file_version','excel_file_version'),('cpu_id','cpu'),
                          ('excel_operating_system','os_version')):
        require(provenance['environment'][canonical] == worker[raw], 'worker environment differs')
    bitness = {'x86-64':'64-bit','x86':'32-bit','arm64':'64-bit'}.get(worker['excel_machine'])
    require(bitness is None or provenance['environment']['excel_bitness'] == bitness, 'worker bitness differs')
    asked = {r['id']:r for r in request['witnesses']}
    witnesses = {r['id']:r for r in doc['witnesses']}
    audits = {r['id']:r for r in doc['input_audits']}
    calls = {r['id']:r for r in submitted['rows']}
    answers = {r['id']:r for r in worker['rows']}
    predictions = {r['id']:r for r in frozen['predictions']}
    require(len(asked) == len(request['witnesses']) == len(doc['witnesses']) == len(doc['input_audits'])
            and asked.keys() == witnesses.keys() == audits.keys(), 'duplicate, missing or extra reading IDs')
    require(all(len(rows) == len(asked) and {r['id'] for r in rows} == asked.keys()
                for rows in (submitted['rows'],worker['rows'],frozen['predictions'])), 'original worker coverage differs')
    require(all([r['id'] for r in rows] == [r['id'] for r in doc['witnesses']]
                for rows in (doc['input_audits'],submitted['rows'],worker['rows'])), 'original worker order differs')
    results = {}
    for qid,q in asked.items():
        r,a = witnesses[qid],audits[qid]
        require(q['args'] == r['args'] and r.get('result_index') is None and a['reading_failure'] is None,
                'changed arguments or failed reading')
        expected_args = [dict(number=v[2:]) for v in q['args']]
        audit = a['audit']
        call,answer,prediction = calls[qid],answers[qid],predictions[qid]
        require(call['args'] == expected_args and call['argument_literals'] is None
                and call['formula'] == audit['formula'] and answer['input_kept'] is True
                and answer['input_audit'] == audit and answer['value'] == a['observed_raw'],
                'original raw worker/audit differs')
        require(prediction['args'] == q['args'] and prediction.get('result_index') is None
                and prediction['calc_func_prediction'] == r['calc_func_prediction'], 'frozen prediction differs')
        require(audit['requested'] == audit['submitted_value2'] == audit['read_back'] == expected_args,
                'Value2 operand change')
        require(audit['transport'] == 'cell_value2_bulk' and audit['argument_literals'] is None
                and audit['formula'] == '='+request['function']+'(RC[1],RC[2])'
                and audit['effective_callee_inputs'] == 'not directly observable', 'formula/surface changed')
        raw = a['observed_raw']
        if set(raw) == {'number'}:
            expected = '0x'+raw['number']
            require(math.isfinite(number(expected)), 'nonfinite reading')
        elif raw == {'error_code':2036}:
            expected = 'error:Num'
        else:
            raise ValueError('unexpected raw result kind')
        require(r['expected_bits'] == expected, 'canonical/raw result differs')
        results[tuple(q['args'])] = expected
    return results


def search(plan, ask, journal):
    """Pure adaptive state machine; ask returns verified observed outcomes."""
    states = [dict(r) for r in plan['brackets']]
    known = {}
    def query(calls, stage):
        # De-duplicate requests, including endpoint rechecks shared by regions.
        unique = list(dict.fromkeys((f,tuple(args)) for f,args in calls))
        answers = ask(unique,stage)
        require(set(answers) == set(unique), 'incomplete adaptive batch')
        for key,value in answers.items():
            require(value == 'error:Num' or value.startswith('0x') and math.isfinite(number(value)),
                    'unexpected domain result')
            require(key not in known or known[key] == value, 'conflicting repeated reading')
            known[key] = value
        journal(stage,states,answers)
        return answers
    def key(s,magnitude):
        x = number(magnitude)
        return s['function'],(bits(-x if s['negative'] else x),s['order'])
    initial = [key(s,v) for s in states for v in (s['low'],s['high'])]
    initial.extend((r['function'],tuple(r['args'])) for r in plan['controls'])
    query(initial,'endpoints')
    # Use observed interior controls only. This saves reader sessions when an
    # explicitly unverified hypothesis happens to locate the transition.
    for s in states:
        values=[]
        for (function,args),result in known.items():
            if function == s['function'] and args[1] == s['order'] and (number(args[0]) < 0) == s['negative']:
                magnitude=bits(abs(number(args[0])))
                if int(s['low'],16) <= int(magnitude,16) <= int(s['high'],16):values.append((magnitude,result))
        numeric=[v for v,result in values if result.startswith('0x')]
        refused=[v for v,result in values if result=='error:Num']
        if numeric and refused:
            s['low']=max(numeric,key=lambda v:int(v,16))
            s['high']=min(refused,key=lambda v:int(v,16))
            require(int(s['low'],16)<int(s['high'],16),'observed nonmonotone initial controls')
    for expansion in range(plan['maximum_endpoint_expansions']+1):
        more = []
        for s in states:
            require(known[key(s,s['low'])].startswith('0x'), 'numeric discovery endpoint is no longer numeric')
            if known[key(s,s['high'])] != 'error:Num':
                require(expansion < plan['maximum_endpoint_expansions'], 'no NUM endpoint within expansion budget')
                low,high = number(s['low']),number(s['high'])
                s['high'] = bits(low+2*(high-low))
                more.append(key(s,s['high']))
        if not more:
            break
        query(more,f'expand-{expansion:02d}')
    for iteration in range(plan['maximum_rounds']):
        active = [(s,midpoint(s['low'],s['high'])) for s in states if int(s['high'],16)-int(s['low'],16) > 1]
        if not active:
            break
        answers = query([key(s,mid) for s,mid in active],f'bisect-{iteration:02d}')
        for s,mid in active:
            s['high' if answers[key(s,mid)] == 'error:Num' else 'low'] = mid
    require(all(int(s['high'],16)-int(s['low'],16) == 1 for s in states), 'incomplete bisection')
    query([key(s,f'0x{int(s[side],16)+offset:016x}') for s in states
           for side in ('low','high') for offset in range(-plan['neighbor_radius'],plan['neighbor_radius']+1)],'neighbors')
    for s in states:
        for (function,args),expected in known.items():
            if function != s['function'] or args[1] != s['order'] or (number(args[0]) < 0) != s['negative']:
                continue
            x = abs(number(args[0]))
            require((x <= number(s['low'])) == expected.startswith('0x'), 'observed nonmonotone domain')
    return dict(complete=True,brackets=states,observations=len(known),
        limits=plan['transition_assumption'],new_excel_answers=True)


def run(checkout, evidence, output):
    packet = Path(__file__).resolve().parent
    checkout,evidence,output = [Path(p).resolve() for p in (checkout,evidence,output)]
    require(evidence.name == 'OxFunc-Evidence', 'evidence must be OxFunc-Evidence')
    require(not any(output == p or p in output.parents for p in (packet,checkout,evidence)), 'output must be outside protected directories')
    manifest = json.loads((packet/'manifest.json').read_text(encoding='utf-8'))
    require(all(sha(packet/p) == digest for p,digest in manifest['files_sha256'].items()), 'packet hash differs')
    output.mkdir(parents=True,exist_ok=False)
    shutil.copytree(packet,output/'packet')
    plan = json.loads((packet/'plan.json').read_text(encoding='utf-8'))
    roots = {'OxFunc-Evidence/readings':evidence/'readings', 'worker/read':checkout/'target/witness/read',
             'worker/frozen':checkout/'target/witness/frozen'}
    snapshots = {name:{p for p in root.rglob('*') if p.is_file()} for name,root in roots.items()}
    batches = []
    def ask(calls,stage):
        results = {}
        for function in sorted({f for f,_ in calls}):
            rows = [dict(id=f'r19-{stage}-{function.lower()}-{i:03d}',args=list(args),input_transport='cell_value2_bulk')
                    for i,(f,args) in enumerate(calls) if f == function]
            path = output/f'requests-{stage}-{function}.json'
            request = dict(function=function,selection='R19 adaptive binary64 domain bisection; endpoint rechecks intentional',witnesses=rows)
            write(path,request)
            request['_sha256'] = sha(path)
            env = dict(os.environ,DNACALC_LEGACY_ROOT=str(evidence.parent))
            before = {p for p in (evidence/'readings').rglob('*.json')}
            workers_before = {p for p in (checkout/'target/witness/read').glob('*') if p.is_dir()}
            with (output/f'read-{stage}-{function}.log').open('x',encoding='utf-8') as log:
                # Operator-selected local checkout and Cargo, argv only, no shell.
                subprocess.run(['cargo','xtask','witness','read','--func',function,'--inputs',str(path),'--timeout','900'],
                               cwd=checkout,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
            found = []
            for candidate in sorted(set((evidence/'readings').rglob('*.json'))-before):
                doc = json.loads(candidate.read_text(encoding='utf-8-sig'))
                if doc.get('capture_provenance',{}).get('frozen_model_predictions',{}).get('inputs_sha256') == request['_sha256']:
                    found.append((candidate,doc))
            require(len(found) == 1,'missing or duplicate canonical batch')
            reading,doc = found[0]
            model = doc['capture_provenance']['frozen_model_predictions']
            model_name = model['frozen_path'].replace('\\','/').split('/')[-1]
            frozen_path = checkout/'target/witness/frozen'/model_name
            require(sha(frozen_path) == model['frozen_sha256'], 'original frozen prediction hash differs')
            # Retries remain separate raw attempts. Join only a unique complete
            # worker result with the original request IDs, never by a prediction.
            worker_paths = []
            workers_after = {p for p in (checkout/'target/witness/read').glob('*') if p.is_dir()}
            for directory in sorted(workers_after-workers_before):
                result_path = directory/'result.json'
                submitted_path = directory/'request.json'
                if not result_path.is_file() or not submitted_path.is_file():
                    continue
                raw = json.loads(result_path.read_text(encoding='utf-8-sig'))
                submitted = json.loads(submitted_path.read_text(encoding='utf-8-sig'))
                if (submitted.get('function') == function and raw.get('fatal') is None and not raw.get('failures')
                        and len(raw.get('rows',[])) == len(rows)
                        and {r['id'] for r in raw['rows']} == {r['id'] for r in rows}):
                    worker_paths.append(directory)
            require(len(worker_paths) == 1,'missing or duplicate complete raw worker attempt')
            directory = worker_paths[0]
            bound_paths = (frozen_path,directory/'request.json',directory/'result.json',reading)
            bound_hashes = {str(p):sha(p) for p in bound_paths}
            frozen,submitted,worker = [json.loads(p.read_text(encoding='utf-8-sig')) for p in bound_paths[:3]]
            answers = verify_reading(request,doc,frozen,submitted,worker)
            require(bound_hashes == {str(p):sha(p) for p in bound_paths},'worker files changed during decision')
            results.update({(function,args):value for args,value in answers.items()})
            batches.append(dict(function=function,path=path.name,calls=len(rows),reading=str(reading.relative_to(evidence)),
                                reading_sha256=sha(reading),request_sha256=request['_sha256'],
                                original_worker_files_sha256=bound_hashes))
        return results
    def journal(stage,states,answers):
        write(output/f'journal-{stage}.json',dict(states=states,answers=[dict(function=f,args=args,expected=v)for(f,args),v in answers.items()]))
    failure = None
    try:
        result = search(plan,ask,journal)
        write(output/'summary.json',result)
    except BaseException as error:
        failure = str(error)
        raise
    finally:
        for name,root in roots.items():
            for path in sorted({p for p in root.rglob('*') if p.is_file()}-snapshots[name]):
                dest = output/name/path.relative_to(root)
                dest.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(path,dest)
        write(output/'batches.json',batches)
        write(output/'run-status.json',dict(complete=failure is None,failure=failure,canonical_batches=len(batches)))
        write(output/'results-sha256.json',{p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file()})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command',required=True)
    create_parser = sub.add_parser('create')
    create_parser.add_argument('trace')
    create_parser.add_argument('output')
    runner = sub.add_parser('run')
    for option in ('checkout','evidence','output'):
        runner.add_argument('--'+option,required=True)
    args = parser.parse_args()
    if args.command == 'create':
        create(args.trace,args.output)
    else:
        run(args.checkout,args.evidence,args.output)
