"""Freeze simple function calls from existing corpus comparison records.

This does not capture Excel or calculate expected values. It keeps the original
capture, evaluation, formula, environment and hashes. Proposed typed arguments
must subsequently pass the application's workbook parser/registered-route audit.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import struct
import sys
from decimal import Decimal, localcontext
from fractions import Fraction


def error_size(expected, actual):
    """Report finite binary64 error without overflowing subtraction/division.

    Decimal strings are diagnostics, rounded to 18 significant digits from
    exact rational differences. They never participate in the verdict.
    Relative error uses abs(expected); it is undefined at an expected zero.
    """
    import math
    x, y = (struct.unpack('>d', bytes.fromhex(b))[0] for b in (expected, actual))
    if not (math.isfinite(x) and math.isfinite(y)):
        return dict(absolute_error=None, relative_error=None)
    difference = abs(Fraction(y)-Fraction(x))
    def decimal(q):
        with localcontext() as context:
            context.prec = 18
            return str(Decimal(q.numerator)/Decimal(q.denominator))
    return dict(absolute_error=decimal(difference),
                relative_error=decimal(difference/abs(Fraction(x))) if x else None)


def hypergeom_symmetry(executable, input_path, output):
    """Check fixed exact-distribution symmetries through the current invoke route.

    Conditional arithmetic diagnostics, never a replacement numerical kernel.
    Integer finite-support rows only; CDF permits marginal interchange alone.
    """
    import math
    import subprocess
    executable = Path(executable).resolve()
    input_path, output = Path(input_path), Path(output)
    output.mkdir()
    originals = [json.loads(line) for line in input_path.read_text().splitlines()]
    bits = lambda v: '0x'+struct.pack('>d', float(v)).hex()
    def transform(s, n, k, total, identity):
        if identity & 2:
            s, k = n-s, total-k
        if identity & 4:
            s, n = k-s, total-n
        if identity & 1:
            n, k = k, n
        return s, n, k, total
    # Independent integer combinatorial controls for the proposed mappings.
    controls = 0
    for total in range(1, 9):
        for n in range(total+1):
            for k in range(total+1):
                for s in range(max(0, n+k-total), min(n, k)+1):
                    exact = Fraction(math.comb(k, s)*math.comb(total-k, n-s), math.comb(total, n))
                    for identity in range(8):
                        ss, nn, kk, _ = transform(s, n, k, total, identity)
                        other = Fraction(math.comb(kk, ss)*math.comb(total-kk, nn-ss), math.comb(total, nn))
                        if exact != other:
                            raise ValueError('incorrect hypergeometric identity')
                        controls += 1
    questions, excluded = [], Counter()
    for index, row in enumerate(originals):
        if row['function'] != 'HYPGEOM.DIST':
            excluded['other_function'] += 1
            continue
        try:
            a = [struct.unpack('>d', bytes.fromhex(v[2:]))[0] for v in row['args']]
            s, n, k, total, flag = a
            if not all(math.isfinite(v) and v == math.trunc(v) for v in a[:4]):
                raise ValueError('not finite integer')
        except (ValueError, TypeError, struct.error):
            excluded['outside_scalar_integer_shape'] += 1
            continue
        if not (0 <= s <= min(n, k) and n+k-s <= total and 0 <= n <= total <= 10000 and 0 <= k <= total):
            excluded['outside_integer_finite_support'] += 1
            continue
        for identity in range(2 if flag else 8):
            args = (*transform(s, n, k, total, identity), flag)
            questions.append(dict(id=f'{index}:{identity}', args=list(map(bits, args))))
    request = output/'inputs.json'
    request.write_text(json.dumps(dict(function='HYPGEOM.DIST', witnesses=questions))+'\n')
    predictions = output/'predictions.jsonl'
    # ubs:ignore — trusted local executable, argument vector only; failed invocation aborts the receipt.
    subprocess.run([str(executable), 'invoke-inputs', str(request), str(predictions)], check=True, timeout=300)
    results = [json.loads(line) for line in predictions.read_text().splitlines()]
    if len(results) != len(questions):
        raise ValueError('missing symmetry predictions')
    counts, first_losses = {}, {}
    for question, result in zip(questions, results):
        if question['id'] != result['id'] or question['args'] != result['args']:
            raise ValueError('changed symmetry request')
        index, identity = map(int, question['id'].split(':'))
        original = originals[index]
        was = original['actual'] == original['expected']
        now = result['calc_func_prediction'] == original['expected']
        if identity == 0 and result['calc_func_prediction'] != original['actual']:
            raise ValueError('baseline route changed')
        tally = counts.setdefault(str(identity), [0]*5)
        for j, value in enumerate((1, now, was, now and not was, was and not now)):
            tally[j] += value
        if was and not now:
            first_losses.setdefault(str(identity), dict(original=original, transformed=result))
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    summary = dict(complete=True, original_rows=len(originals), predictions=len(results),
        combinatorial_controls=controls, counts_all_exact_baseline_gains_losses=counts,
        excluded=excluded, first_losses=first_losses,
        input_sha256={str(p):sha(p) for p in (input_path, executable, Path(__file__), request, predictions)},
        source='Hypergeometric combination definition and exact marginal/sample/success-complement identities; no third-party code.',
        limits='Current invoke is a conditional numerical leaf. No fitted parameter, kernel identification, full-bank improvement or product change. CDF is tested only under marginal interchange.')
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('input_sha256', 'first_losses')}))


def domain_orders(inputs, output):
    """Eliminate fixed validation orders using recorded result types only.

    This certificate identifies branch order on the stated domain; numeric
    values still require the application/bank comparison. No answer enters a
    proposed computation and no numerical parameter is inferred.
    """
    import itertools
    import math
    paths = [Path(p) for p in inputs.split(',')]
    groups = {'POISSON.DIST': [], 'HYPGEOM.DIST': []}
    omitted = Counter()
    for path in paths:
        for index, row in enumerate(map(json.loads, path.read_text().splitlines())):
            name = row['function']
            if name not in groups:
                omitted['other_function'] += 1
                continue
            try:
                a = [struct.unpack('>d', bytes.fromhex(v[2:]))[0] for v in row['args'][:-1]]
                if not all(math.isfinite(v) for v in a):
                    raise ValueError('nonfinite')
                if name == 'POISSON.DIST':
                    # EXP overflow is deliberately outside this type-only proof.
                    if not (a[0] < 0 or -709 <= a[1] < 0):
                        omitted['poisson_ordinary_or_overflow'] += 1
                        continue
                elif not (any(v < 0 for v in a) or a[3] == 0):
                    omitted['hypergeometric_ordinary'] += 1
                    continue
                expected = row['expected']
                if expected == 'error:Num':
                    kind = 'num_error'
                elif expected.startswith('0x'):
                    kind = 'number'
                else:
                    raise ValueError('outside_number_or_num')
            except (ValueError, TypeError, struct.error):
                omitted['input_or_output_shape'] += 1
                continue
            groups[name].append(dict(source=str(path), index=index, id=row['id'], args=a,
                                     expected=kind))
    results = {}
    for name, rows in groups.items():
        if not rows:
            raise ValueError('empty domain-order scope '+name)
        profiles = (list(itertools.product((False, True), repeat=3)) if name == 'POISSON.DIST'
                    else list(itertools.product((False, True), repeat=2)))
        records, survivors = [], []
        for identity, profile in enumerate(profiles):
            for order, row in enumerate(rows):
                a = row['args']
                if name == 'POISSON.DIST':
                    raw_sign, zero_before_mean, truncate_zero = profile
                    x, mean = a
                    count = math.trunc(x)
                    zero = (count == 0) if truncate_zero else (x == 0)
                    if (x if raw_sign else count) < 0:
                        got = 'num_error'
                    elif zero_before_mean and zero:
                        got = 'number'
                    else:
                        got = 'num_error' if mean < 0 else 'number'
                else:
                    raw_sign, allow_empty = profile
                    s, n, k, total = map(math.trunc, a)
                    negative = any(v < 0 for v in (a if raw_sign else (s,n,k,total)))
                    got = ('num_error' if negative or total == 0 and not allow_empty
                           or n > total or k > total else 'number')
                if got != row['expected']:
                    records.append(dict(candidate_id=identity, profile=profile, order=order,
                                        witness=row, got=got))
                    break
            else:
                survivors.append(dict(candidate_id=identity, profile=profile))
        results[name] = dict(profiles=len(profiles), observations=len(rows),
                            first_kills=records, survivors=survivors)
    result = dict(complete=True, functions=results, omitted=omitted,
                  inputs_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
                  limits='Result types only on finite numeric scalar inputs. Poisson EXP overflow and '
                         'all ordinary numerical bodies remain outside this certificate. Discovery '
                         'observations, no fitted parameter; full application bank replay supplies value evidence.')
    with Path(output).open('x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps({k:dict(profiles=v['profiles'],observations=v['observations'],survivors=v['survivors'])
                      for k,v in results.items()}))
    return result


def discrete_questions(output, round_number='20'):
    """Small fixed neighbor/symmetry questions, without predicted answers."""
    import math
    if round_number not in ('20', '21'):
        raise ValueError('unknown discrete question round')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    bits = lambda v: '0x'+struct.pack('>d',float(v)).hex()
    groups = {'BINOM.DIST': [], 'HYPGEOM.DIST': []}
    reasons = []
    def add(name, arguments, purpose):
        identity = 'r'+round_number+'-'+name.lower().replace('.','-')+'-'+str(len(groups[name])).zfill(3)
        groups[name].append(dict(id=identity,args=list(map(bits,arguments)),input_transport='cell_value2_bulk'))
        reasons.append(dict(id=identity,function=name,purpose=purpose))
    # These paired tails distinguish independent saddle-point evaluation from
    # symmetric or recurrent consumers. Neighboring probabilities distinguish
    # complement rounding without tuning a coefficient or branch boundary.
    for n, probability, counts in (((1000,.5,(1,999)),(100,.3,(1,99)),(1000,.02,(333,667))) if round_number == '20' else ()):
        for p in (math.nextafter(probability,0.0),probability,math.nextafter(probability,1.0)):
            for which, k in enumerate(counts):
                for cumulative in (0,1):
                    add('BINOM.DIST',(k,n,p if which==0 else 1.0-p,cumulative),
                        'Mirrored counts at rounded complementary probabilities and adjacent probability bits; PMF/CDF distinguish independent masses from recurrence or beta consumers. A rounded complement need not be the exact real complement.')
    for a in (((0,8.5,3,8.1),(0,8.1,3,8.5),(0,3,8.5,8.1),(0,3,8.1,8.5)) if round_number == '20' else ()):
        for cumulative in (0,1):
            add('HYPGEOM.DIST',(*a,cumulative),
                'Raw upper-domain comparisons versus truncation before comparison; both cumulative flags. Negative signs/empty populations are already recorded.')
    if round_number == '21':
        for s, n, k, total in ((60,200,300,1000), (2,5,7,20)):
            for identity in range(8):
                ss, nn, kk = s, n, k
                if identity & 2:
                    ss, kk = nn-ss, total-kk
                if identity & 4:
                    ss, nn = kk-ss, total-nn
                if identity & 1:
                    nn, kk = kk, nn
                add('HYPGEOM.DIST', (ss, nn, kk, total, 0),
                    'Integer finite-support PMF under marginal interchange and exact sample/success complements; fixed identity '+str(identity))
            for nn, kk in ((n,k), (k,n)):
                add('HYPGEOM.DIST', (s,nn,kk,total,1),
                    'CDF marginal interchange; paired with the corresponding PMF. Rounded Excel symmetry is a question, never assumed.')
    batches = []
    for name, rows in groups.items():
        if not rows:
            continue
        filename = 'requests-'+name+'.json'
        (output/filename).write_text(json.dumps(dict(function=name,
            selection='R'+round_number+' fixed discriminating questions, with intentional recorded controls; no --new-only',
            witnesses=rows),indent=2)+'\n')
        batches.append(dict(function=name,path=filename,calls=len(rows)))
    (output/'batches.json').write_text(json.dumps(batches,indent=2)+'\n')
    (output/'selection.json').write_text(json.dumps(dict(calls=len(reasons),questions=reasons,
        new_excel_answers=0,rule='Fixed mathematical symmetry, adjacent binary64 probabilities and validation-order alternatives. No coefficient fitting or expected-answer file.'),indent=2)+'\n')
    (output/'README.txt').write_text(('''R20 discrete distribution questions: 36 BINOM.DIST calls and8 HYPGEOM.DIST calls.
From the published DnaCalc Windows checkout, use the existing witness reader:
  cargo xtask witness read --func BINOM.DIST --inputs <packet>/requests-BINOM.DIST.json
  cargo xtask witness read --func HYPGEOM.DIST --inputs <packet>/requests-HYPGEOM.DIST.json
Do not use --new-only. Center points intentionally repeat recorded controls.
Only requests-*.json are reader inputs. Predictions and route audits are current
application calculations, never Excel answers or expected values.

Preserve the full requested/submitted/Value2 bit audits, worker settings,
canonical WitnessSets, original worker requests/results and frozen predictions.
Return those files, the unchanged packet, file hashes and partial attempts in
~/codex-runs/excel-results-r20. Preserve canonical readings in OxFunc-Evidence.
Input changes, missing/unread results and automation failures remain failures.
This packet uses the established witness reader and capture lineage only.
''' if round_number == '20' else '''R21 hypergeometric questions: 20 scalar calls; no expected answers.
From the published DnaCalc Windows checkout, use the existing witness reader:
  cargo xtask witness read --func HYPGEOM.DIST --inputs <packet>/requests-HYPGEOM.DIST.json
Do not use --new-only: recorded center controls are intentional. Only requests-*.json
are reader inputs. Predictions are application calculations, never Excel readings.
Preserve the canonical WitnessSet in OxFunc-Evidence, requested/submitted/Value2
bit audits, worker settings, original requests/results, frozen files, packet hashes
and partial attempts. Return them under ~/codex-runs/excel-results-r21.
No second oracle or capture lineage is introduced. Missing/unread outcomes fail.
'''))
    print(json.dumps(dict(calls=len(reasons),batches=batches,new_excel_answers=0)))


def round21_questions(executable, output):
    """Extend the existing scalar packet and catalogue section question formats."""
    import subprocess
    executable, output = Path(executable).resolve(), Path(output)
    discrete_questions(str(output), '21')
    questions = []
    def add(formula, reason, **metadata):
        row = len(questions)+1
        questions.append(dict(id=f'r21-sheet-{row:04d}', row=row, column=4,
                              formula=formula, reason=reason, **metadata))
    literal = lambda x: format(float(x), '.17g')
    array = lambda rows: '{'+';'.join(','.join(literal(x) for x in row) for row in rows)+'}'
    for y in (1,2,4,7):
        add(literal(y), 'response operand control')
    shapes = [('exact-control', 0.0, 1.0, 1.0, 0.0, False, False),
              ('resolved-control', 2**-16, 1.0, 1.0, 0.0, False, False)]
    for label, scale1, scale2, shift, reverse, swap in (
        ('base', 1,1,0,False,False), ('reverse',1,1,0,True,False),
        ('small-both',2**-20,2**-20,0,False,False),
        ('large-both',2**20,2**20,0,False,False),
        ('small-second',1,2**-20,0,False,False),
        ('large-second',1,2**20,0,False,False),
        ('translated',1,1,2**20,False,False),
        ('swap-columns',1,1,0,False,True)):
        shapes.append((label,2**-18,scale1,scale2,shift,reverse,swap))
    for label, delta, scale1, scale2, shift, reverse, swap in shapes:
        x = [[v*scale1+shift, (2*v+(delta if i==3 else 0))*scale2+shift]
             for i,v in enumerate((1,2,5,4))]
        y = [[v] for v in (1,2,4,7)]
        if reverse:
            x, y = x[::-1], y[::-1]
        if swap:
            x = [row[::-1] for row in x]
        xs, ys = array(x), array(y)
        for r in range(1,5):
            for c in range(1,3):
                add(f'INDEX({xs},{r},{c})', 'matrix operand echo; Excel parser bits are observed, never assumed',
                    shape=label, matrix_index=[r,c])
        for function in ('LINEST','LOGEST'):
            for constant in ('FALSE','TRUE'):
                for r,c in ((1,1),(1,2),(4,2)):
                    add(f'INDEX({function}({ys},{xs},{constant},TRUE),{r},{c})',
                        'fixed rank/removed-column question under scaling, translation or permutation; no threshold fitted',
                        shape=label, function=function, result_index=[r,c], use_const=constant=='TRUE')
    sheet = dict(sheet='dcj9k_r21', questions=questions)
    sheet_path = output/'catalogue-questions.json'
    sheet_path.write_text(json.dumps(sheet,indent=2)+'\n')
    section = output/'sheet-dcj9k_r21.txt'
    section.write_text('# R21: numeric rank transformations and operand echoes. Questions only.\n'+
                       '\n'.join(q['formula'] for q in questions)+'\n',encoding='ascii')
    # ubs:ignore — trusted local executable and argv; missing application values abort publication.
    subprocess.run([str(executable),'sheet-predictions',str(sheet_path),str(output/'application-predictions-catalogue.jsonl')],check=True,timeout=300)
    predictions = [json.loads(line)for line in(output/'application-predictions-catalogue.jsonl').read_text().splitlines()]
    if len(predictions)!=len(questions)or any(q['id']!=p['id']or q['formula']!=p['formula']for q,p in zip(questions,predictions)):
        raise ValueError('missing or changed catalogue question')
    request = output/'requests-HYPGEOM.DIST.json'
    prediction = output/'application-predictions-HYPGEOM.DIST.jsonl'
    # ubs:ignore — same operator-selected executable, argument vector only.
    subprocess.run([str(executable),'invoke-inputs',str(request),str(prediction),'1900'],check=True,timeout=300)
    batch = json.loads(request.read_text())['witnesses']
    calculated = [json.loads(line)for line in prediction.read_text().splitlines()]
    if len(batch)!=len(calculated):
        raise ValueError('missing scalar question prediction')
    audit = output/'application-route-inputs.jsonl'
    with audit.open('x')as f:
        for q,p in zip(batch,calculated):
            if q['id']!=p['id']or q['args']!=p['args']or p['date_system']!='1900':
                raise ValueError('changed scalar question input/environment')
            f.write(json.dumps(dict(**q,function='HYPGEOM.DIST',result_index=None,
                actual=p['calc_func_prediction'],numeric_cell_inputs=True,date_system_1900=True,
                not_excel_evidence=True))+'\n')
    # ubs:ignore — same executable and argv; audit assertions must succeed.
    subprocess.run([str(executable),'function-routes',str(audit),str(output/'application-route-audit.jsonl')],check=True,timeout=300)
    with(output/'README.txt').open('a')as f:
        f.write('''
The same packet asks 204 scalar formula cells on catalogue sheet dcj9k_r21:
120 LINEST/LOGEST selected outputs and84 operand echoes. They extend the one
catalogue oracle; the existing generator registers sheet-dcj9k_r21.txt.
From the published DnaCalc checkout on Windows:
  cargo xtask capture --pending --list
  cargo xtask capture --pending
Keep sections/dcj9k_r21.capture.json and the original worker capture outputs.
Do not place application predictions into a capture or a WitnessSet. The operand
echoes record Excel's literal parser results; requested binary64 arithmetic and
parsed literals need not have identical bits. These are discovery questions,
not fitted rank rules or a coefficient program. Return all failures and attempts.
''')
    sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    summary = dict(format='dnacalc-research-questions-r21/1',scalar_calls=20,
        catalogue_formula_cells=len(questions),selected_regression_outputs=120,operand_echoes=84,
        new_excel_answers=0,catalogue_sheet='dcj9k_r21',catalogue_section_sha256=sha(section),
        prediction_executable_sha256=sha(executable),source_sha256=sha(Path(__file__)),
        limits='Predictions from the frozen R21 pre-merge research executable; all are application calculations, never Excel evidence. Existing standard scalar reader and catalogue capture only; no alternative oracle.')
    summary['files_sha256']={p.name:sha(p)for p in sorted(output.iterdir())if p.is_file()}
    (output/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if 'sha256'not in k}))


def round22_questions(executable, output):
    """Public I0/neighbor controls for the complete fixed ITK program.

    Scalar Value2 transport only; this extends the established packet format
    and creates no catalogue sheet or new capture lineage.
    """
    import math
    import subprocess
    executable, output = Path(executable).resolve(), Path(output)
    # ubs:ignore — trusted local executable, argv only. Its compiled unit binds
    # the application sources; the working checkout may have advanced since.
    binding=json.loads(subprocess.check_output([str(executable),'provenance'],timeout=60))
    output.mkdir(parents=True, exist_ok=False)
    bits = lambda v: '0x'+struct.pack('>d', float(v)).hex()
    questions, reasons = [], []
    def add(x, order, reason):
        identity = f'r22-besseli-{len(questions):03d}'
        questions.append(dict(id=identity, args=[bits(x),bits(order)], input_transport='cell_value2_bulk'))
        reasons.append(dict(id=identity, function='BESSELI', purpose=reason))
    for x in (1.8661547234119729,3.2741177381715616,5.432532676165661,
              6.577498665623793,12.514410892637624,13.469754413420809):
        for order in (0,2,7):
            add(x,order,'Public I0 and fixed neighboring integer orders at prior arithmetic discriminators; expose normalization independently of a fitted recurrence.')
        add(-x,7,'Odd-order parity control at the same magnitude.')
    for x in (math.nextafter(3.75,0),3.75,math.nextafter(3.75,math.inf)):
        for order in (0,2,17):
            add(x,order,'Adjacent binary64 inputs at the published I0 branch boundary; no boundary is estimated from answers.')
    for order in range(1568,1573):
        for x in (709.7,-709.7):
            add(x,order,'Neighbors of an independently recorded normalization-underflow case, including even/odd parity. Existing center readings are intentional controls.')
    for x,order in ((0,0),(0,2),(.1,2),(3,2),(10,20),(100,100)):
        add(x,order,'Public zero, small-argument and rescaled-recurrence controls for the same fixed program.')
    request=output/'requests-BESSELI.json'
    request.write_text(json.dumps(dict(function='BESSELI',selection='R22 fixed complete-program controls; intentional repeats, no --new-only',witnesses=questions),indent=2)+'\n')
    predictions=output/'application-predictions-BESSELI.jsonl'
    # ubs:ignore — trusted local executable, argument vector; failures abort.
    subprocess.run([str(executable),'invoke-inputs',str(request),str(predictions),'1900'],check=True,timeout=300)
    rows=[json.loads(line) for line in predictions.read_text().splitlines()]
    if len(rows)!=len(questions):
        raise ValueError('missing R22 application prediction')
    route=output/'application-route-inputs.jsonl'
    with route.open('x') as out:
        for question,prediction in zip(questions,rows):
            if any(prediction[key]!=question[key] for key in ('id','args')):
                raise ValueError('changed R22 question input')
            out.write(json.dumps(dict(**question,function='BESSELI',result_index=None,
                actual=prediction['calc_func_prediction'],numeric_cell_inputs=True,
                date_system_1900=True,not_excel_evidence=True))+'\n')
    # ubs:ignore — same trusted executable and established route audit.
    subprocess.run([str(executable),'function-routes',str(route),str(output/'application-route-audit.jsonl')],check=True,timeout=300)
    (output/'batches.json').write_text(json.dumps([dict(function='BESSELI',path=request.name,calls=len(questions))],indent=2)+'\n')
    (output/'selection.json').write_text(json.dumps(dict(calls=len(questions),questions=reasons,new_excel_answers=0,
        rule='Published branch, integer-order/parity neighbors and public I0 controls. No fitted coefficient, threshold or expected answer.'),indent=2)+'\n')
    (output/'README.txt').write_text('''R22 BESSELI scalar questions; no Excel answers included.
From the published Windows DnaCalc checkout:
  cargo xtask witness read --func BESSELI --inputs <packet>/requests-BESSELI.json
Do not use --new-only: independently recorded controls intentionally repeat.
Only requests-BESSELI.json is reader input. Application predictions and route
audits are calculations, never Excel answers or expected values.
Preserve canonical WitnessSets in OxFunc-Evidence, requested/submitted/Value2
bit audits, worker settings, original requests/results, file hashes and partial
attempts. Return those files and the unchanged packet under
~/codex-runs/excel-results-r22. Missing/unread outcomes remain failures.

Do not append research questions to the captured catalogue workbook. R21's
204 cells now live in crates/calc-func/tests/fixtures/dcj9k_r21/oracle.xlsx;
fn2 questions live in fixtures/fn2_questions/oracle.xlsx. The original catalogue
has been restored byte-for-byte and needs no recapture for those1633 old cells.
Use cargo xtask capture --pending for the separate pending workbooks.
This supersedes only R21's instruction to capture its questions in catalogue;
the historical R21 request/answer formats and hashes remain unchanged.
''')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest=dict(format='dnacalc-research-questions-r22/1',scalar_calls=len(questions),
        catalogue_formula_cells=0,new_excel_answers=0,application_revision=binding['application_revision'],
        application_working_blobs=binding.get('application_working_blobs',{}),prediction_unit=binding['id'],
        prediction_executable_sha256=sha(executable),source_sha256=sha(Path(__file__)),
        limits='Fixed discovery questions; registered and workbook-parser routes audited. No new Excel evidence, catalogue mutation or replacement oracle.')
    manifest['files_sha256']={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if 'sha256' not in k}))


def round22_qr_questions(executable, output):
    """Add an isolated regression workbook packet; preserve the scalar packet."""
    import subprocess
    executable=Path(executable).resolve();output=Path(output)/'qr'
    output.mkdir(parents=True,exist_ok=False)
    # ubs:ignore — selected local executable and argv; immutable compiled binding.
    binding=json.loads(subprocess.check_output([str(executable),'provenance'],timeout=60))
    questions=[]
    def add(formula,kind,shape):
        questions.append(dict(id=f'r22-qr-sheet-{len(questions)+1:04d}',row=len(questions)+1,
            column=4,formula=formula,kind=kind,shape=shape))
    literal=lambda v:format(float(v),'.17g')
    array=lambda rows:'{'+';'.join(','.join(literal(x)for x in row)for row in rows)+'}'
    def fit(label,y,x,constant,outputs):
        ys,xs=array([[v]for v in y]),array(x)
        for i in range(len(y)):
            add(f'INDEX({ys},{i+1},1)','operand_echo',label)
            for j in range(len(x[0])):add(f'INDEX({xs},{i+1},{j+1})','operand_echo',label)
        for row,column in outputs:
            add(f'INDEX(LINEST({ys},{xs},{str(constant).upper()},TRUE),{row},{column})','selected_output',label)
    y=[31.78,5.61,-9.09]
    x=[[53.56,30.22,32],[-9,-17.977,47.17],[9.4,-11.83,39.39]]
    coefficients=[(1,c)for c in range(1,4)]
    fit('recorded_control',y,x,False,coefficients)
    fit('reverse_rows',y[::-1],x[::-1],False,coefficients)
    fit('swap_columns',y,[[r[1],r[0],r[2]]for r in x],False,coefficients)
    fit('dyadic_column_scale',y,[[r[0]*2**20,r[1],r[2]]for r in x],False,coefficients)
    y=[-17.979,-.126,8,-5,1.7,37,12.4,14.46,-8,24.2]
    x=[[20.8,50],[34,23.41],[38.875,-17.7],[4.863,36.9],[42.1,6],[-14.55,44.41],[-4.57,-12.96],[-9.909,2.589],[47.549,42],[54.9,-17.63]]
    fit('reverse_statistical_fit',y[::-1],x[::-1],True,[(r,c)for r in range(1,6)for c in range(1,4)])
    sheet=output/'workbook-questions.json';sheet.write_text(json.dumps(dict(sheet='dcj9k_r22',questions=questions),indent=2)+'\n')
    section=output/'sheet-dcj9k_r22.txt';section.write_text('# R22: full-rank reflector transforms and operand echoes. Questions only.\n'+'\n'.join(q['formula']for q in questions)+'\n',encoding='ascii')
    predictions=output/'application-predictions-workbook.jsonl'
    # ubs:ignore — same trusted executable, registered workbook parser route.
    subprocess.run([str(executable),'sheet-predictions',str(sheet),str(predictions)],check=True,timeout=300)
    rows=list(map(json.loads,predictions.read_text().splitlines()))
    if len(rows)!=len(questions)or any(p['id']!=q['id']or p['formula']!=q['formula']for p,q in zip(rows,questions)):
        raise ValueError('missing or changed R22 regression question')
    (output/'README.txt').write_text('''R22 regression workbook questions:105 formulas,27 selected LINEST outputs and78 operand echoes.
These are discovery questions, not Excel answers. The original scalar R22 packet stays frozen.
The section belongs only in crates/calc-func/tests/fixtures/dcj9k_r22/oracle.xlsx.
It must never be appended to the captured catalogue. The existing
make-catalogue-fixtures.ps1 generates this separate workbook from the identical
committed sheet-dcj9k_r22.txt. From the published Windows checkout:
  cargo xtask capture --pending --list
  cargo xtask capture --pending
Preserve all original worker outputs, saved packages, settings, attempts and
canonical captures; return them with the unchanged packet under
~/codex-runs/excel-results-r22. Predictions are application calculations only.
The original3x3 center is an intentional control. Reversed rows, swapped columns
and a fixed dyadic column scale ask about arithmetic invariance; no Excel
invariance is assumed. Operand echoes measure literal-parser bits separately.
No coefficient, threshold, rank tolerance or per-row arithmetic choice is fitted.
''')
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    manifest=dict(format='dnacalc-research-questions-r22-qr/1',formula_cells=len(questions),
        selected_outputs=sum(q['kind']=='selected_output'for q in questions),operand_echoes=sum(q['kind']=='operand_echo'for q in questions),
        new_excel_answers=0,application_revision=binding['application_revision'],application_working_blobs=binding.get('application_working_blobs',{}),
        prediction_unit=binding['id'],prediction_executable_sha256=sha(executable),source_sha256=sha(Path(__file__)),
        section_sha256=sha(section),files_sha256={p.name:sha(p)for p in sorted(output.iterdir())if p.is_file()},
        limits='Separate new workbook via the existing generator/capture lineage. All prior captures remain unchanged; every new cell awaits capture and fails the inclusive denominator.')
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items()if 'sha256'not in k}))


def split_arguments(text, separator=','):
    depth = 0
    quoted = False
    start = 0
    result = []
    i = 0
    while i < len(text):
        c = text[i]
        if c == '"':
            if quoted and i+1 < len(text) and text[i+1] == '"':
                i += 2
                continue
            quoted = not quoted
        elif not quoted:
            if c in '({': depth += 1
            if c in ')}': depth -= 1
            if c == separator and depth == 0:
                result.append(text[start:i].strip())
                start = i+1
        i += 1
    if quoted or depth:
        raise ValueError('unbalanced argument')
    result.append(text[start:].strip())
    return result


def scalar(text):
    if text in ('TRUE', 'FALSE'):
        return text == 'TRUE'
    if text.startswith('"') and text.endswith('"'):
        return text[1:-1].replace('""', '"')
    if not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', text):
        raise ValueError('argument requires formula evaluation')
    return '0x'+struct.pack('>d', float(text)).hex()


def argument(text):
    if text.startswith('{') and text.endswith('}'):
        return [[scalar(s) for s in split_arguments(row)]
                for row in split_arguments(text[1:-1], ';')]
    return scalar(text)


def outcome(value):
    if set(value) == {'number'}:
        return '0x'+value['number']
    if set(value) == {'error'}:
        # Capture and evaluation formats use different casing for the same
        # worksheet-error enum. Keep the original objects in each frozen row.
        names = {name.lower(): name for name in ('Null','Div0','Value','Ref','Name','Num','NA','GettingData','Spill','Calc','Blocked','Field','Unknown','Connect','Busy')}
        return 'error:'+names[value['error'].lower()]
    if set(value) == {'logical'}:
        return 'bool:'+str(value['logical']).lower()
    if set(value) == {'text'}:
        return 'text:'+value['text']
    raise ValueError('unsupported result encoding')


def classify(corpus, output, after=None, original_catalogue_after=None, scope_receipt=None):
    """Account for every failed formula without changing the verdict policy.

    Last-bits means two finite numeric values at most 16 binary64 steps
    apart. Larger differences are kept as value bugs, even if their relative
    error is small. Capture failures/refusals/missing results are separately
    named gross obligations; none are asserted to be proven product bugs.
    """
    corpus, output = Path(corpus), Path(output)
    compared = None
    if any(p is not None for p in (after, original_catalogue_after, scope_receipt)):
        if any(p is None for p in (after, original_catalogue_after, scope_receipt)):
            raise ValueError('after classification needs both comparison roots and the original-scope receipt')
        after, original_catalogue_after, scope_receipt = map(Path,
            (after, original_catalogue_after, scope_receipt))
        compared = json.loads(scope_receipt.read_text())
        if compared.get('complete') is not True:
            raise ValueError('original-scope comparison is incomplete')
    output.mkdir(parents=True, exist_ok=False)
    summary_path = corpus/'summary.json'
    summary = json.loads(summary_path.read_text())
    inputs = {str(summary_path): hashlib.sha256(summary_path.read_bytes()).hexdigest()}
    if compared is not None:
        if compared['counts']['formula_cells'] != summary['totals']['counts']['formula_cells']:
            raise ValueError('classification formula inventory differs from the comparison receipt')
        if compared['inputs_sha256'].get(str(summary_path)) != inputs[str(summary_path)]:
            raise ValueError('classification baseline differs from the comparison receipt')
        inputs[str(scope_receipt)] = hashlib.sha256(scope_receipt.read_bytes()).hexdigest()
    rows, counts, groups = [], Counter(), Counter()
    key = lambda r: (r['sheet'], r['row'], r['column'])
    runs = []
    pending_workbooks = []
    for workbook in summary['workbooks']:
        label = workbook['name'].replace('/', '_').replace('\\', '_')
        directory = corpus/'workbooks'/label
        captures = sorted(directory.glob('*.capture.json'))
        if captures:
            runs.extend((capture, directory/capture.name.removesuffix('.capture.json'), workbook['name'])
                        for capture in captures)
        elif workbook['capture'].startswith('committed capture '):
            runs.append((Path(workbook['capture'].removeprefix('committed capture ')), directory, workbook['name']))
        elif workbook['capture'].startswith('awaiting capture:'):
            c = workbook['counts']
            if c['awaiting_capture'] != c['formula_cells'] or c['pass'] or c['fail']:
                raise ValueError('inconsistent wholly pending workbook '+workbook['name'])
            pending_workbooks.append(dict(workbook=workbook['name'], cells=c['formula_cells'],
                                          source_sha256=workbook['sha256']))
        else:
            raise ValueError('unaccounted workbook '+workbook['name'])
    def number(bits):
        return struct.unpack('>d', bytes.fromhex(bits))[0]
    def ordered(bits):
        n = int(bits, 16)
        return (~n & ((1 << 64)-1)) if n >> 63 else n | (1 << 63)
    for capture_path, directory, workbook in runs:
        if compared is not None:
            relative = directory.relative_to(corpus)
            directory = (original_catalogue_after if 'calc-func_catalogue_oracle.xlsx' in relative.parts
                         else after/relative)
            for name in ('report.json', 'evaluations.json', 'triage.json'):
                path = directory/name
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                if compared['inputs_sha256'].get(str(path)) != digest:
                    raise ValueError('classification input differs from comparison receipt: '+str(path))
        paths = (capture_path, directory/'evaluations.json', directory/'triage.json')
        for p in paths:
            inputs[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
        capture, evaluations, triage = [json.loads(p.read_text()) for p in paths]
        expected = {key(r): r for r in capture['observations']}
        formulas = {key(r): r for r in capture['formulas']}
        actual = {key(r): r for r in evaluations['evaluations']}
        for verdict in triage['cells']:
            if verdict['verdict'] in ('matched', 'external_link', 'circular_rule', 'awaiting_capture'):
                continue
            if verdict['verdict'] not in ('mismatched', 'not_modeled', 'not_run', 'failed', 'missing', 'stale', 'unproven', 'file_cached'):
                raise ValueError('unknown nonpass '+str(verdict))
            k = key(verdict)
            formula = formulas.get(k, {}).get('formula', '')
            names = re.findall(r'([A-Za-z][A-Za-z0-9_.]*)\s*\(', formula)
            function = (names[1] if len(names) > 1 and names[0] == 'INDEX' else names[0]) if names else ('Infix &' if '&' in formula else formula)
            if verdict.get('class') == 'not-modelled':
                function = verdict['key']
            e, a = expected.get(k), actual.get(k)
            x = e.get('excel_raw') if e else None
            y = a.get('returned', {}).get('value') if a else None
            reason, distance = 'unavailable', None
            magnitude = dict(absolute_error=None, relative_error=None)
            if verdict.get('class') == 'capture':
                reason = 'capture_failure'
            elif verdict['verdict'] != 'mismatched':
                reason = verdict['verdict']
            elif isinstance(x, dict) and isinstance(y, dict):
                if set(x) == set(y) == {'number'}:
                    import math
                    if math.isfinite(number(x['number'])) and math.isfinite(number(y['number'])):
                        distance = abs(ordered(x['number'])-ordered(y['number']))
                        magnitude = error_size(x['number'], y['number'])
                        reason = 'last_bits' if distance <= 16 else 'numeric_value'
                    else:
                        reason = 'nonfinite'
                elif set(x) != set(y):
                    reason = 'type_or_error'
                elif set(x) == {'error'}:
                    reason = 'error_kind'
                elif set(x) == {'text'}:
                    reason = 'text_value'
                elif set(x) == {'array'}:
                    reason = 'shape_or_array_value'
                else:
                    reason = 'value'
            category = 'last_bits_only' if reason == 'last_bits' else 'gross_or_unavailable'
            counts[category] += 1
            groups[category, function, reason] += 1
            rows.append(dict(workbook=workbook, formula=formula, function=function,
                category=category, cause=reason, ulp_distance=distance,
                **magnitude,
                original_verdict=verdict, expected=e, actual=a,
                capture=str(capture_path), evaluation=str(paths[1])))
    require_count = compared['after_fail'] if compared is not None else summary['totals']['counts']['fail']
    if len(rows) != require_count:
        raise ValueError(f'nonpass accounting differs: {len(rows)} != {require_count}')
    (output/'cells.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    totals = summary['totals']['counts']
    matched = compared['counts']['after_matched'] if compared is not None else totals['pass']
    inclusive_failed = len(rows)+totals['awaiting_capture']+totals['external_links']+totals['circular_rule']
    if matched+inclusive_failed != totals['formula_cells']:
        raise ValueError('inclusive formula accounting differs')
    result = dict(complete=True, failed=len(rows), counts=dict(counts),
        matched=matched, formula_cells=totals['formula_cells'], inclusive_failed=inclusive_failed,
        pass_rate=matched/totals['formula_cells'] if totals['formula_cells'] else 0,
        awaiting_capture=totals['awaiting_capture'], wholly_pending_workbooks=pending_workbooks,
        groups=[dict(category=c, function=f, cause=r, cells=n)
                for (c,f,r), n in groups.most_common()], inputs_sha256=inputs,
        cells_sha256=hashlib.sha256((output/'cells.jsonl').read_bytes()).hexdigest(),
        rule='Last-bits-only: finite numbers <=16 ULP apart. Larger errors are a triage group, not a mathematical root-cause proof. '
             'Capture failures, refusals and missing values remain failures, explicitly distinguished from proven product bugs.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('inputs_sha256','groups')}))
    return result


def function_table(classification, output):
    """Reconcile every mismatch by function, retaining nonnumeric obligations."""
    source, output = Path(classification), Path(output)
    summary = json.loads((source/'summary.json').read_text())
    rows = [json.loads(line) for line in (source/'cells.jsonl').read_text().splitlines()]
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    if summary['complete'] is not True or len(rows) != summary['failed'] or sha(source/'cells.jsonl') != summary['cells_sha256']:
        raise ValueError('classification coverage or hash differs')
    groups, outside = {}, Counter()
    errors = {}
    for row in rows:
        verdict = row['original_verdict']['verdict']
        if verdict != 'mismatched':
            outside[verdict] += 1
            continue
        counts = groups.setdefault(row['function'], Counter())
        counts[row['cause']] += 1
        if row['ulp_distance'] is not None:
            errors.setdefault(row['function'], []).append(row)
    functions = [dict(function=f, mismatches=sum(c.values()),
                      numeric=c['last_bits']+c['numeric_value'], within_16_ulp=c['last_bits'],
                      larger_numeric=c['numeric_value'],
                      other={k:n for k,n in sorted(c.items()) if k not in ('last_bits','numeric_value')})
                 for f,c in groups.items()]
    functions.sort(key=lambda r: (-r['mismatches'], r['function']))
    for function in functions:
        numeric = errors.get(function['function'], [])
        ulps = sorted(r['ulp_distance'] for r in numeric)
        relatives = sorted(Decimal(r['relative_error']) for r in numeric if r.get('relative_error') is not None)
        function['gross'] = function['mismatches']-function['within_16_ulp']
        function['ulp_min'] = ulps[0] if ulps else None
        function['ulp_median_upper'] = ulps[len(ulps)//2] if ulps else None
        function['ulp_max'] = ulps[-1] if ulps else None
        function['relative_error_max'] = str(relatives[-1]) if relatives else None
        function['relative_error_undefined'] = sum(r.get('relative_error') is None for r in numeric)
    output.mkdir(parents=True, exist_ok=False)
    result = dict(complete=True, mismatches=sum(r['mismatches'] for r in functions),
        nonpasses_in_classification=len(rows), other_verdicts=dict(outside), functions=functions,
        inputs_sha256={str(source/name):sha(source/name) for name in ('cells.jsonl','summary.json')},
        limits='Exactly the mismatched verdicts, including capture failures. Refusals and missing cells remain failures outside this table. '
               'External and circular-policy nonpasses are outside the classifier; the corpus denominator includes them. '
               'The 16-ULP split is triage, not an algorithm or accuracy concession.')
    if result['mismatches']+sum(outside.values()) != len(rows):
        raise ValueError('function totals do not reconcile')
    (output/'table.json').write_text(json.dumps(result, indent=2)+'\n')
    columns = ('function','mismatches','numeric','within_16_ulp','larger_numeric','gross',
               'ulp_min','ulp_median_upper','ulp_max','relative_error_max','relative_error_undefined')
    lines = ['\t'.join(columns)+'\tother']
    for r in functions:
        lines.append('\t'.join(str(r[k]) for k in columns)+
                     '\t'+json.dumps(r['other'],sort_keys=True))
    (output/'table.tsv').write_text('\n'.join(lines)+'\n')
    columns = ('function','workbook','sheet','row','column','category','cause','ulp_distance',
               'absolute_error','relative_error','expected','actual','formula')
    lines = ['\t'.join(columns)]
    for row in rows:
        if row['original_verdict']['verdict'] != 'mismatched':
            continue
        cell = dict(row, **{k:row['original_verdict'][k] for k in ('sheet','row','column')})
        cell['expected'] = row['expected'].get('excel_raw') if row['expected'] else None
        cell['actual'] = row['actual'].get('returned', {}).get('value') if row['actual'] else None
        lines.append('\t'.join(json.dumps(cell.get(k),ensure_ascii=True) for k in columns))
    (output/'mismatches.tsv').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'functions'}))
    return result


def freeze(corpus, output, functions):
    corpus, output = Path(corpus), Path(output)
    output.mkdir()
    wanted = set(functions.split(','))
    counts = Counter()
    inputs = {}
    runs = [(capture, capture.parent/capture.name.removesuffix('.capture.json')/'evaluations.json', capture.parent.name)
            for capture in sorted(corpus.glob('workbooks/*/*.capture.json'))]
    # Sectioned catalogues have merged capture copies in the corpus output.
    # Ordinary oracles instead name their committed capture in summary.json.
    summary_path=corpus/'summary.json'
    if summary_path.is_file():
        inputs[str(summary_path)]=hashlib.sha256(summary_path.read_bytes()).hexdigest()
        for workbook in json.loads(summary_path.read_text())['workbooks']:
            source=workbook['capture']
            if not source.startswith('committed capture '):continue
            capture=Path(source.removeprefix('committed capture '))
            label=workbook['name'].replace('/','_').replace('\\','_')
            evaluation=corpus/'workbooks'/label/'evaluations.json'
            if not evaluation.is_file():raise ValueError('missing corpus evaluation '+str(evaluation))
            runs.append((capture,evaluation,label))
    with (output/'cases.jsonl').open('x') as out, (output/'excluded.jsonl').open('x') as excluded:
        for capture_path,evaluation_path,label in runs:
            capture = json.loads(capture_path.read_text())
            triage_path = evaluation_path.with_name('triage.json')
            evaluations = json.loads(evaluation_path.read_text())
            triage = json.loads(triage_path.read_text())
            key = lambda r: (r['sheet'], r['row'], r['column'])
            expected = {key(r): r for r in capture['observations']}
            actual = {key(r): r for r in evaluations['evaluations']}
            verdicts = {key(r): r for r in triage['cells']}
            for f in capture['formulas']:
                formula = f['formula'].lstrip('=')
                selected = None
                if formula.startswith('INDEX(') and formula.endswith(')'):
                    parts = split_arguments(formula[6:-1])
                    if len(parts) == 3 and all(re.fullmatch(r'\d+', x) for x in parts[1:]):
                        formula = parts[0]
                        selected = list(map(int, parts[1:]))
                m = re.fullmatch(r'(@)?(?:_xlfn\.)?([A-Za-z.]+)\((.*)\)', formula)
                if not m or m[2].upper() not in wanted:
                    continue
                function = m[2].upper()
                e, a, v = expected[key(f)], actual[key(f)], verdicts[key(f)]
                raw = dict(id=f"corpus:{label}:{f['sheet']}:{f['row']}:{f['column']}",
                           function=function, result_index=selected, formula=f['formula'],
                           source=str(capture_path), reading=capture['environment'],
                           original_capture=e, original_evaluation=a,
                           original_formula=f, verdict=v['verdict'])
                counts[function+':'+v['verdict']] += 1
                try:
                    if m[1]:raise ValueError('implicit intersection surface outside adapter')
                    sources = split_arguments(m[3])
                    raw['argument_sources'] = sources
                    raw['args'] = list(map(argument, sources))
                    raw['expected'] = outcome(e['excel_raw'])
                    raw['actual'] = outcome(a['returned']['value'])
                except (KeyError, ValueError) as error:
                    excluded.write(json.dumps(dict(raw=raw, reason=str(error)))+'\n')
                    counts['excluded'] += 1
                    continue
                out.write(json.dumps(raw)+'\n')
                counts['frozen'] += 1
            for p in (capture_path, evaluation_path, triage_path):
                inputs[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
    result = dict(counts=counts, inputs=inputs,
                  condition='Typed arguments require application parser/workbook audit; no new readings or expected values.')
    with (output/'summary.json').open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(dict(counts=counts)))


def variance_questions(input_path, output):
    """Freeze scalar transport/component audits after the variance exclusions.

    Original captures use arrays; these requests deliberately use three direct
    Value2 numbers. They do not assert typed-input identity with an array call.
    Adjacent inputs are fixed one-ULP perturbations, never chosen by residuals.
    """
    from fractions import Fraction
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    original = [json.loads(line) for line in Path(input_path).read_text().splitlines()]
    def flatten(value):
        return [v for item in value for v in flatten(item)] if isinstance(value, list) else [value]
    selected = []
    seen = set()
    for row in original:
        args = tuple(flatten(row['args']))
        if args not in seen:
            selected.append((args, row))
            seen.add(args)
        if len(selected) == 6:
            break
    if len(selected) != 6:
        raise ValueError('six original triples are required')
    names = ('VAR.S', 'VAR.P', 'STDEV.S', 'STDEV.P', 'DEVSQ', 'AVERAGE')
    requests = {name: [] for name in names}
    selection = []
    controls = []
    def add(function, args, source, role):
        identity = f'r14-{function.lower().replace(".","-")}-{len(requests[function]):02d}'
        requests[function].append(dict(id=identity, args=list(args), input_transport='cell_value2_bulk'))
        selection.append(dict(id=identity, function=function, role=role, source_observation=source,
                              requested_args=list(args),
                              typed_surface='Three direct numeric Value2 arguments; source used one numeric array.',
                              exposure='All old corpus/bank observations are discovery. No complete candidate survives.'))
    tuples = []
    for index, (args, row) in enumerate(selected):
        for function in names:
            add(function, args, row, f'original numeric tuple {index}: surface/component audit')
        tuples.append((args, 'original numeric tuple'))
    for index, (args, row) in enumerate(selected[:4]):
        adjacent = list(args)
        bits = int(adjacent[0], 16)
        bits += -1 if bits >> 63 else 1
        adjacent[0] = f'0x{bits:016x}'
        for function in ('VAR.S', 'DEVSQ'):
            add(function, adjacent, row, f'tuple {index}: first operand one ULP toward positive infinity')
        tuples.append((adjacent, 'new adjacent numeric tuple; no Excel reading'))
    for args, role in tuples:
        values = [Fraction(struct.unpack('>d', bytes.fromhex(x[2:]))[0]) for x in args]
        mean = sum(values)/3
        devsq = sum((x-mean)**2 for x in values)
        encode = lambda x: [str(x.numerator), str(x.denominator)]
        controls.append(dict(args=list(args), role=role, exact_mean=encode(mean),
                             exact_devsq=encode(devsq), exact_sample_variance=encode(devsq/2),
                             exact_population_variance=encode(devsq/3),
                             not_excel_answers=True))
    for function, witnesses in requests.items():
        with (output/f'requests-{function}.json').open('x') as out:
            json.dump(dict(function=function, witnesses=witnesses), out, indent=2)
            out.write('\n')
    for name, data in [('selection.json', selection), ('mathematical-controls.json', controls)]:
        with (output/name).open('x') as out:
            json.dump(data, out, indent=2)
            out.write('\n')
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(output.iterdir())}
    manifest = dict(format='dnacalc-research-questions-r14-variance/1', requests=sum(map(len, requests.values())),
                    counts={k: len(v) for k,v in requests.items()}, input_sha256=hashlib.sha256(Path(input_path).read_bytes()).hexdigest(),
                    files=files, new_excel_readings=0,
                    purpose='Distinguish direct numeric transport from the original array surface and observe mean/DEVSQ companions. No complete numerical program is validated.',
                    predictions='Application predictions are separate files, not recorded Excel answers. Exact rational mathematical controls are also separate.',
                    commands=[f'cargo xtask witness read --func {f} --inputs requests-{f}.json' for f in names],
                    limits='Do not use --new-only. This is a mixed recheck/adjacent-input packet. Refused, missing or changed-input responses remain failed obligations. No annuity request is added.')
    with (output/'manifest.json').open('x') as out:
        json.dump(manifest, out, indent=2)
        out.write('\n')
    print(json.dumps(dict(requests=manifest['requests'], counts=manifest['counts'])))


def variance_scope(bank, corpus, observations, output, names=None):
    """Keep all discovery constraints, including duplicate and differing surfaces."""
    from collections import defaultdict
    functions = {'VAR', 'VAR.S', 'VAR.P', 'VARP', 'VARA', 'VARPA',
                 'STDEV', 'STDEV.S', 'STDEV.P', 'STDEVP', 'STDEVA', 'STDEVPA',
                 'DEVSQ', 'AVERAGE'}
    if names is not None:
        functions = set(names.split(','))
    output = Path(output)
    output.mkdir()
    paths = list(map(Path, (observations, corpus, bank)))
    counts = Counter()
    by_flat_input = defaultdict(list)
    def flatten(value):
        return [v for item in value for v in flatten(item)] if isinstance(value, list) else [value]
    sources = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths[:2]}
    with (output/'constraints.jsonl').open('x') as out:
        for path in paths:
            hasher = hashlib.sha256()
            with path.open('rb') as stream:
                for line in stream:
                    hasher.update(line)
                    row = json.loads(line)
                    if row['function'] not in functions:
                        continue
                    surface = 'array_argument' if any(isinstance(a, list) for a in row['args']) else 'separate_arguments'
                    out.write(json.dumps(dict(row, research_argument_surface=surface))+'\n')
                    counts[row['function']+':'+surface] += 1
                    values = flatten(row['args'])
                    if all(isinstance(x, str) and x.startswith('0x') and len(x) == 18 for x in values):
                        by_flat_input[(row['function'], tuple(values), tuple(row.get('result_index') or ()))].append(dict(
                            id=row['id'], surface=surface, expected=row['expected'],
                            source=row['source'], reading=row['reading']))
            sources[str(path)] = hasher.hexdigest()
    contradictions = []
    for (function, values, selected), rows in by_flat_input.items():
        if len({r['expected'] for r in rows}) > 1:
            contradictions.append(dict(function=function, flattened_input=list(values), result_index=list(selected) or None, observations=rows,
                interpretation='Different original surfaces/environments are retained; flattening is not proof of identical Excel calls. A surface-independent program must satisfy every row.'))
    (output/'different-answers.json').write_text(json.dumps(contradictions, indent=2)+'\n')
    result = dict(complete=True, counts=counts, different_flattened_answers=len(contradictions),
                  input_sha256=sources, source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  limits='All selected canonical full-bank calls, all parsable corpus calls and all R14 observations, including repeats. Array shape remains in args; scalar/reference provenance remains in raw metadata. Nonnumeric, domain and unparsed calls remain obligations outside the finite grammar. No surface switch is inferred.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


def static_reference_cases(corpus, package, output):
    """Extend the corpus adapter for the measured GEOMEAN Data!A1:A3 references.

    Keep the original reference formula. The numeric array is a proposed input
    to arithmetic models, justified by exact static package cells and the
    separately tested application collector, never a new scalar Excel reading.
    """
    import zipfile
    import xml.etree.ElementTree as ET
    from fractions import Fraction
    corpus,package,output=map(Path,(corpus,package,output));output.mkdir()
    ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    rid='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id'
    with zipfile.ZipFile(package)as z:
        rel={r.attrib['Id']:r.attrib['Target']for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets={s.attrib['name']:'xl/'+rel[s.attrib[rid]]for s in ET.fromstring(z.read('xl/workbook.xml')).find('m:sheets',ns)}
        data={c.attrib['r']:c for c in ET.fromstring(z.read(sheets['Data'])).findall('.//m:c',ns)}
        values=[]
        for name in('A1','A2','A3'):
            cell=data[name]
            if cell.find('m:f',ns)is not None or cell.attrib.get('t','n')!='n':raise ValueError('reference input is not a static number')
            text=cell.find('m:v',ns).text;value=float(text)
            if Fraction(value)!=Fraction(text):raise ValueError('reference input is not exactly representable')
            values.append('0x'+struct.pack('>d',value).hex())
    rows=[];inputs={str(package):hashlib.sha256(package.read_bytes()).hexdigest()}
    for capture in corpus.glob('workbooks/*/*.capture.json'):
        doc=json.loads(capture.read_text())
        evaluated=capture.parent/capture.name.removesuffix('.capture.json')/'evaluations.json'
        observations={(r['sheet'],r['row'],r['column']):r for r in doc['observations']}
        actual={(r['sheet'],r['row'],r['column']):r for r in json.loads(evaluated.read_text())['evaluations']}
        for formula in doc['formulas']:
            if formula['formula'].lstrip('=')not in('GEOMEAN(Data!A1:A3)','_xlfn.GEOMEAN(Data!A1:A3)'):continue
            package_source='sha256:'+inputs[str(package)]
            if doc['source']!=package_source or doc['session']['packages']['input']!=package_source:
                raise ValueError('reference capture/package lineage differs')
            key=(formula['sheet'],formula['row'],formula['column']);e=observations[key];a=actual[key]
            rows.append(dict(id=f"corpus:{capture.parent.name}:{key[0]}:{key[1]}:{key[2]}",function='GEOMEAN',
                result_index=None,formula=formula['formula'],argument_sources=['Data!A1:A3'],
                args=[[[v]for v in values]],expected=outcome(e['excel_raw']),actual=outcome(a['returned']['value']),
                source=str(capture),reading=doc['environment'],original_formula=formula,original_capture=e,original_evaluation=a,
                reference_binding=dict(package=str(package),package_sha256=inputs[str(package)],sheet='Data',cells=['A1','A2','A3'],
                    values=values,capture_sha256=hashlib.sha256(capture.read_bytes()).hexdigest()),
                typed_surface='Original static reference; proposed numeric operands for arithmetic models, not a scalar/array Excel re-reading.'))
            inputs[str(capture)]=hashlib.sha256(capture.read_bytes()).hexdigest()
            inputs[str(evaluated)]=hashlib.sha256(evaluated.read_bytes()).hexdigest()
    if len(rows)!=2:raise ValueError('expected the two recorded GEOMEAN reference obligations')
    (output/'cases.jsonl').write_text(''.join(json.dumps(r)+'\n'for r in rows))
    result=dict(complete=True,reference_observations=2,input_sha256=inputs,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limits='Recorded reference calls and exact static package cells only. Proposed numeric operands retain their original reference provenance. No new Excel answer, changed expectation, or product pass is created.')
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


def round16_questions(readings, executable, output, section):
    """Freeze surface, domain and rank questions; keep every answer field absent."""
    import itertools
    import math
    import subprocess
    readings, output, section = map(Path, (readings, output, section))
    output.mkdir()
    def bits(x): return '0x'+struct.pack('>d', float(x)).hex()
    def number(x): return struct.unpack('>d', bytes.fromhex(x[2:]))[0]
    def literal(x): return format(float(x), '.17g')
    requests, selection, sheet = {}, [], []
    seen = set()
    def ask(function, values, reason, origin=None):
        args = [bits(x) for x in values]
        key = function, tuple(args)
        if key in seen: return
        seen.add(key)
        group = requests.setdefault(function, [])
        qid = 'r16-'+function.lower().replace('.', '-')+f'-{len(group):03d}'
        group.append(dict(id=qid, args=args, input_transport='cell_value2_bulk'))
        selection.append(dict(id=qid, function=function, reason=reason, source_id=origin,
                              requested_args=args, surface='separate numeric Value2 references'))
    def formula(text, reason, **metadata):
        row = len(sheet)+1
        sheet.append(dict(id=f'r16-sheet-{row:04d}', row=row, column=4,
                          formula=text, reason=reason, **metadata))
        return f'D{row}'
    def surfaces(function, values, reason):
        literals = list(map(literal, values))
        refs = [formula(x, reason+': operand echo', input_bits=bits(v))
                for x, v in zip(literals, values)]
        cells = f'{refs[0]}:{refs[-1]}'
        for body, kind in [(','.join(refs), 'separate references'),
                (cells, 'vertical reference'),
                ('{'+','.join(literals)+'}', 'horizontal array literal'),
                ('{'+ ';'.join(literals)+'}', 'vertical array literal'),
                (f'TRANSPOSE({cells})', 'transposed reference'),
                (','.join(literals), 'separate literals')]:
            formula(f'{function}({body})', reason, function=function, argument_surface=kind,
                    intended_input_bits=list(map(bits, values)), operand_echo_cells=refs)
    docs = {f: json.loads((readings/(f+'.json')).read_text()) for f in
            ('VAR.S', 'VAR.P', 'STDEV.S', 'STDEV.P', 'DEVSQ', 'AVERAGE')}
    triples = [list(map(number, w['args'])) for w in docs['VAR.S']['witnesses']]
    names = list(docs)
    for index, values in enumerate(triples):
        for function in names:
            ask(function, values, 'R14 tuple; complete component recheck', docs['VAR.S']['witnesses'][index]['id'])
        refs = [formula(literal(v), 'variance operand echo', input_bits=bits(v), tuple_index=index) for v in values]
        for function in names if index < 6 else ('VAR.S', 'DEVSQ'):
            literals = list(map(literal, values))
            for body, surface in [(','.join(refs), 'separate references'),
                    (f'{refs[0]}:{refs[-1]}', 'vertical reference'),
                    ('{'+','.join(literals)+'}', 'horizontal array literal'),
                    ('{'+ ';'.join(literals)+'}', 'vertical array literal'),
                    (f'TRANSPOSE({refs[0]}:{refs[-1]})', 'transposed reference'),
                    (','.join(literals), 'separate literals')]:
                formula(f'{function}({body})', 'variance surface/component audit', function=function,
                        argument_surface=surface, tuple_index=index, intended_input_bits=list(map(bits, values)),
                        operand_echo_cells=refs)
    limit = float(2**27)
    offsets = (-2, -1, -.75, -.5, -.25, 0, .25, .5, .75, 1, 2)
    xs = [limit+v for v in offsets]+[math.nextafter(limit, 0), math.nextafter(limit, math.inf)]
    for x, order in itertools.product(xs, (0, 1, 2, 3, 4)):
        ask('BESSELJ', [x, order], 'fixed binary phase-boundary grid; no fitted cutoff')
    for x, order in itertools.product((limit-1, limit-.5, limit, limit+.5), (0, 1, 2)):
        ask('BESSELJ', [-x, order], 'negative argument parity/domain control')
    for x, order in itertools.product((limit-1, limit-.5, math.nextafter(limit, 0), limit, math.nextafter(limit, math.inf), limit+1), (0, 1, 2)):
        ask('BESSELY', [x, order], 'same-domain first/second-kind comparison')
    for x, order in itertools.product((709, 709.5, math.nextafter(709.7, 0), 709.7,
            math.nextafter(709.7, math.inf), 709.78, 710, 713, 714), (0, 1, 2, 3, 10, 100)):
        ask('BESSELI', [x, order], 'fixed exponent-domain grid across orders; no inferred admission rule')
    for x, order in itertools.product((709.7, math.nextafter(709.7, math.inf), 710), (0, 1, 2)):
        ask('BESSELI', [-x, order], 'negative argument parity/domain control')
    for x in (709, 709.5, 709.7, math.nextafter(709.7, math.inf), 709.78,
              math.log(float.fromhex('0x1.fffffffffffffp+1023')), 710):
        ask('EXP', [x], 'independent exponential admission; not a hidden Bessel leaf')
    for x in (limit-1, limit-.5, math.nextafter(limit, 0), limit, math.nextafter(limit, math.inf), limit+.5):
        ref = formula(literal(x), 'phase input echo', input_bits=bits(x))
        for function in ('SIN', 'COS'):
            ask(function, [x], 'independent trigonometric admission; no Bessel call identity assumed')
        for phase in (1, 3):
            p = formula(f'{ref}-{phase}*PI()/4', 'published Bessel asymptotic phase; worksheet-computed echo')
            for function in ('SIN', 'COS'):
                formula(f'{function}({p})', 'stored phase trigonometric probe', argument_surface='reference')
                formula(f'{function}({ref}-{phase}*PI()/4)', 'expression phase trigonometric probe', argument_surface='expression')
    y = '{1;2;4;7}'
    def matrix(delta, scale1, scale2, affine=False):
        rows = [[x*scale1, (2*x+(1 if affine else 0)+(delta if i==3 else 0))*scale2]
                for i, x in enumerate((1, 2, 5, 4))]
        return '{'+';'.join(','.join(literal(v) for v in r) for r in rows)+'}'
    matrices = set()
    def regression(delta, scale1=1, scale2=1, affine=False):
        x = matrix(delta, scale1, scale2, affine)
        for constant, function in itertools.product(('FALSE', 'TRUE'), ('LINEST', 'LOGEST')):
            for r, c in ((1, 1), (1, 2), (4, 2)):
                text = f'INDEX({function}({y},{x},{constant},TRUE),{r},{c})'
                if text in matrices: continue
                matrices.add(text)
                formula(text, 'rank and removed-column audit', function=function,
                        delta=literal(delta), scale1=literal(scale1), scale2=literal(scale2),
                        use_const=constant=='TRUE', affine=affine, result_index=[r,c])
    for delta in (0, 2**-30, 2**-26, 2**-24, 2**-22, 2**-20, 1e-6, 2**-18, 2**-16, 2**-12):
        regression(delta)
    for delta, scales in itertools.product((0, 2**-24, 2**-20, 2**-16),
            ((2**-20, 2**-20), (2**20, 2**20), (1, 2**-20), (1, 2**20))):
        regression(delta, *scales)
    regression(0, affine=True)
    # Unperturbed affine and multi-column dependencies are structural controls.
    exact_matrices = ('{1,3;2,5;5,11;4,9}', '{1,0,1;0,1,1;1,1,2;2,1,3}', '{1,0;2,0;5,0;4,0}')
    for x, constant, function in itertools.product(exact_matrices, ('FALSE', 'TRUE'), ('LINEST', 'LOGEST')):
        formula(f'INDEX({function}({y},{x},{constant},TRUE),4,2)',
                'exact dependency/zero-column degrees-of-freedom control', function=function,
                use_const=constant=='TRUE', result_index=[4,2])
    for values in list(itertools.permutations((1, 5, 3)))+[(15,), (1,15), (3,5), (15,1,1), (1,15,1), (1,1,15), (1,3,5,1)]:
        ask('GEOMEAN', values, 'permutation/count/product control for the two reference misses')
        surfaces('GEOMEAN', values, 'GEOMEAN reference/array/literal audit')
    batches = []
    executable = str(Path(executable).resolve())
    for function, rows in requests.items():
        request = dict(function=function, selection='R16 fixed surface/domain questions; intentional rechecks, no --new-only', witnesses=rows)
        name = f'requests-{function}.json'
        (output/name).write_text(json.dumps(request, indent=2)+'\n')
        # ubs:ignore — the operator supplies the trusted local executable; arguments never enter a shell.
        subprocess.run([executable, 'invoke-inputs', str(output/name), str(output/f'application-predictions-{function}.jsonl')], check=True, timeout=300)
        batches.append(dict(function=function, path=name, calls=len(rows)))
    (output/'batches.json').write_text(json.dumps(batches, indent=2)+'\n')
    (output/'selection.json').write_text(json.dumps(selection, indent=2)+'\n')
    (output/'catalogue-questions.json').write_text(json.dumps(dict(sheet='dcj9k_r16', questions=sheet), indent=2)+'\n')
    text = '# R16: surface, domain and rank questions. No expected answers.\n'+'\n'.join(q['formula'] for q in sheet)+'\n'
    section.write_text(text, encoding='ascii')
    # This is a mirror of question text, not a second workbook or capture lineage.
    (output/section.name).write_text(text, encoding='ascii')
    result = dict(format='dnacalc-research-questions-r16/1', scalar_calls=sum(len(x) for x in requests.values()),
        calls_by_function={f:len(rows) for f,rows in requests.items()}, catalogue_formula_cells=len(sheet),
        catalogue_sheet='dcj9k_r16', catalogue_section=str(section.resolve()), executable_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
        sources_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(readings.glob('*.json'))},
        new_excel_answers=0, limits='Questions only. Scalar Value2 references use CompatibilityVersion/settings recorded by the standard reader; catalogue questions use the existing section capture. Literal operand echoes expose worksheet parsing, not private callee inputs. No fitted cutoff/rank tolerance or validated numerical candidate. Preserve all settings, failures and original raw outputs. No R15 annuity question or answer.')
    (output/'question-summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))


def par_yield_questions(executable, output, bank_trace, text_failures=None):
    """Par identity and recorded text ties; extend the standard scalar packet."""
    import datetime
    import math
    import subprocess
    executable, output, bank_trace = Path(executable).resolve(), Path(output), Path(bank_trace)
    output.mkdir()
    def bits(value):
        return '0x'+struct.pack('>d', float(value)).hex()
    def date(year, month, day):
        return float((datetime.date(year, month, day)-datetime.date(1899, 12, 30)).days)
    bank = [json.loads(line) for line in bank_trace.read_text().splitlines()]
    known = {json.dumps(row['args'], separators=(',', ':')): row['id'] for row in bank}
    questions, selection, identities = [], [], set()
    def ask(settlement, maturity, rate, price, frequency, basis, reason):
        args = list(map(bits, (settlement, maturity, rate, price, 100.0, frequency, basis)))
        identity = json.dumps(args, separators=(',', ':'))
        require = identity not in identities
        if not require:
            raise ValueError('duplicate par-yield question')
        identities.add(identity)
        qid = f'r18-yield-par-{len(questions):03d}'
        questions.append(dict(id=qid, args=args, input_transport='cell_value2_bulk'))
        selection.append(dict(id=qid, reason=reason, existing_bank_id=known.get(identity),
                              role='recorded recheck' if identity in known else 'new typed identity'))
    shapes = [(1, date(2021, 1, 1), date(2024, 1, 1)),
              (2, date(2020, 7, 1), date(2022, 1, 1)),
              (4, date(2020, 10, 1), date(2022, 1, 1))]
    for frequency, settlement, maturity in shapes:
        for rate in (0.0, 0.05, 0.0575, 0.125):
            ask(settlement, maturity, rate, 100.0, frequency, 0, 'multi-period on-coupon mathematical par identity')
            ask(settlement, maturity, rate, math.nextafter(100.0, math.inf), frequency, 0,
                'one binary64 price step above par; independent boundary control')
    for basis in (2, 3):
        ask(date(2020, 7, 1), date(2022, 1, 1), 0.0575, 100.0, 2, basis,
            'multi-period par identity with Actual/360 or Actual/365')
        ask(date(2020, 7, 1), date(2021, 1, 1), 0.05, 100.0, 2, basis,
            'single-period documented day-count formula; excluded from the identity')
    ask(date(2020, 7, 1), date(2028, 1, 1), 0.05, 100.0, 2, 0,
        'recorded long-par recheck')
    if len(questions) != 29:
        raise ValueError('incomplete par-yield question plan')
    batches = {'YIELD': questions}
    if text_failures is not None:
        from number_conversion import observations
        failures = [json.loads(line) for line in Path(text_failures).read_text().splitlines()]
        selected = {(Path(r['workbook']).parts[1], r['expected']['sheet'], r['expected']['row'], r['expected']['column'])
                    for r in failures if r['cause'] == 'text_value' and r['function'] == 'Infix &'}
        originals = {}
        for row in observations(bank_trace.resolve().parents[2]):
            if (row['case'], row['sheet'], row['row'], row['column']) in selected:
                originals.setdefault(row['input'], []).append({k:v for k,v in row.items() if k not in ('exact_input','expected_decimal')})
        if len(originals) != 5:
            raise ValueError('recorded five text-tie operands changed')
        text_questions = []
        for source, bindings in sorted(originals.items()):
            center = float.fromhex(source)
            for value, direction in ((math.nextafter(center, -math.inf), 'below'), (center, 'center'),
                                     (math.nextafter(center, math.inf), 'above')):
                qid = f'r18-concat-tie-{len(text_questions):03d}'
                text_questions.append(dict(id=qid, args=[bits(value)], input_transport='cell_value2_bulk'))
                selection.append(dict(id=qid, reason='CONCAT numeric coercion with exact Value2 readback, '+direction,
                    original_infix_observations=bindings, role='new coercion-surface question; original infix answers are discovery'))
        batches['CONCAT'] = text_questions
    audit_input = output/'application-route-inputs.jsonl'
    with audit_input.open('x') as out:
        for function, batch in batches.items():
            request = dict(function=function, selection='R18 mathematical/boundary questions; intentional rechecks, no --new-only', witnesses=batch)
            request_path = output/('requests-'+function+'.json')
            request_path.write_text(json.dumps(request, indent=2)+'\n')
            prediction_path = output/('application-predictions-'+function+'.jsonl')
            # ubs:ignore — operator-selected local executable and argv, no shell; predictions are explicitly not Excel evidence.
            subprocess.run([str(executable), 'invoke-inputs', str(request_path), str(prediction_path), '1900'], check=True, timeout=300)
            predictions = [json.loads(line) for line in prediction_path.read_text().splitlines()]
            if len(predictions) != len(batch):
                raise ValueError('missing packet predictions')
            for question, prediction in zip(batch, predictions):
                if question['id'] != prediction['id'] or question['args'] != prediction['args'] or prediction.get('date_system') != '1900':
                    raise ValueError('packet prediction input/environment differs')
                out.write(json.dumps(dict(**question, function=function, result_index=None,
                    actual=prediction['calc_func_prediction'], numeric_cell_inputs=True,
                    date_system_1900=True, not_excel_evidence=True))+'\n')
    # ubs:ignore — same selected local executable, argv only; a failed application audit aborts packet publication.
    subprocess.run([str(executable), 'function-routes', str(audit_input), str(output/'application-route-audit.jsonl')], check=True, timeout=300)
    (output/'selection.json').write_text(json.dumps(selection, indent=2)+'\n')
    (output/'batches.json').write_text(json.dumps([dict(function=f, path='requests-'+f+'.json', calls=len(q)) for f,q in batches.items()], indent=2)+'\n')
    (output/'README.txt').write_text('''R18 scalar questions, no catalogue additions; counts in batches.json.

From the published DnaCalc checkout on Windows, use the existing reader:
  cargo xtask witness read --func YIELD --inputs <packet>/requests-YIELD.json
  cargo xtask witness read --func CONCAT --inputs <packet>/requests-CONCAT.json
Do not use --new-only. The packet deliberately rechecks two recorded par calls.
Only requests-*.json are inputs to the reader. Predictions and route audits
are application calculations, never Excel readings or expected-answer files.

These calls ask whether Excel follows the exact mathematical par identity over
more coupon rates/frequencies, and retain adjacent-price and single-period
day-count controls. Mathematics alone does not settle Excel's iteration bits.
No fitted stopping tolerance, cutoff, residual selector or seed is requested.
The optional CONCAT batch asks the five unresolved infix number-text ties and
their adjacent binary64 values through numeric-cell coercion. Exact Value2
readback separates submitted operands from authored worksheet literals. CONCAT
is a new surface; matching it would not by itself prove the infix route fixed.
The Linux predictions declare date system 1900; keep the reader's complete
Excel environment, requested/submitted/Value2 operand audits and failures.

Return the canonical new WitnessSet under OxFunc-Evidence/readings, the new
target/witness/read run directory and target/witness/frozen files, plus this
unchanged packet and file hashes. Preserve partial attempts. Missing, unread
or changed-input obligations fail. Use the existing capture lineage only.
''')
    sha = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
    summary = dict(format='dnacalc-research-questions-r18/1', scalar_calls=sum(map(len, batches.values())),
        calls_by_function={f:len(q) for f,q in batches.items()}, catalogue_formula_cells=0,
        role='implementer', new_excel_answers=0, executable_sha256=sha(executable),
        sources_sha256={str(p):sha(p) for p in (Path(__file__), bank_trace,
            *([Path(text_failures), Path(__file__).with_name('number_conversion.py')] if text_failures else []))},
        exposure='Predictions frozen before new answers. Two recorded par rechecks are discovery; no full YIELD algorithm or held-out performance claim.',
        limits='Tests the mathematical par identity, neighboring prices, the documented single-period day-count exception and optional text-coercion ties. Ordinary product routes only; standalone calls do not expose Excel private solver or conversion state.')
    summary['files_sha256'] = {p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'manifest.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if 'sha256' not in k}))


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'round22-qr-questions':
        round22_qr_questions(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'round22-questions':
        round22_questions(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'round21-questions':
        round21_questions(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'hypergeom-symmetry':
        hypergeom_symmetry(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'discrete-questions':
        discrete_questions(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'domain-orders':
        domain_orders(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'function-table':
        function_table(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'par-yield-questions':
        par_yield_questions(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'classify':
        classify(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'static-reference-cases':
        static_reference_cases(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'round16-questions':
        round16_questions(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'variance-scope':
        variance_scope(*sys.argv[2:])
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == 'variance-questions':
        parser = argparse.ArgumentParser()
        parser.add_argument('command')
        parser.add_argument('input')
        parser.add_argument('output')
        args = parser.parse_args()
        variance_questions(args.input, args.output)
        sys.exit(0)
    parser = argparse.ArgumentParser()
    parser.add_argument('corpus')
    parser.add_argument('output')
    parser.add_argument('functions')
    args = parser.parse_args()
    freeze(args.corpus, args.output, args.functions)
