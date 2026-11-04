"""Bounded, deterministic questions from actual local corpus calls.

No expected answer is synthesized. The existing corpus freezer, workbook
generator, application parser and Excel capture command remain the only routes.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import subprocess
import tempfile

if __package__:
    from .corpus_cases import freeze
else:
    from corpus_cases import freeze


# Owner's parked numerical producers. They remain in the problem report.
PARKED = {'PMT', 'IPMT', 'PPMT', 'CUMIPMT', 'CUMPRINC', 'FV', 'PV', 'NPER',
          'VAR', 'VAR.S', 'VARP', 'VAR.P', 'VARA', 'VARPA', 'STDEV', 'STDEV.S',
          'STDEVP', 'STDEV.P', 'STDEVA', 'STDEVPA', 'GEOMEAN'}


def bits(value):
    return '0x'+struct.pack('>d', value).hex()


def number(value):
    return struct.unpack('>d', bytes.fromhex(value[2:]))[0]


def exact_expression(value):
    """Two small integer limbs avoid decimal-literal collapse of adjacent inputs.

    Powers and integer limbs are exact. Subnormal outcomes deliberately probe
    the application's/Excel's underflow policies and require operand echoes.
    """
    if not math.isfinite(value):
        raise ValueError('finite operand required')
    if value == 0:
        return '-0' if math.copysign(1, value) < 0 else '0'
    numerator, denominator = abs(value).as_integer_ratio()
    exponent = -(denominator.bit_length()-1)
    while numerator >= 2**53:
        if numerator & 1:
            raise ValueError('nonbinary64 significand')
        numerator //= 2
        exponent += 1
    high, low = divmod(numerator, 2**26)
    expression = f'({high}*67108864+{low})'
    if exponent > 0:
        expression = f'({expression}*2^{exponent})'
    elif exponent < -1022:
        expression = f'(({expression}/2^1022)/2^{-exponent-1022})'
    elif exponent < 0:
        expression = f'({expression}/2^{-exponent})'
    return '-'+expression if value < 0 else expression


def scalar_source(value):
    if isinstance(value, bool):
        return str(value).upper()
    if isinstance(value, str) and re.fullmatch(r'0x[0-9a-fA-F]{16}', value):
        return exact_expression(number(value))
    if isinstance(value, str):
        return '"'+value.replace('"', '""')+'"'
    raise ValueError('scalar seed required')


def choose_functions(problem, limit):
    candidates = [r for r in problem['functions'] if r['function'] not in PARKED
                  and re.fullmatch(r'[A-Z][A-Z0-9_.]*', r['function'])
                  and any(v for source in ('corpus', 'bank') for v in r[source].values())]
    candidates.sort(key=lambda r: (-r['corpus']['larger_numeric'], -r['corpus']['last_bit'],
                                  -r['bank']['larger_numeric'], -r['bank']['last_bit'], r['function']))
    return [r['function'] for r in candidates[:limit]]


def build_questions(seeds, per_function=50, round_number=23, excluded_formulas=()):
    """One-coordinate mutations; each call retains its exact source seed."""
    questions, counts, seen = [], Counter(), set()
    excluded_formulas = set(excluded_formulas)
    seed_counts = Counter(seed['function'] for seed in seeds)
    limits = {}
    def emit(expression, seed, mutation, **extra):
        row = len(questions)+1
        questions.append(dict(id=f'r{round_number}-fuzz-{row:05d}', row=row, column=4, formula=expression,
                              seed_id=seed['id'], function=seed['function'], mutation=mutation, **extra))
        return f'D{row}'
    def call(seed, arguments, mutation, **extra):
        body = seed['function']+'('+','.join(arguments)+')'
        if seed['result_index']:
            body = f"INDEX({body},{seed['result_index'][0]},{seed['result_index'][1]})"
        key = (seed['function'], body)
        if body in excluded_formulas or key in seen or counts[seed['function']] >= limits[seed['function']]:
            return
        seen.add(key)
        counts[seed['function']] += 1
        emit(body, seed, mutation, role='question', **extra)
    for seed in seeds:
        function = seed['function']
        allowance = (per_function-counts[function])//seed_counts[function]
        seed_counts[function] -= 1
        limits[function] = counts[function]+allowance
        original = seed['argument_sources']
        call(seed, original, 'original_control')
        if round_number > 23 and function in ('LINEST', 'LOGEST') and len(seed['args']) >= 2:
            # Paired data mutations preserve the response/predictor relationship.
            # They ask whether rank depends on count, order or exact column scale.
            # Every changed array element is echoed; no invariance is asserted.
            y, x = seed['args'][:2]
            if isinstance(y, list) and y and all(isinstance(r, list) for r in y):
                y = y[0] if len(y) == 1 else [r[0] for r in y] if all(len(r) == 1 for r in y) else []
            else:
                y = []
            if isinstance(x, list) and x and all(isinstance(r, list) for r in x):
                if len(x) == 1 and len(x[0]) == len(y):
                    x = [[v] for v in x[0]]
            else:
                x = []
            numeric = lambda v: isinstance(v, str) and re.fullmatch(r'0x[0-9a-fA-F]{16}', v) and math.isfinite(number(v))
            if y and len(x) == len(y) and all(numeric(v) for v in y) and all(numeric(v) for row in x for v in row):
                def paired_array(rows, index):
                    text = '{'+';'.join(','.join(format(number(v), '.17g') for v in row) for row in rows)+'}'
                    echoes = [emit(f'INDEX({text},{r+1},{c+1})', seed, 'paired_array_operand',
                                   role='operand_echo', argument_index=index, intended_bits=v)
                              for r, row in enumerate(rows) for c, v in enumerate(row)]
                    return text, echoes
                variants = [('reverse_observations', list(reversed(y)), list(reversed(x)))]
                for copies in (2, 4):
                    if len(y)*copies <= 64:
                        variants.append((f'repeat_observations_{copies}', y*copies, x*copies))
                for exponent in (-20, 20):
                    try:
                        transformed = [[bits(math.ldexp(number(v), exponent if c == 0 else 0))
                                        for c, v in enumerate(row)] for row in x]
                    except OverflowError:
                        continue
                    variants.append((f'first_predictor_scale_2^{exponent}', y, transformed))
                for tag, response, predictors in variants:
                    if counts[function] >= limits[function]:
                        break
                    changed = original.copy()
                    changed[0], y_echo = paired_array([[v] for v in response], 0)
                    changed[1], x_echo = paired_array(predictors, 1)
                    call(seed, changed, tag, operand_echoes=y_echo+x_echo,
                         transport='Paired decimal arrays; interpret only the captured operand echoes.')
        for index, value in enumerate(seed['args']):
            if counts[seed['function']] >= limits[function]:
                break
            # Re-use the real formula's array shape and values; no synthetic answer.
            flat = value if isinstance(value, list) else [[value]]
            if not all(isinstance(row, list) for row in flat):
                flat = [flat]
            scalars = [v for row in flat for v in row]
            if not scalars:
                continue
            numeric = next((number(v) for v in scalars if isinstance(v, str)
                            and re.fullmatch(r'0x[0-9a-fA-F]{16}', v)), None)
            replacements = []
            if numeric is not None:
                neighbors = [math.nextafter(numeric, -math.inf), math.nextafter(numeric, math.inf)]
                boundaries = [0.0, -1.0, 1.0, float.fromhex('0x1p-1022')]
                for tag, values in (('neighbor', neighbors), ('boundary', boundaries)):
                    for new in values:
                        if math.isfinite(new):
                            expression = exact_expression(new)
                            echo = emit(expression, seed, tag+'_operand', role='operand_echo',
                                        argument_index=index, intended_bits=bits(new))
                            replacements.append((expression, tag, dict(intended_bits=bits(new), operand_echo=echo)))
                replacements.append(('"'+format(numeric, '.17g')+'"', 'text_number', {}))
            replacements += [(v, tag, {}) for v, tag in
                             [('TRUE', 'boolean_true'), ('FALSE', 'boolean_false'), ('', 'omitted'),
                              ('E1', 'blank_reference'), ('E1:E2', 'blank_range'), ('""', 'empty_text'),
                              ('NA()', 'error_na'), ('1/0', 'error_div0')]]
            # A vertical reference and array use exactly the same original cells.
            # A horizontal/matrix seed is explicitly reshaped for this mutation.
            if len(scalars) <= 12:
                refs = [emit(scalar_source(v), seed, 'range_operand', role='operand_echo',
                             argument_index=index, source_value=v) for v in scalars]
                replacements.append((refs[0]+':'+refs[-1], 'range_values', dict(operand_echoes=refs)))
                literals = []
                for v in scalars:
                    if isinstance(v, str) and re.fullmatch(r'0x[0-9a-fA-F]{16}', v):
                        literals.append(format(number(v), '.17g'))
                    else:
                        literals.append(scalar_source(v))
                array = '{'+';'.join(literals)+'}'
                array_echoes = [emit(f'INDEX({array},{i+1},1)', seed, 'array_literal_operand', role='operand_echo',
                                    argument_index=index, intended_value=value)
                                for i, value in enumerate(scalars)]
                replacements.append((array, 'array_values', dict(operand_echoes=array_echoes,
                    transport='decimal array literals; check Excel echoes before interpreting neighbor/shape behavior')))
            for replacement, tag, details in replacements:
                changed = original.copy()
                changed[index] = replacement
                call(seed, changed, tag, argument_index=index, **details)
    return questions, counts


def generate(corpus, problems, executable, powershell, output, functions=12, per_function=50,
             round_number=23, seeds_per_function=2, extend_from=None):
    corpus, problems, executable, powershell, output = map(lambda p: Path(p).resolve(),
                                                         (corpus, problems, executable, powershell, output))
    if not isinstance(round_number, int) or round_number < 23 or not 1 <= seeds_per_function <= 12:
        raise ValueError('bounded round and seed count required')
    case = f'dcztvv_r{round_number}'
    checkout = Path(__file__).resolve().parents[2]/'DnaCalc'
    generator = checkout/'crates/calc-func/tests/fixtures/make-catalogue-fixtures.ps1'
    inputs = (executable, problems, generator, Path(__file__).resolve())
    input_hashes = {str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    predecessor, excluded_formulas = None, set()
    if extend_from is not None:
        prior = Path(extend_from).resolve()
        prior_manifest = json.loads((prior/'manifest.json').read_text())
        if 'workbook-questions.json' not in prior_manifest['files_sha256']:
            raise ValueError('predecessor question document is not hash-bound')
        for relative, expected in prior_manifest['files_sha256'].items():
            file = (prior/relative).resolve()
            if not file.is_relative_to(prior) or hashlib.sha256(file.read_bytes()).hexdigest() != expected:
                raise ValueError('predecessor packet changed')
        previous = json.loads((prior/'workbook-questions.json').read_text())['questions']
        excluded_formulas = {q['formula'] for q in previous if q['role'] == 'question'}
        predecessor = dict(manifest_sha256=hashlib.sha256((prior/'manifest.json').read_bytes()).hexdigest(),
                           selected_calls=len(excluded_formulas), repeated_question_formulas=0)
    output.mkdir(parents=True, exist_ok=False)
    steering = output/'problem-functions.json'
    steering.write_bytes(problems.read_bytes())
    problem = json.loads(steering.read_text())
    if problem.get('complete') is not True:
        raise ValueError('complete maintained problem list required')
    wanted = choose_functions(problem, functions)
    frozen = output/'corpus-seeds'
    frozen.mkdir()
    # Keep the original source bindings and exclusions, plus selected records.
    # The full candidate pool is reproducible from the corpus and need not be
    # copied into each question packet.
    with tempfile.TemporaryDirectory(prefix='dnacalc-fuzz-seeds-') as scratch:
        pool = Path(scratch)/'pool'
        freeze(corpus, pool, ','.join(wanted))
        raw = (pool/'cases.jsonl').read_bytes()
        rows = [json.loads(line) for line in raw.splitlines()]
        (frozen/'summary.json').write_bytes((pool/'summary.json').read_bytes())
        excluded = (pool/'excluded.jsonl').read_bytes()
        (frozen/'excluded.jsonl.gz').write_bytes(gzip.compress(excluded, mtime=0))
        pool_receipt = dict(cases_sha256=hashlib.sha256(raw).hexdigest(), cases=len(rows),
                            exclusions_sha256=hashlib.sha256(excluded).hexdigest())
        (frozen/'pool-receipt.json').write_text(json.dumps(pool_receipt, indent=2)+'\n')
    seeds, shape_counts = [], Counter()
    for function in wanted:
        candidates = [r for r in rows if r['function'] == function]
        candidates.sort(key=lambda r: (r['verdict'] == 'matched', r['id']))
        shapes = set()
        for row in candidates:
            shape = (tuple(type(v).__name__ for v in row['args']), row['verdict'])
            if round_number > 23:
                shape += (tuple(tuple(map(len, v)) if isinstance(v, list) and all(isinstance(r, list) for r in v)
                                else len(v) if isinstance(v, list) else v if isinstance(v, bool) else None
                                for v in row['args']), tuple(row['result_index'] or ()))
            if shape in shapes:
                continue
            shapes.add(shape)
            if 'date1904=false' not in row['original_evaluation'].get('context', ''):
                raise ValueError('first fuzzer slice requires the recorded 1900 workbook context')
            row['date_system_1900'] = True
            seeds.append(row)
            shape_counts[function] += 1
            if shape_counts[function] == seeds_per_function:
                break
    if not seeds:
        raise ValueError('no supported real corpus seed')
    seed_path = output/'seeds.jsonl'
    seed_path.write_text(''.join(json.dumps(r)+'\n' for r in seeds))
    binding = subprocess.check_output([str(executable), 'provenance'], timeout=60)
    (output/'prediction-binding.json').write_bytes(binding)
    # The existing audit checks each seed operand and complete formula through
    # the registered and Workbook parser routes, preserving the original answer.
    subprocess.run([str(executable), 'function-routes', str(seed_path),
                    str(output/'seed-route-audit.jsonl')], check=True, timeout=600)
    questions, counts = build_questions(seeds, per_function, round_number, excluded_formulas)
    sheet = output/'workbook-questions.json'
    sheet.write_text(json.dumps(dict(sheet=case, allow_refusals=True, questions=questions), indent=2)+'\n')
    section = output/f'sheet-{case}.txt'
    section.write_text(f'# R{round_number} corpus-seeded questions and operand echoes. No Excel answers.\n'+
                       '\n'.join(q['formula'] for q in questions)+'\n', encoding='ascii')
    predictions = output/'application-predictions-workbook.jsonl'
    subprocess.run([str(executable), 'sheet-predictions', str(sheet), str(predictions)], check=True, timeout=600)
    calculated = [json.loads(line) for line in predictions.read_text().splitlines()]
    if len(calculated) != len(questions) or any(a['id'] != b['id'] or a['formula'] != b['formula']
                                              for a, b in zip(questions, calculated)):
        raise ValueError('question predictions changed inputs or omitted formulas')
    subprocess.run([str(powershell), '-NoProfile', '-File', str(generator), '-OutDir', str(output/'workbooks'),
                    '-QuestionSectionFile', str(section), '-QuestionCase', case], check=True, timeout=300)
    package = output/f'workbooks/{case}/oracle.xlsx'
    import zipfile
    import xml.etree.ElementTree as ET
    ns = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(package) as archive:
        cells = ET.fromstring(archive.read('xl/worksheets/sheet1.xml')).findall('.//s:c', ns)
        if len(cells) != len(questions) or any(c.find('s:v', ns) is not None for c in cells):
            raise ValueError('new question workbook has unexpected cells or cached values')
    if any(hashlib.sha256(p.read_bytes()).hexdigest() != input_hashes[str(p)] for p in inputs):
        raise ValueError('packet producer or steering input changed during generation')
    manifest = dict(format='dnacalc-corpus-fuzz-questions/1', new_excel_answers=0, round=round_number,
        predecessor=predecessor, seeds_per_function=seeds_per_function,
        stable_producer_hashes=input_hashes,
        seeds=len(seeds), seeded_functions=dict(shape_counts), calls_by_function=dict(counts),
        selected_calls=sum(counts.values()), formula_cells=len(questions),
        prediction_provenance=dict(Counter(r['prediction_provenance'] for r in calculated)),
        mutations=dict(Counter(q['mutation'] for q in questions if q['role'] == 'question')),
        sources=dict(problem_list=steering.name, problem_list_sha256=hashlib.sha256(steering.read_bytes()).hexdigest(),
            steering_report_original_path=str(problems),
            corpus_summary_sha256=hashlib.sha256((corpus/'summary.json').read_bytes()).hexdigest(),
            generator_sha256=hashlib.sha256(generator.read_bytes()).hexdigest(),
            prediction_executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
            fuzzer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        limits='Bounded direct scalar/literal-array calls accepted by the existing corpus freezer. '
            'R24 adds array-shape, boolean-flag and selected-result seed classes; predecessor question formulas are omitted. '
            'Reference/nested-expression seeds are excluded and counted by that freezer. '
            'Generic mutations change one argument; paired regression mutations change response and predictor arrays together. Arrays may be explicitly reshaped. '
            'Finite binary64 neighbors use exact small-integer expressions; underflow, signed zero and array decimal transport require echoes. '
            'Application predictions are not Excel evidence. Original controls are intentional rechecks. '
            'No corpus download, captured workbook edit, fitted numerical constant or answer-bearing request.')
    (output/'README.txt').write_text(f'''R{round_number} corpus-seeded discovery questions; no Excel answers.
The manifest lists seeds, mutation classes, counts and immutable hashes.
Only the new workbooks/{case}/oracle.xlsx is a capture input. Never append
these questions to any captured workbook. From the updated Windows checkout:
  cargo xtask capture --fixture <packet>/workbooks/{case}/oracle.xlsx --out <results>/r{round_number}.capture.json
  cargo xtask compare --fixture <packet>/workbooks/{case}/oracle.xlsx --capture <results>/r{round_number}.capture.json --out <results>/comparison
Return the unchanged packet, complete canonical capture, saved package and all
worker/settings/partial-attempt artefacts under ~/codex-runs/excel-results-r{round_number}.
Use the existing capture lineage. Operand echoes establish what Excel consumed;
neighbor intent is not an asserted observation. No expected values are supplied.
Predictions are application calculations only. Missing, unread and refused calls
remain failures; report all formula cells and selected calls separately.
''')
    manifest['files_sha256'] = {str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(output.rglob('*')) if p.is_file()}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k not in ('sources', 'files_sha256')}))
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('corpus', 'problems', 'executable', 'powershell', 'output'):
        parser.add_argument(name)
    parser.add_argument('--functions', type=int, default=12)
    parser.add_argument('--per-function', type=int, default=50)
    parser.add_argument('--round-number', type=int, default=23)
    parser.add_argument('--seeds-per-function', type=int, default=2)
    parser.add_argument('--extend-from')
    options = parser.parse_args()
    generate(**vars(options))
