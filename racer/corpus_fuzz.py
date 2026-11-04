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
                  and re.fullmatch(r'[A-Z][A-Z0-9_.]*', r['function'])]
    candidates.sort(key=lambda r: (-r['corpus']['larger_numeric'], -r['corpus']['last_bit'],
                                  -r['bank']['larger_numeric'], -r['bank']['last_bit'], r['function']))
    return [r['function'] for r in candidates[:limit]]


def build_questions(seeds, per_function=50):
    """One-coordinate mutations; each call retains its exact source seed."""
    questions, counts, seen = [], Counter(), set()
    seed_counts = Counter(seed['function'] for seed in seeds)
    limits = {}
    def emit(expression, seed, mutation, **extra):
        row = len(questions)+1
        questions.append(dict(id=f'r23-fuzz-{row:05d}', row=row, column=4, formula=expression,
                              seed_id=seed['id'], function=seed['function'], mutation=mutation, **extra))
        return f'D{row}'
    def call(seed, arguments, mutation, **extra):
        body = seed['function']+'('+','.join(arguments)+')'
        if seed['result_index']:
            body = f"INDEX({body},{seed['result_index'][0]},{seed['result_index'][1]})"
        key = (seed['function'], body)
        if key in seen or counts[seed['function']] >= limits[seed['function']]:
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


def generate(corpus, problems, executable, powershell, output, functions=12, per_function=50):
    corpus, problems, executable, powershell, output = map(lambda p: Path(p).resolve(),
                                                         (corpus, problems, executable, powershell, output))
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
            if shape in shapes:
                continue
            shapes.add(shape)
            if 'date1904=false' not in row['original_evaluation'].get('context', ''):
                raise ValueError('first fuzzer slice requires the recorded 1900 workbook context')
            row['date_system_1900'] = True
            seeds.append(row)
            shape_counts[function] += 1
            if shape_counts[function] == 2:
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
    questions, counts = build_questions(seeds, per_function)
    sheet = output/'workbook-questions.json'
    sheet.write_text(json.dumps(dict(sheet='dcztvv_r23', allow_refusals=True, questions=questions), indent=2)+'\n')
    section = output/'sheet-dcztvv_r23.txt'
    section.write_text('# R23 corpus-seeded questions and operand echoes. No Excel answers.\n'+
                       '\n'.join(q['formula'] for q in questions)+'\n', encoding='ascii')
    predictions = output/'application-predictions-workbook.jsonl'
    subprocess.run([str(executable), 'sheet-predictions', str(sheet), str(predictions)], check=True, timeout=600)
    calculated = [json.loads(line) for line in predictions.read_text().splitlines()]
    if len(calculated) != len(questions) or any(a['id'] != b['id'] or a['formula'] != b['formula']
                                              for a, b in zip(questions, calculated)):
        raise ValueError('question predictions changed inputs or omitted formulas')
    checkout = Path(__file__).resolve().parents[2]/'DnaCalc'
    generator = checkout/'crates/calc-func/tests/fixtures/make-catalogue-fixtures.ps1'
    subprocess.run([str(powershell), '-NoProfile', '-File', str(generator), '-OutDir', str(output/'workbooks'),
                    '-QuestionSectionFile', str(section), '-QuestionCase', 'dcztvv_r23'], check=True, timeout=300)
    package = output/'workbooks/dcztvv_r23/oracle.xlsx'
    import zipfile
    import xml.etree.ElementTree as ET
    ns = {'s':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(package) as archive:
        cells = ET.fromstring(archive.read('xl/worksheets/sheet1.xml')).findall('.//s:c', ns)
        if len(cells) != len(questions) or any(c.find('s:v', ns) is not None for c in cells):
            raise ValueError('new question workbook has unexpected cells or cached values')
    manifest = dict(format='dnacalc-corpus-fuzz-questions/1', new_excel_answers=0,
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
        limits='First bounded slice: direct scalar/literal-array calls accepted by the existing corpus freezer, '
            'up to two argument-class/verdict seeds per function. Reference/nested-expression seeds are excluded and counted by that freezer. '
            'One argument changes at a time; arrays may be explicitly reshaped to a column. '
            'Finite binary64 neighbors use exact small-integer expressions; underflow, signed zero and array decimal transport require echoes. '
            'Application predictions are not Excel evidence. Original controls are intentional rechecks. '
            'No corpus download, captured workbook edit, fitted numerical constant or answer-bearing request.')
    (output/'README.txt').write_text('''R23 corpus-seeded discovery questions; no Excel answers.
The manifest lists seeds, mutation classes, counts and immutable hashes.
Only the new workbooks/dcztvv_r23/oracle.xlsx is a capture input. Never append
these questions to any captured workbook. From the updated Windows checkout:
  cargo xtask capture --fixture <packet>/workbooks/dcztvv_r23/oracle.xlsx --out <results>/r23.capture.json
  cargo xtask compare --fixture <packet>/workbooks/dcztvv_r23/oracle.xlsx --capture <results>/r23.capture.json --out <results>/comparison
Return the unchanged packet, complete canonical capture, saved package and all
worker/settings/partial-attempt artefacts under ~/codex-runs/excel-results-r23.
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
    options = parser.parse_args()
    generate(**vars(options))
