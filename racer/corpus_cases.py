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


def freeze(corpus, output, functions):
    corpus, output = Path(corpus), Path(output)
    output.mkdir()
    wanted = set(functions.split(','))
    counts = Counter()
    inputs = {}
    with (output/'cases.jsonl').open('x') as out, (output/'excluded.jsonl').open('x') as excluded:
        for capture_path in sorted(corpus.glob('workbooks/*/*.capture.json')):
            capture = json.loads(capture_path.read_text())
            evaluation_path = capture_path.parent/capture_path.name.removesuffix('.capture.json')/'evaluations.json'
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
                m = re.fullmatch(r'(?:_xlfn\.)?([A-Za-z.]+)\((.*)\)', formula)
                if not m or m[1].upper() not in wanted:
                    continue
                function = m[1].upper()
                e, a, v = expected[key(f)], actual[key(f)], verdicts[key(f)]
                raw = dict(id=f"corpus:{capture_path.parent.name}:{f['sheet']}:{f['row']}:{f['column']}",
                           function=function, result_index=selected, formula=f['formula'],
                           source=str(capture_path), reading=capture['environment'],
                           original_capture=e, original_evaluation=a,
                           original_formula=f, verdict=v['verdict'])
                counts[function+':'+v['verdict']] += 1
                try:
                    sources = split_arguments(m[2])
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


if __name__ == '__main__':
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
