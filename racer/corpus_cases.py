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
    for workbook in summary['workbooks']:
        label = workbook['name'].replace('/', '_').replace('\\', '_')
        directory = corpus/'workbooks'/label
        captures = sorted(directory.glob('*.capture.json'))
        if captures:
            runs.extend((capture, directory/capture.name.removesuffix('.capture.json'), workbook['name'])
                        for capture in captures)
        elif workbook['capture'].startswith('committed capture '):
            runs.append((Path(workbook['capture'].removeprefix('committed capture ')), directory, workbook['name']))
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
            if verdict.get('class') == 'capture':
                reason = 'capture_failure'
            elif verdict['verdict'] != 'mismatched':
                reason = verdict['verdict']
            elif isinstance(x, dict) and isinstance(y, dict):
                if set(x) == set(y) == {'number'}:
                    import math
                    if math.isfinite(number(x['number'])) and math.isfinite(number(y['number'])):
                        distance = abs(ordered(x['number'])-ordered(y['number']))
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
                original_verdict=verdict, expected=e, actual=a,
                capture=str(capture_path), evaluation=str(paths[1])))
    require_count = compared['after_fail'] if compared is not None else summary['totals']['counts']['fail']
    if len(rows) != require_count:
        raise ValueError(f'nonpass accounting differs: {len(rows)} != {require_count}')
    (output/'cells.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    result = dict(complete=True, failed=len(rows), counts=dict(counts),
        groups=[dict(category=c, function=f, cause=r, cells=n)
                for (c,f,r), n in groups.most_common()], inputs_sha256=inputs,
        cells_sha256=hashlib.sha256((output/'cells.jsonl').read_bytes()).hexdigest(),
        rule='Last-bits-only: finite numbers <=16 ULP apart. Larger errors are a triage group, not a mathematical root-cause proof. '
             'Capture failures, refusals and missing values remain failures, explicitly distinguished from proven product bugs.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('inputs_sha256','groups')}))
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


if __name__ == '__main__':
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
