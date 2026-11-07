"""Exact elimination of ordinary digit-extraction conversion programs.

Sources: IEEE754 and decimal positional notation. Only powers of ten, 0, 1,
and the half-unit rounding point occur. No parameter recovery or exceptions.
Apache POI source was inspected; this is not a complete POI implementation.
Existing Excel observations are discovery, never fresh validation.
"""
from decimal import Decimal
from fractions import Fraction as Q
from pathlib import Path
import ast
import itertools
import json
import sys
import zipfile
import xml.etree.ElementTree as E
from functools import lru_cache
from collections import Counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from checker.exact import Lattice


def arithmetic(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    if isinstance(node, ast.UnaryOp):
        v = arithmetic(node.operand)
        return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.BinOp):
        a, b = arithmetic(node.left), arithmetic(node.right)
        if isinstance(node.op, ast.Add): return a + b
        if isinstance(node.op, ast.Sub): return a - b
        if isinstance(node.op, ast.Mult): return a * b
        if isinstance(node.op, ast.Div): return a / b
    raise ValueError(ast.dump(node))


def observations(main):
    """Inspect existing canonical capture values and their original workbook inputs."""
    ns = {'x': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    rows = []
    for case, selected_sheet in [('dc_q9rb', 'number_text_for_concatenation'), ('catalogue', 'dch9u_functions')]:
        capture = main / 'crates/calc-func/tests/captures' / case
        if case == 'catalogue': capture = capture / 'sections/dch9u_functions.capture.json'
        else: capture = capture / 'oracle.capture.json'
        doc = json.loads(capture.read_text())
        workbook = main / 'crates/calc-func/tests/fixtures' / case / 'oracle.xlsx'
        with zipfile.ZipFile(workbook) as z:
            sheets = E.fromstring(z.read('xl/workbook.xml')).find('x:sheets', ns)
            index = next(i for i, s in enumerate(sheets, 1) if s.attrib['name'] == selected_sheet)
            xml = E.fromstring(z.read(f'xl/worksheets/sheet{index}.xml'))
            cells = {c.attrib['r']: c for c in xml.findall('.//x:c', ns)}
        formulas = {(r['row'], r['column']): r for r in doc['formulas'] if r['sheet'] == selected_sheet}
        for obs in doc['observations']:
            key = obs['row'], obs['column']
            if obs['sheet'] != selected_sheet or key not in formulas or 'text' not in obs.get('excel_raw', {}): continue
            formula = formulas[key]['formula'].lstrip('=')
            if not formula.endswith('&""'): continue
            source = formula[:-3]
            try:
                if source.startswith('A') and source[1:].isdigit():
                    x = float(cells[source].find('x:v', ns).text)
                else:
                    x = arithmetic(ast.parse(source, mode='eval').body)
            except (ValueError, AttributeError, SyntaxError): continue
            text = bytes.fromhex(obs['excel_raw']['text']).decode('utf-16-be')
            if x == 0 or abs(x) < sys.float_info.min: continue
            rows.append(dict(case=case, sheet=selected_sheet, row=obs['row'], column=obs['column'],
                formula=formula, input=x.hex(), expected_text=text,
                exact_input=Q(abs(x)), expected_decimal=abs(Q(Decimal(text)))))
    return rows


def ten(e):
    return Q(10**e) if e >= 0 else Q(1, 10**-e)


@lru_cache(maxsize=None)
def power64(exponent):
    return Lattice(64, -20000, 20000).round(ten(exponent))


def radix16_scale(value, exponent):
    """Fixed low-nibble-first powers-of-ten scaling, all operations PC64."""
    rnd = Lattice(64, -20000, 20000).round
    magnitude = abs(exponent)
    for shift in range(0, magnitude.bit_length(), 4):
        part = ((magnitude >> shift) & 15) << shift
        if part: value = rnd(value*power64(-part if exponent<0 else part))
    return value


def rounding_radix_controls(checkout, output, cache, offline=False):
    """Paired calls that change coefficient length while keeping decimal value.

    Every selected native mining result is independently recomputed using exact
    rational PC64 operations before Excel is asked. No answer selects a case.
    """
    import numpy as np
    import struct
    try:
        from .fingerprint import load_service, CachedService
    except ImportError:
        from fingerprint import load_service, CachedService
    if np.finfo(np.longdouble).nmant != 63: raise ValueError('PC64 native miner required')
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    output = Path(output)
    output.mkdir(exist_ok=False)
    encode = lambda x: struct.pack('>d', float(x)).hex()
    def native_scale(values, exponent):
        values = values.astype(np.longdouble)
        magnitude = abs(exponent)
        for shift in range(0, magnitude.bit_length(), 4):
            part = ((magnitude >> shift)&15)<<shift
            if part: values = values*np.longdouble('1e'+str(-part if exponent<0 else part))
        return values.astype(np.float64)
    cases, seen = [], set()
    rng = np.random.default_rng(340016)
    for digits in range(1, 15):
        accepted = Counter()
        for exponent in range(-307, 295, 7):
            side = exponent < 0
            if accepted[side] >= 12: continue
            coefficients = rng.integers(10**(digits-1), 10**digits, size=2000, dtype=np.int64)
            coefficients = coefficients[coefficients%10 != 0]
            padding = 15-digits
            normalized = native_scale(coefficients, exponent)
            fixed = native_scale(coefficients*10**padding, exponent-padding)
            selected = np.flatnonzero((normalized != fixed) & np.isfinite(normalized)
                                     & (np.abs(normalized) >= sys.float_info.min))
            for index in selected:
                c = int(coefficients[index])
                identity = c, exponent
                if identity in seen: continue
                seen.add(identity)
                exact = c*ten(exponent)
                short = radix16_scale(Q(c), exponent)
                full = radix16_scale(Q(c*10**padding), exponent-padding)
                if encode(short) != encode(normalized[index]) or encode(full) != encode(fixed[index]):
                    raise ValueError('native/rational scaling disagreement')
                for sign in (1, -1):
                    cases.append(dict(coefficient=sign*c, exponent=exponent, digits=digits,
                        input_bits=encode(sign*exact), normalized_bits=encode(sign*short),
                        fixed15_bits=encode(sign*full)))
                accepted[side] += 1
                if accepted[side] >= 12: break
    if not cases: raise ValueError('empty representation discrimination cohort')
    (output/'questions.json').write_text(json.dumps(cases, separators=(',', ':'))+'\n')
    counts, replies = {}, {}
    arguments = [[client.number(c['input_bits']), client.number(digits)]
        for c in cases for digits in (32000., float(-c['exponent']))]
    for function in ('ROUNDUP', 'ROUNDDOWN', 'TRUNC', 'ROUND'):
        answers = service.ask_rows(function, arguments, chunk_size=2000)
        replies[function] = answers
        counts[function] = dict(rows=len(arguments),
            whole_fixed15=sum(a.get('bits') == c['fixed15_bits'] for c,a in zip(cases,answers[::2])),
            whole_normalized=sum(a.get('bits') == c['normalized_bits'] for c,a in zip(cases,answers[::2])),
            position_fixed15=sum(a.get('bits') == c['fixed15_bits'] for c,a in zip(cases,answers[1::2])),
            position_normalized=sum(a.get('bits') == c['normalized_bits'] for c,a in zip(cases,answers[1::2])))
    (output/'answers.json').write_text(json.dumps(replies, separators=(',', ':'))+'\n')
    result = dict(complete=True, cases=len(cases), counts=counts, receipts=service.receipts,
        source='Decimal positional arithmetic; fixed radix16 PC64 scaling and independently selected native/rational disagreements.',
        limits='Paired whole-reading and decimal-position publications only. No initial-reading algorithm or broader scaling table is installed.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(counts))


def rounding_publications(checkout, output, cache, offline=False):
    """Preselect exact double-rounding discriminants before querying Excel."""
    import hashlib
    import struct
    try:
        from .fingerprint import load_service, CachedService
    except ImportError:
        from fingerprint import load_service, CachedService
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    output = Path(output)
    output.mkdir(exist_ok=False)
    lattice = Lattice(64, -20000, 20000)
    bits = lambda x: struct.pack('>d', float(x)).hex()
    cases = []
    # An input-independent integer sequence and complete exponent range. Neither
    # the Excel answers nor the bank's match counts select a coefficient.
    seed = 340011
    for exponent in range(-321, 295):
        power = ten(exponent)
        for _ in range(1000):
            seed = (6364136223846793005*seed+1442695040888963407) % 2**64
            coefficient = 100000000000000 + seed % 900000000000000
            exact = coefficient*power
            if not Q(sys.float_info.min) <= exact <= Q(sys.float_info.max):
                continue
            direct, stored = bits(exact), bits(lattice.round(exact))
            if direct != stored:
                for sign in (1, -1):
                    cases.append(dict(coefficient=coefficient*sign, exponent=exponent,
                        input_bits=bits(sign*exact), strict_bits=bits(sign*exact),
                        stored_bits=bits(sign*lattice.round(exact))))
    if not cases:
        raise ValueError('no arithmetic discriminants')
    (output/'questions.json').write_text(json.dumps(cases, separators=(',', ':'))+'\n')
    arguments = [[client.number(c['input_bits']), client.number(32000.)] for c in cases]
    counts, all_answers = {}, {}
    for function in ('ROUNDUP', 'ROUNDDOWN', 'TRUNC', 'ROUND'):
        answers = service.ask_rows(function, arguments)
        counts[function] = dict(rows=len(cases), stored=sum(a.get('bits') == c['stored_bits']
            for a, c in zip(answers, cases)), strict=sum(a.get('bits') == c['strict_bits']
            for a, c in zip(answers, cases)))
        all_answers[function] = answers
    (output/'answers.json').write_text(json.dumps(all_answers, separators=(',', ':'))+'\n')
    result = dict(complete=True, counts=counts, receipts=service.receipts,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        source='IEEE754 exact decimal rational, RN64 followed by RN53; independently selected arithmetic disagreements.',
        limits='Whole fifteen-digit decimal publication only. Initial digit extraction and final digit rounding remain separate unresolved stages.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(counts))


def rounding_scale128(checkout, output, cache, offline=False):
    """Discriminate one negative power-table operand with no other table factor.

    The declared program scales x by a resident10^-128, then by exact10, and
    rounds the resulting15-digit coefficient half up. Binary exponent selects
    the upper binade's decimal exponent. No value is installed by this study.
    """
    import numpy as np
    import struct
    from checker.exact import EXTENDED64
    try:
        from .fingerprint import load_service, CachedService
    except ImportError:
        from fingerprint import load_service, CachedService
    if np.finfo(np.longdouble).nmant != 63:
        raise ValueError('control mining requires the native64-bit significand')
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    output = Path(output)
    output.mkdir(exist_ok=False)
    p = EXTENDED64.round(ten(-128))
    exponent = p.numerator.bit_length()-p.denominator.bit_length()
    if p < Q(2)**exponent: exponent -= 1
    unit = Q(2)**(exponent-63)
    native_p = np.longdouble('1e-128')
    if Q(str(native_p)) == 0:
        raise ValueError('native seed underflow')
    cases = []
    seen = set()
    for batch in range(8):
        coefficients = np.arange(620000000000000+batch*100000,
                                 620000000100000+batch*100000, dtype=np.int64)
        xs = ((coefficients.astype(np.longdouble)+np.longdouble('.5'))/10/native_p).astype(np.float64)
        predictions = []
        for offset in range(-1, 4):
            candidate_p = native_p
            for _ in range(abs(offset)):
                candidate_p = np.nextafter(candidate_p, np.longdouble('-inf' if offset<0 else 'inf'))
            z = xs.astype(np.longdouble)*candidate_p
            z = z*10
            predictions.append(np.floor(z+np.longdouble('.5')).astype(np.int64))
        selected = np.flatnonzero(np.any(np.array(predictions) != predictions[0], axis=0))
        for index in selected:
            x = float(xs[index])
            encoded = struct.pack('>d', x).hex()
            if encoded in seen: continue
            seen.add(encoded)
            # Recheck every selected native mining result with exact rationals.
            b = Q(x).numerator.bit_length()-Q(x).denominator.bit_length()
            if Q(x) < Q(2)**b: b -= 1
            if decade(Q(2)**(b+1)) != 142 or decade(Q(x)) != 141:
                raise ValueError('control does not isolate the declared exponent path')
            exact = [int(EXTENDED64.round(EXTENDED64.round(Q(x)*(p+offset*unit))*10)+Q(1, 2))
                     for offset in range(-1, 4)]
            if exact != [int(v[index]) for v in predictions]:
                raise ValueError('native miner differs from rational stage replay')
            for sign in (1, -1):
                cases.append(dict(input_bits=struct.pack('>d', sign*x).hex(),
                                  sign=sign, predictions=exact))
            if len(cases) >= 256: break
        if len(cases) >= 256: break
    if len(cases) != 256: raise ValueError('not enough independent table controls')
    (output/'questions.json').write_text(json.dumps(cases, separators=(',', ':'))+'\n')
    counts, replies = {}, {}
    arguments = [[client.number(c['input_bits']), client.number(32000.)] for c in cases]
    for function in ('ROUNDUP', 'ROUNDDOWN', 'TRUNC'):
        replies[function] = service.ask_rows(function, arguments)
        matches = [0]*5
        for c, answer in zip(cases, replies[function]):
            if answer['type'] != 'number': raise ValueError('unread numeric table control')
            got = abs(Q(struct.unpack('>d', bytes.fromhex(answer['bits']))[0]))
            coefficient = round(got/ten(127))
            for k, expected in enumerate(c['predictions']): matches[k] += coefficient == expected
        counts[function] = matches
    (output/'answers.json').write_text(json.dumps(replies, separators=(',', ':'))+'\n')
    result = dict(complete=True, cases=len(cases), offsets=list(range(-1, 4)),
        matching_initial_coefficients=counts, receipts=service.receipts,
        limits='Conditional constant-operand discrimination in the stated two-multiply initial-reading program. No count-selected value or complete conversion kernel is installed.')
    (output/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(counts))


def decade(x):
    e = len(str(x.numerator)) - len(str(x.denominator))
    return e if x >= ten(e) else e - 1


def evaluate(x, precision, scale, power, normalization, scale_round='nearest', digit_round='nearest', ties='away', extraction='digits'):
    lattice = Lattice(precision, -20000, 20000)
    rnd = lattice.round
    e = decade(x)
    shift = e + normalization if extraction == 'digits' else e - 14
    if power == 'exact': p = ten(-shift if scale == 'multiply' else shift)
    elif power == 'nearest': p = rnd(ten(-shift if scale == 'multiply' else shift))
    else:
        p0 = ten(-shift if scale == 'multiply' else shift)
        p = lattice.value(lattice.floor_coordinate(p0))
    rounded = lambda v, mode: rnd(v) if mode == 'nearest' else lattice.value(lattice.floor_coordinate(v))
    y = rounded(x*p if scale == 'multiply' else x/p, scale_round)
    if extraction == 'scaled_integer':
        coefficient = int(y // 1)
        fraction = y-coefficient
        up = fraction > Q(1,2) or (fraction == Q(1,2) and (ties=='away' or (ties=='even' and coefficient%2)))
        return (coefficient+up)*ten(e-14)
    if normalization: y = rounded(y*10, digit_round)
    digits = []
    for _ in range(15):
        d = y // 1
        digits.append(int(d))
        y = rounded(rnd(y-d)*10, digit_round)
    coefficient = 0
    for d in digits: coefficient = coefficient*10+d
    coefficient += y > 5 or (y == 5 and (ties=='away' or (ties=='even' and coefficient%2)))
    return coefficient*ten(e-14)


def run(main, out, stage='nearest'):
    # The existing scientific renderer rounds again to14 digits for |exponent|>=99.
    # Keep those observations outside this first-rounding necessary constraint.
    all_rows = observations(Path(main))
    rows = [r for r in all_rows if abs(decade(r['exact_input'])) < 99]
    programs = []
    rounding = ['nearest','truncated'] if stage in ('digit-staging','tie-rules') else ['nearest']
    ties = ['away','down','even'] if stage == 'tie-rules' else ['away']
    extraction = ['digits','scaled_integer'] if stage == 'tie-rules' else ['digits']
    for precision, scale, power, normalization, scale_round, digit_round, tie, extract in itertools.product([53,64,96], ['multiply','divide'], ['exact','nearest','truncated'], [0,1], rounding, rounding, ties, extraction):
        if extract == 'scaled_integer' and (normalization or digit_round != 'nearest'): continue
        p = dict(precision=precision, scale=scale, power=power, normalization=normalization,scale_round=scale_round,digit_round=digit_round,ties=tie,extraction=extract)
        first = None; matches = 0
        for r in rows:
            got = evaluate(r['exact_input'], **p)
            passed = got == r['expected_decimal']
            matches += passed
            if not passed and first is None:
                first = {k:v for k,v in r.items() if k not in ['exact_input','expected_decimal']}
                first['predicted_decimal'] = str(got)
        programs.append(dict(program=p, matches=matches, failures=len(rows)-matches, first_kill=first))
    result = dict(rows=len(rows), excluded_second_rounding=len(all_rows)-len(rows), candidates=len(programs), survivors=sum(p['failures']==0 for p in programs),
        programs=programs, scope='Decimal value only; notation/sign rules are outside this necessary constraint.',
        exposure='All existing observations are discovery; no fitted parameter or fresh validation.',
        inputs=[{k:v for k,v in r.items() if k not in ['exact_input','expected_decimal']} for r in rows])
    with Path(out).open('x') as f: json.dump(result, f, indent=2); f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['programs','inputs']}))


def chunked_evaluate(x, precision, scale, power, chunk_width, scale_round, digit_round, ties):
    """Decimal positional extraction in fixed chunks, with explicit rounding."""
    lattice = Lattice(precision, -20000, 20000)
    trace = []
    def rounded(value, mode):
        result = lattice.round(value) if mode == 'nearest' else lattice.value(lattice.floor_coordinate(value))
        trace.append(dict(exact=str(value), rounded=str(result), mode=mode))
        return result
    exponent = decade(x)
    shift = exponent + 1
    exact_power = ten(-shift if scale == 'multiply' else shift)
    p = exact_power if power == 'exact' else rounded(exact_power, 'nearest' if power == 'nearest' else 'truncated')
    y = rounded(x*p if scale == 'multiply' else x/p, scale_round)
    coefficient = 0
    remaining = 15
    while remaining:
        take = min(chunk_width, remaining)
        base = 10**take
        y = rounded(y*base, digit_round)
        digits = y // 1
        coefficient = coefficient*base + digits
        y = rounded(y-digits, digit_round)
        remaining -= take
    up = y > Q(1,2) or (y == Q(1,2) and (ties == 'away' or (ties == 'even' and coefficient % 2)))
    return (coefficient+up)*ten(exponent-14), trace


def chunked_run(main, output):
    import hashlib
    rows = observations(Path(main))
    scoped = [row for row in rows if abs(decade(row['exact_input'])) < 99]
    programs = []
    for precision, scale, power, width, scale_round, digit_round, ties in itertools.product(
            (53,64,96), ('multiply','divide'), ('exact','nearest','truncated'), range(1,16),
            ('nearest','truncated'), ('nearest','truncated'), ('away','down','even')):
        program = dict(precision=precision, scale=scale, power=power, chunk_width=width,
                       scale_round=scale_round, digit_round=digit_round, ties=ties)
        first_kill = None
        for order, row in enumerate(scoped):
            got, trace = chunked_evaluate(row['exact_input'], **program)
            if got != row['expected_decimal']:
                first_kill = dict(order=order, observation={k:v for k,v in row.items()
                    if k not in ['exact_input','expected_decimal']}, exact_input=str(row['exact_input']),
                    expected_decimal=str(row['expected_decimal']), prediction=str(got), trace=trace)
                break
        programs.append(dict(id=len(programs), program=program, first_kill=first_kill))
    assert len(programs) == 3240
    result = dict(complete=True, programs=len(programs), survivors=sum(p['first_kill'] is None for p in programs),
        scoped_rows=len(scoped), excluded_second_rounding=len(rows)-len(scoped), profiles=programs,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        arithmetic_sha256=hashlib.sha256((Path(__file__).resolve().parents[1]/'checker/exact.py').read_bytes()).hexdigest(),
        limits='Fixed declared positional extraction family only, original authored operand bits and recorded text. Existing observations are discovery. No fitted switch, selected intermediate, full Excel conversion algorithm or independent validation.')
    with Path(output).open('x') as sink:
        json.dump(result, sink, indent=2)
        sink.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'profiles'}))


def published_poi(main, out, jdk_bin, poi_source, application_records):
    """Execute the complete published converter without translating its arithmetic."""
    import hashlib
    import struct
    import subprocess
    out, jdk_bin, poi_source = map(Path, (out, jdk_bin, poi_source))
    manifest_path = poi_source / 'source-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    package = poi_source / 'org/apache/poi/ss/util'
    for name, identity in manifest['files'].items():
        source = package / name if name.endswith('.java') else poi_source / name
        if hashlib.sha256(source.read_bytes()).hexdigest() != identity['sha256']:
            raise ValueError('Published source changed: ' + name)
    wrapper = Path(__file__).with_name('NumberTextProbe.java')
    classes = out / 'classes'
    out.mkdir()
    classes.mkdir()
    compile_command = [str(jdk_bin/'javac'), '-d', str(classes), str(wrapper),
                       *(str(p) for p in sorted(package.glob('*.java')))]
    compiled = subprocess.run(compile_command, capture_output=True, text=True)  # ubs:ignore — explicitly chosen local compiler; argv only, no shell
    (out/'compile.log').write_text(compiled.stdout + compiled.stderr)
    compiled.check_returncode()
    rows = observations(Path(main))
    bits = [struct.pack('>d', float.fromhex(row['input'])).hex() for row in rows]
    command = [str(jdk_bin/'java'), '-cp', str(classes), 'NumberTextProbe']
    result = subprocess.run(command, input='\n'.join(bits)+'\n', capture_output=True,  # ubs:ignore — explicitly chosen local JVM; argv only, no shell
                            text=True, check=True)
    (out/'output.tsv').write_text(result.stdout)
    (out/'stderr.log').write_text(result.stderr)
    lines = result.stdout.splitlines()
    if len(lines) != len(rows):
        raise ValueError('Published converter row count differs')
    records = []
    for row, tag, line in zip(rows, bits, lines):
        echoed, got = line.split('\t')
        if echoed != tag:
            raise ValueError('Published converter operand echo differs')
        records.append({**{k:v for k,v in row.items() if k not in ['exact_input','expected_decimal']},
                        'bits':tag, 'POI':got, 'matched':got == row['expected_text']})
    key = lambda row: (row['case'], row['sheet'], row['row'], row['column'])
    indexed = {key(row):row for row in records}
    application_path = Path(application_records)
    application = json.loads(application_path.read_text())
    delta = []
    for actual in application:
        original = actual['source']
        predicted = indexed[key(original)]
        if predicted['bits'] != actual['admitted_bits'] or predicted['expected_text'] != original['expected_text']:
            raise ValueError('Application binding changed')
        delta.append(dict(source=original, application=actual['calculated_text'], POI=predicted['POI'],
                          application_match=actual['matches_original_text'], POI_match=predicted['matched']))
    summary = dict(rows=len(records), matches=sum(row['matched'] for row in records),
        failures=sum(not row['matched'] for row in records),
        audited_overlap=dict(rows=len(delta), application_matches=sum(row['application_match'] for row in delta),
            POI_matches=sum(row['POI_match'] for row in delta),
            gains=sum(not row['application_match'] and row['POI_match'] for row in delta),
            losses=sum(row['application_match'] and not row['POI_match'] for row in delta)),
        inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                (Path(__file__), wrapper, manifest_path, application_path)},
        commands=[compile_command, command], records=records, application_delta=delta,
        limits='Unmodified Apache POI sources on original authored operand bits. Existing captures are discovery; effective Excel input readback is unavailable. No tuned rule or product implementation.')
    with (out/'result.json').open('x') as output:
        json.dump(summary, output, indent=2)
        output.write('\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['records','application_delta']}))


if __name__ == '__main__':
    if sys.argv[1] == 'rounding-radix-controls': rounding_radix_controls(*sys.argv[2:])
    elif sys.argv[1] == 'rounding-scale128': rounding_scale128(*sys.argv[2:])
    elif sys.argv[1] == 'rounding-publications': rounding_publications(*sys.argv[2:])
    elif sys.argv[1] == 'published-poi': published_poi(*sys.argv[2:])
    elif sys.argv[1] == 'chunked': chunked_run(*sys.argv[2:])
    else: run(*sys.argv[1:])
