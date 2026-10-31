"""Check research certificates, independently of search code (stdlib only).

All accepted observations must already be decoded by the existing xtask witness
loader. This reads its frozen records; it is not a new Excel capture/importer.
Failures abort without a success summary. No tolerance, fitting or match ranking.
"""
import argparse
from collections import Counter
from fractions import Fraction as Rat
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import time

from exact import (BINARY64 as B, EXTENDED64 as E, Lattice, Span, binary64,
                   inverse, operation, round_model, two)


def require(condition, detail):
    if not condition:
        raise ValueError(detail)


def records(path):
    with open(path) as source:
        for line in source:
            yield json.loads(line)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as source:
        for chunk in iter(lambda: source.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def source_hashes():
    return {name: digest(Path(__file__).parent/name) for name in ('exact.py', 'check.py')}


# Capture before processing observations and refuse a receipt if source changes.
LOADED_SOURCE_SHA256 = source_hashes()


def validate_row(row):
    """Bind derived research leaves to the original typed witness arguments."""
    raw = row['raw']
    args = raw['args']
    require(raw['function'] == 'PMT' and raw['result_index'] is None,
            'unsupported surface or selected result element')
    require(3 <= len(args) <= 5 and len(row['leaves']) == 5, 'wrong PMT input arity')
    require(all(isinstance(a, str) and a.startswith('0x') and len(a) == 18 for a in args),
            'only original numeric argument encodings are in this scope')
    bits = [int(a, 16) for a in args]
    values = list(map(binary64, bits))
    require(all(v == 0 or abs(v) >= two(-1022) for v in values), 'subnormal input admission')
    fv_bits = bits[3] if len(bits) >= 4 else 0
    timing = int(values[4]) != 0 if len(values) >= 5 else False  # truncate toward zero
    require(row['leaves'][:4] == [bits[0], bits[2], fv_bits, 0x3ff0000000000000 if timing else 0],
            'research leaf differs from the original typed argument/default/coercion')
    require(values[0] != 0 and values[1] != 0 and B.round(1+values[0]) > 0,
            'outside nonzero-rate/period positive-base scope')
    future = values[3] if len(values) >= 4 else Rat(0)
    require(not (values[2] != 0 and future == -values[2]), 'opposing-balance control in scope')
    require(binary64(row['neg_log']) < 0 and binary64(row['leaves'][4]) != 0,
            'outside frozen finite discount helper scope')


def write_result(path, result):
    require(source_hashes() == LOADED_SOURCE_SHA256, 'checker source changed during this run')
    result['checker_sha256'] = LOADED_SOURCE_SHA256
    with open(path, 'x') as out:
        json.dump(result, out, indent=2)
        out.write('\n')
    print(json.dumps(result), flush=True)


@lru_cache(maxsize=65536)
def product_preimage(target):
    return E.pullback(B.pullback(Span(target, target)))


def quotient(raw, certificates, output):
    start = time.monotonic()
    expected_gaps = list(records(certificates))
    by_id = {r['id']: r for r in expected_gaps}
    require(len(by_id) == len(expected_gaps), 'duplicate legacy certificate ID')
    counts = Counter()
    reproduced = set()
    with open(raw) as source, open(str(output)+'.gaps.jsonl', 'x') as gaps:
        for line_number, line in enumerate(source, 1):
            columns = line.strip().split('\t')
            if columns[0] != 'Q':
                require(line_number <= 2, 'unexpected non-row in raw quotient trace')
                continue
            require(len(columns) == 6, 'malformed quotient row')
            # Original instrumenter writes Q,id,rate,quotient,EXPECTED,ACTUAL.
            _, identity, rate_bits, old_q, target_bits, actual_bits = columns
            r, y = abs(binary64(rate_bits)), abs(binary64(target_bits))
            require(r > 0 and y >= B.value(B.half), 'scope needs nonzero rate and normal target')
            interval = product_preimage(y).linear(1/r)
            possible = B.coordinates(interval)
            require(round_model(binary64(old_q)*binary64(rate_bits), 'stored_pc64') ==
                    binary64(actual_bits), 'original quotient/product forward replay failed')
            counts['rows'] += 1
            counts['original_quotient_product_replays'] += 1
            counts['old_match_controls'] += actual_bits == target_bits
            if possible is not None:
                q = B.value(possible[0])
                require(round_model(q*r, 'stored_pc64') == y, 'constructive reachability replay failed')
                counts['reachable_by_binary64_quotient'] += 1
                continue
            require(identity in by_id, f'unrecorded gap {identity}')
            require(identity not in reproduced, 'duplicate gap ID in raw records')
            reproduced.add(identity)
            old = by_id[identity]
            require(old['rate_bits'] == rate_bits and old['expected_bits'] == target_bits and
                    old['old_actual_bits'] == actual_bits and old['old_quotient_bits'] == old_q,
                    'certificate/raw bit identity mismatch')
            require(interval.low == Rat(*map(int, old['lower'])) and
                    interval.high == Rat(*map(int, old['upper'])) and
                    interval.left == interval.right == old['closed'], 'inverse endpoint mismatch')
            lower, upper = map(binary64, old['neighbor_magnitudes'])
            k = B.floor_coordinate(lower)
            require(B.value(k) == lower and B.value(k+1) == upper, 'nonadjacent legacy neighbors')
            require(lower <= interval.low and upper >= interval.high and
                    not interval.contains(lower) and not interval.contains(upper), 'gap bracket failed')
            for q in (lower, upper):
                require(round_model(q*r, 'stored_pc64') != y, 'legacy neighbor reaches target')
            finer = Lattice(54, -1022, 1023).coordinates(interval)
            require(finer is not None, '54-bit existence diagnostic failed')
            q54 = Lattice(54, -1022, 1023).value(finer[0])
            require(round_model(q54*r, 'stored_pc64') == y, '54-bit constructive diagnostic failed')
            counts['unreachable_by_binary64_quotient'] += 1
            counts['reachable_by_54_bit_quotient_diagnostic'] += 1
            counts['also_unreachable_by_single_binary64_product'] += (
                B.coordinates(B.pullback(Span(y, y)).linear(1/r)) is None)
            gaps.write(json.dumps(dict(id=identity, rate_bits=rate_bits, target_bits=target_bits,
                quotient_magnitude_interval=interval.json(), adjacent_magnitudes=old['neighbor_magnitudes'],
                constructive_54_bit_quotient=[str(q54.numerator), str(q54.denominator)]))+'\n')
    require(reproduced == set(by_id), 'missing legacy gaps')
    write_result(output, dict(complete=True, counts=counts, seconds=time.monotonic()-start,
        inputs={str(p): digest(p) for p in (raw, certificates)},
        claim='Independent exact rational verification of the original final shell RN53(RN64(q*r)), '
              'quantifying over every finite binary64 q, including signed rates by reflection. '
              'Reachability is not an Excel match or an executable quotient producer.',
        limits='Normal nonzero targets only; no assertion about other shells, growth arm, helper '
               'programs or Excel precision. Original observation transport remains an assumption.'))


def compile_tree(expression, models):
    position = 1  # site zero is the timing-factor add

    def walk(tree):
        nonlocal position
        if set(tree) == {'leaf'}:
            require(tree['leaf'] in range(4), 'unknown leaf')
            return tree['leaf']
        require(len(tree) == 1, 'malformed tree')
        op, children = next(iter(tree.items()))
        require(op in ('mul', 'div') and len(children) == 2, 'unsupported combine')
        left, right = map(walk, children)
        result = (op, left, right, models[position])
        position += 1
        return result

    tree = walk(expression)
    require(position == len(models) == 4, 'incorrect model arity')
    return tree


def helper_occurrences(tree):
    return int(tree == 2) if isinstance(tree, int) else sum(map(helper_occurrences, tree[1:3]))


def constant(tree, args):
    if isinstance(tree, int):
        require(tree != 2, 'helper used as a constant')
        return args[tree]
    op, left, right, model = tree
    return operation(constant(left, args), constant(right, args), op, model)


def undo_tree(tree, args, span):
    if span is None:
        return None
    if isinstance(tree, int):
        require(tree == 2, 'inverse did not end at helper')
        return span
    op, left, right, model = tree
    where = helper_occurrences(left) != 0
    require(helper_occurrences(tree) == 1, 'only one helper occurrence supported')
    span = inverse(span, constant(right if where else left, args), op, where, model)
    return undo_tree(left if where else right, args, span)


def group_row(row, tree, models):
    rate, pv, fv, timing, _ = map(binary64, row['leaves'])
    require(fv == 0 and pv != 0 and rate != 0 and timing in (0, 1), 'wrong projection scope')
    target = binary64(row['raw']['expected'])
    require(abs(target) >= B.value(B.half), 'requires a normal finite nonzero target')
    tf = operation(Rat(1), rate*timing, 'add', models[0])
    span = B.pullback(Span(target, target))  # final publication, including resident roots
    return undo_tree(tree, [pv, rate, None, tf], span)


def validate_combine_catalogue(expressions, profiles, count):
    """Reconstruct the finite grammar from leaf exponents, not search IDs."""
    def trees(leaves):
        if len(leaves) == 1:
            yield {'leaf': leaves[0]}, tuple(int(i == leaves[0]) for i in range(4))
            return
        for split in range(1, len(leaves)):
            for a, ae in trees(leaves[:split]):
                for b, be in trees(leaves[split:]):
                    for op, sign in (('mul', 1), ('div', -1)):
                        yield {op: [a, b]}, tuple(x+sign*y for x, y in zip(ae, be))
    canonical = lambda x: json.dumps(x, sort_keys=True, separators=(',', ':'))
    expected = {canonical(t) for order in itertools.permutations(range(4))
                for t, powers in trees(order) if powers == (1, 1, -1, -1)}
    require(len(expected) == len(expressions) == 60 and
            set(map(canonical, expressions)) == expected, 'four-leaf association grammar incomplete')
    names = ('strict', 'stored_pc64', 'resident_pc64')
    if count == 37500:
        names += ('stored_pc53', 'resident_pc53')
    require(count == 60*len(names)**4, 'undeclared grammar size')
    expected_profiles = {(e, ms) for e in range(60) for ms in itertools.product(names, repeat=4)}
    actual = {(c['expression'], tuple(c['models'])) for c in profiles}
    require(actual == expected_profiles and len(actual) == len(profiles), 'model assignment coverage incomplete')


def groups(scoped, catalogue, certificates, output, count):
    start = time.monotonic()
    expressions = json.loads((Path(catalogue)/'expressions.json').read_text())
    profiles = json.loads((Path(catalogue)/'combines.json').read_text())
    validate_combine_catalogue(expressions, profiles, count)
    lookup = {c['id']: c for c in profiles}
    require(len(lookup) == len(profiles) == count, 'catalogue coverage mismatch')
    needed = {w['trace_index'] for p in certificates for c in records(p) for w in c['proof']['witnesses']}
    rows, helper_bits, group_counts = {}, {}, Counter()
    frozen_group_violations = set()
    scoped_rows = 0
    for row in records(scoped):
        validate_row(row)
        if row['trace_index'] in needed:
            require(row['trace_index'] not in rows, 'duplicate trace index')
            rows[row['trace_index']] = row
        if binary64(row['leaves'][2]) != 0:
            continue
        raw = row['raw']
        key = (row['leaves'][0], raw['args'][1], json.dumps(raw['reading'], sort_keys=True))
        if key in helper_bits and helper_bits[key] != row['leaves'][4]:
            frozen_group_violations.add(key)
        helper_bits[key] = row['leaves'][4]
        group_counts[key] += 1
        scoped_rows += 1
    require(set(rows) == needed, 'certificate witness missing from original frozen scope')
    results = []
    for path in certificates:
        seen, reasons, witness_sets = set(), Counter(), set()
        multi_source_certificates = 0
        for cert in records(path):
            cid = cert['candidate_id']
            require(cid not in seen and cert['combine'] == lookup[cid], 'profile identity mismatch')
            require(cert['expression'] == expressions[lookup[cid]['expression']], 'expression mismatch')
            seen.add(cid)
            models = cert['combine']['models']
            tree = compile_tree(cert['expression'], models)
            current, first = None, True
            frozen_allowed = True
            indices = []
            for entry in cert['proof']['witnesses']:
                row = rows[entry['trace_index']]
                require(row['raw'] == entry['witness'], 'original witness bytes/fields changed')
                raw = row['raw']
                require(row['leaves'][0] == int(cert['rate'], 16) and
                        raw['args'][0] == cert['rate'] and raw['args'][1] == cert['nper'] and
                        raw['reading'] == cert['environment'], 'group identity mismatch')
                if cert['canonical_timing'] is not None:
                    require(row['leaves'][3] == cert['canonical_timing'], 'timing group mismatch')
                if Path(path).name.startswith('nonunit'):
                    require(binary64(raw['args'][1]) != 1, 'unit row in nonunit certificate')
                require(binary64(raw['args'][2]) == binary64(row['leaves'][1]), 'PV leaf mismatch')
                span = group_row(row, tree, models)
                frozen_allowed &= span is not None and span.contains(binary64(row['leaves'][4]))
                current = span if first else current.meet(span) if current is not None else None
                first = False
                indices.append(row['trace_index'])
            require(not first and len(indices) <= 3, 'unsupported sufficient-subset size')
            require(current is None or E.coordinates(current) is None, 'certificate does not exclude helper')
            reason = 'empty_real_intersection' if current is None else 'empty_ext80_lattice'
            require(reason == cert['proof']['reason'], 'reason mismatch')
            require((current.json() if current is not None else None) == cert['proof']['remaining_interval'],
                    'remaining interval mismatch')
            require(not frozen_allowed, 'excluded group still allows the frozen helper')
            multi_source_certificates += len({w['witness']['source'] for w in cert['proof']['witnesses']}) > 1
            reasons[reason] += 1
            witness_sets.add(tuple(sorted(indices)))
        require(seen == set(lookup), 'incomplete profile certificate coverage')
        results.append(dict(certificates=str(path), sha256=digest(path), checked_profiles=len(seen),
            reasons=reasons, sufficient_witness_sets=len(witness_sets), survivors=0))
        results[-1]['sufficient_subsets_spanning_multiple_sources'] = multi_source_certificates
        print(json.dumps(results[-1]), flush=True)
    write_result(output, dict(complete=True, results=results, seconds=time.monotonic()-start,
        original_FV0_scope_rows=scoped_rows, exact_rate_nper_environment_groups=len(group_counts),
        frozen_binary64_helper_inconsistent_groups=len(frozen_group_violations),
        inputs={str(p): digest(p) for p in (scoped, Path(catalogue)/'expressions.json', Path(catalogue)/'combines.json')},
        claim='Every declared combine profile has an independently checked sufficient original group '
              'with no common finite Ext80 helper; each also excludes the frozen binary64 helper. '
              'Real-empty intersections additionally exclude arbitrary real helpers.',
        limits='FV=0 projection, specified multiply/divide trees and stores, normal killing targets; '
               'does not exclude amount-dependent helpers, other operations/branches or publication semantics. '
               'Frozen helper sharing is checked as research data, not asserted as Excel fact.'))


def tagged_value(text):
    if text.startswith('binary64:'):
        return binary64(text.split(':')[1])
    require(text.startswith('ext80-le:'), 'unknown full-width value encoding')
    raw = bytes.fromhex(text.split(':')[1])
    require(len(raw) == 10, 'wrong extended byte count')
    significand = int.from_bytes(raw[:8], 'little')
    se = int.from_bytes(raw[8:], 'little')
    exponent = se & 0x7fff
    require(exponent != 0x7fff, 'nonfinite extended value outside this verifier scope')
    require(bool(significand >> 63) == (exponent != 0), 'noncanonical extended encoding')
    return (-1 if se >> 15 else 1)*significand*two((exponent-16383 if exponent else -16382)-63)


def graph_shape(graph):
    nodes = []
    for index, node in enumerate(graph['nodes']):
        op = node['op']
        require(len(op) == 1, 'malformed graph operation')
        kind, operands = next(iter(op.items()))
        if kind == 'arg':
            require(operands in range(5), 'unknown argument')
            shape = ('arg', operands)
        elif kind == 'const':
            require(binary64(operands['bits']) == 1, 'undeclared constant')
            shape = ('one',)
        else:
            require(kind in ('add', 'mul', 'div') and len(operands) == 2 and
                    all(0 <= n < index for n in operands), 'invalid graph references')
            shape = (kind, *(nodes[n] for n in operands))
        nodes.append(shape)
    return nodes[graph['output']]


def expected_shape(numerator, combine):
    r, pv, fv, timing, em = [('arg', n) for n in range(5)]
    one = ('one',)
    add = lambda a, b: ('add', a, b)
    mul = lambda a, b: ('mul', a, b)
    div = lambda a, b: ('div', a, b)
    tf = add(one, mul(r, timing))
    if numerator is None:
        require(combine == 8, 'unknown rearrangement')
        return mul(add(fv, div(add(pv, fv), em)), div(r, tf))
    n = [add(pv, mul(fv, add(one, em))), add(add(pv, fv), mul(fv, em)),
         add(pv, add(fv, mul(fv, em))), add(add(pv, mul(fv, em)), fv)][numerator]
    return [mul(div(div(n, em), tf), r), div(mul(n, r), mul(em, tf)),
            mul(div(n, tf), div(r, em)), div(div(mul(n, r), em), tf),
            div(mul(div(n, em), r), tf), mul(div(n, em), div(r, tf)),
            mul(n, div(r, mul(em, tf))), mul(div(n, mul(em, tf)), r)][combine]


def graph_trace(graph, args):
    values = []
    for node in graph['nodes']:
        kind, operands = next(iter(node['op'].items()))
        if kind == 'arg':
            value = args[operands]
        elif kind == 'const':
            value = binary64(operands['bits'])
        else:
            value = operation(*(values[n] for n in operands), kind, node['model'])
        values.append(value)
    return values


def graphs(scoped, programs, kills, output):
    """Validate grammar coverage and every full-width first-kill trace."""
    start = time.monotonic()
    locations, seen, shapes = [], set(), {}
    models = ('strict', 'stored_pc64', 'resident_pc64')
    with open(programs, 'rb') as source:
        while True:
            location = source.tell()
            line = source.readline()
            if not line:
                break
            c = json.loads(line)
            require(c['id'] == len(locations), 'noncanonical candidate IDs')
            graph = c['graph']
            template_key = (c['numerator'], c['combine'])
            shape = shapes.setdefault(template_key, expected_shape(*template_key))
            require(graph_shape(graph) == shape, 'graph differs from declared association')
            require(c['models'] == [n['model'] for n in graph['nodes']
                if next(iter(n['op'])) in ('add', 'mul', 'div')], 'model metadata mismatch')
            expected_sites = 7 if c['numerator'] is None else 8
            require(len(c['models']) == expected_sites and all(m in models for m in c['models']),
                    'invalid arithmetic annotation')
            key = (*template_key, tuple(c['models']))
            require(key not in seen, 'duplicate annotated template')
            seen.add(key)
            locations.append(location)
    require(len(shapes) == 33 and len(seen) == 32*3**8+3**7 == len(locations), 'grammar coverage incomplete')
    # Every template has exactly its full Cartesian model assignments.
    for template_key in shapes:
        actual = {key[2] for key in seen if key[:2] == template_key}
        require(len(actual) == 3**(7 if template_key[0] is None else 8), 'annotation coverage incomplete')
    del seen
    needed = {c['trace_index'] for c in records(kills)}
    rows = {}
    for row in records(scoped):
        if row['trace_index'] in needed:
            validate_row(row)
            require(row['trace_index'] not in rows, 'duplicate original trace index')
            rows[row['trace_index']] = row
    require(set(rows) == needed, 'missing original first-kill witness')
    checked, node_count, witness_counts = set(), 0, Counter()
    with open(programs, 'rb') as source:
        for cert in records(kills):
            cid = cert['candidate_id']
            require(cid not in checked and 0 <= cid < len(locations), 'invalid/duplicate first kill')
            source.seek(locations[cid])
            candidate = json.loads(source.readline())
            row = rows[cert['trace_index']]
            require(row['raw'] == cert['witness'] and row['raw']['expected'] == cert['expected'],
                    'changed first-kill observation')
            args = list(map(binary64, row['leaves']))
            trace = graph_trace(candidate['graph'], args)
            require(len(trace) == len(cert['trace']), 'incomplete full-width trace')
            for node, expected, encoded in zip(candidate['graph']['nodes'], trace, cert['trace']):
                kind = next(iter(node['op']))
                extended = kind in ('add', 'mul', 'div') and node['model'].startswith('resident')
                require(encoded.startswith('ext80-le:' if extended else 'binary64:'), 'value format mismatch')
                require(tagged_value(encoded) == expected, 'independent arithmetic trace mismatch')
                node_count += 1
            root = trace[candidate['graph']['output']]
            require(tagged_value(cert['root']) == root, 'root mismatch')
            published = B.round(root)
            target = binary64(cert['expected'])
            require(published != target, 'alleged killing witness actually matches')
            if published is None:
                got = 'error:Num'
            elif published != 0 and abs(published) < two(-1022):
                got = 'fault:NotImplemented'
            else:
                require(binary64(cert['got']) == published, 'published result mismatch')
                got = cert['got']
            require(got == cert['got'], 'publication policy mismatch')
            checked.add(cid)
            witness_counts[row['raw']['id']] += 1
    require(len(checked) == len(locations), 'missing candidate contradiction')
    write_result(output, dict(complete=True, grammar_programs=len(locations), checked_first_kills=len(checked),
        checked_full_width_nodes=node_count, witness_counts=witness_counts, seconds=time.monotonic()-start,
        inputs={str(p): digest(p) for p in (scoped, programs, kills)},
        claim='Every annotated ordinary-association program is contradicted by a trusted original witness; '
              'every stored/resident operation in its first-kill trace is independently checked as an exact rational.',
        limits='This bounded grammar only, frozen binary64 helper leaves. Intermediate zero signs are '
               'retained in traces but rational checking ignores them; no zero divisor occurs in the checked '
               'traces and PMT normalizes published zero. No conclusion about Excel internals or other arms.'))


def controls(path, output):
    count, cases = 0, Counter()
    for control in records(path):
        op = next(iter(control['op']))
        args = [tagged_value(control[k]) for k in ('a', 'b')]
        try:
            expected = operation(*args, op, control['model'])
        except ValueError as error:
            require(str(error) in ('nonfinite operand', 'nonfinite operation'), 'unsupported reference control')
            # All overflow controls here have positive operands and positive results.
            require(all(a > 0 for a in args), 'unjustified infinity sign')
            expected = None
        for backend in ('hard', 'soft'):
            if expected is None:
                require(control[backend] == 'binary64:7ff0000000000000', 'positive overflow mismatch')
            else:
                require(tagged_value(control[backend]) == expected, 'constructed reference mismatch')
            count += 1
        cases[control['id']] += 1
    require(len(cases) == 7 and all(n == 5 for n in cases.values()), 'incomplete reference suite')
    write_result(output, dict(complete=True, reference_comparisons=count, controls=cases,
        input_sha256=digest(path), claim='Both numeric backends validated against deliberate independent '
        'exact rational controls for double rounding, input stores, resident cancellation, spills, '
        'extended exponent range, binary64 underflow ties and negative division.',
        limits='Finite suite; zero signs ignored by rational comparison. No global backend equivalence '
               'or transcendental equivalence claim. Does not cover the known PC53 extended-exponent underflow discrepancy.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    q = sub.add_parser('quotient')
    for name in ('raw', 'certificates', 'output'):
        q.add_argument(name, type=Path)
    g = sub.add_parser('groups')
    for name in ('scoped', 'catalogue', 'output'):
        g.add_argument(name, type=Path)
    g.add_argument('certificates', nargs='+', type=Path)
    g.add_argument('--count', type=int, required=True)
    t = sub.add_parser('graphs')
    for name in ('scoped', 'programs', 'kills', 'output'):
        t.add_argument(name, type=Path)
    c = sub.add_parser('controls')
    c.add_argument('input', type=Path)
    c.add_argument('output', type=Path)
    a = parser.parse_args()
    if a.command == 'quotient':
        quotient(a.raw, a.certificates, a.output)
    elif a.command == 'groups':
        groups(a.scoped, a.catalogue, a.certificates, a.output, a.count)
    elif a.command == 'graphs':
        graphs(a.scoped, a.programs, a.kills, a.output)
    else:
        controls(a.input, a.output)


if __name__ == '__main__':
    main()
