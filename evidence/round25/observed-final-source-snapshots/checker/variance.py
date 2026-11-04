"""Independent rational verifier of the R14/R16 finite variance families.

No search generator, native interpreter, product source or numeric library is
imported. The only shared arithmetic is the existing integer/rational checker.
The program grammar is derived here from the unit's mathematical identities.
"""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
import itertools
import json
from pathlib import Path
import time

from exact import BINARY64 as B, binary64, operation, two
from check import require, records, digest, tagged_value
from normal import sqrt_round

MODELS = ('strict', 'stored_pc64', 'resident_pc64', 'stored_pc53', 'resident_pc53')
COUNTS = Counter()


def profiles(function=None, family='r16'):
    if family == 'structural-r17':
        identity = 0
        for algorithm in range(1 if function == 'AVERAGE' else 22):
            if algorithm in (5, 6, 19, 20, 21):
                orders = [(0, 0), (1, 1)]
            elif function == 'AVERAGE':
                orders = [(i, 0) for i in range(13)]
            else:
                orders = list(itertools.product(range(3), repeat=2))+[(i, i) for i in range(3, 13)]
            for linear, squared in orders:
                for code in range(25 if function == 'AVERAGE' else 3125):
                    models = [MODELS[(code//5**i) % 5] if function != 'AVERAGE' or i < 2 else 'strict'
                              for i in range(5)]
                    yield dict(id=identity, algorithm=algorithm, linear_order=linear,
                               square_order=squared, models=models)
                    identity += 1
        return
    require(family == 'r16', 'unknown variance family')
    if function == 'AVERAGE':
        identity = 0
        for order, division, reduction in itertools.product(range(3), MODELS, MODELS):
            yield dict(id=identity, algorithm=0, linear_order=order, square_order=0,
                       models=[reduction, division, 'strict', 'strict', 'strict'])
            identity += 1
        return
    identity = 0
    for algorithm in range(14):
        orders = itertools.product(range(3), repeat=2) if algorithm < 5 or algorithm >= 7 else ((0, 0), (1, 1))
        for linear, square in orders:
            for reversed_models in itertools.product(MODELS, repeat=5):
                yield dict(id=identity, algorithm=algorithm, linear_order=linear,
                           square_order=square, models=list(reversed_models[::-1]))
                identity += 1


@lru_cache(maxsize=262144)
def arithmetic(kind, model, x, y=Q(0)):
    COUNTS['distinct_exact_operations'] += 1
    if kind == 'sqrt':
        # A negative intermediate is a contradiction for a numeric target;
        # no negative value is clipped, tolerated, or replaced by an answer.
        return None if x < 0 else sqrt_round(x, model)
    return operation(x, y, kind, model)


def evaluate(profile, values, function):
    """Re-derive the fourteen mathematical programs; all rounding is rational."""
    modes = profile['models']
    def calc(region, kind, x, y=Q(0)):
        return arithmetic(kind, modes[region], x, y)
    def subtract(region, x, y):
        return calc(region, 'add', x, -y)
    def reduction(items, order, region):
        if order == 2:
            if not items: return Q(0)
            if len(items) == 1: return items[0]
            middle = len(items)//2
            return calc(region, 'add', reduction(items[:middle], 2, region),
                        reduction(items[middle:], 2, region))
        if order in (3, 4):
            return reduction(sorted(items, key=lambda x: abs(B.round(x))), int(order == 4), region)
        if order in (5, 6):
            total = error = Q(0)
            for x in items if order == 5 else items[::-1]:
                adjusted = subtract(region, x, error)
                updated = calc(region, 'add', total, adjusted)
                error = subtract(region, subtract(region, updated, total), adjusted)
                total = updated
            return total
        if order in (7, 8):
            total = error = Q(0)
            for x in items if order == 7 else items[::-1]:
                updated = calc(region, 'add', total, x)
                larger, smaller = (total, x) if abs(B.round(total)) >= abs(B.round(x)) else (x, total)
                lost = calc(region, 'add', subtract(region, larger, updated), smaller)
                error = calc(region, 'add', error, lost)
                total = updated
            return calc(region, 'add', total, error)
        if order in (9, 10, 11):
            lanes = [Q(0)]*(2**(order-8))
            for i, x in enumerate(items):
                lanes[i % len(lanes)] = calc(region, 'add', lanes[i % len(lanes)], x)
            return reduction(lanes, 0, region)
        if order == 12:
            items = list(items)
            while len(items) > 1:
                items = [calc(region, 'add', items[i], items[i+1]) if i+1 < len(items) else items[i]
                         for i in range(0, len(items), 2)]
            return items[0] if items else Q(0)
        total = Q(0)
        for item in items if order == 0 else items[::-1]:
            total = calc(region, 'add', total, item)
        return total
    if function == 'AVERAGE':
        return calc(1, 'div', reduction(values, profile['linear_order'], 0), Q(len(values)))
    algorithm = profile['algorithm']
    linear, squared = profile['linear_order'], profile['square_order']
    size = Q(len(values))
    if algorithm >= 14:
        if algorithm == 21:
            def merge(data):
                if len(data) == 1: return data[0], Q(0)
                middle = len(data)//2
                left, lss = merge(data[:middle])
                right, rss = merge(data[middle:])
                difference = subtract(1, right, left)
                weight = calc(1, 'div', Q(len(data)-middle), Q(len(data)))
                center = calc(1, 'add', left, calc(1, 'mul', difference, weight))
                square = calc(2, 'mul', difference, difference)
                extra = calc(4, 'div', calc(2, 'mul', square, Q(middle*(len(data)-middle))), Q(len(data)))
                return center, calc(3, 'add', calc(3, 'add', lss, rss), extra)
            dispersion = merge(values if linear == 0 else values[::-1])[1]
        elif algorithm in (19, 20):
            data = values if linear == 0 else values[::-1]
            mean, dispersion = data[0], Q(0)
            for count, x in enumerate(data[1:], 2):
                difference = subtract(1, x, mean)
                if algorithm == 19:
                    mean = calc(1, 'div', calc(1, 'add', calc(1, 'mul', mean, Q(count-1)), x), Q(count))
                    term = calc(2, 'mul', difference, subtract(1, x, mean))
                else:
                    mean = calc(1, 'add', mean, calc(1, 'div', difference, Q(count)))
                    term = calc(4, 'div', calc(2, 'mul', calc(2, 'mul', difference, difference), Q(count-1)), Q(count))
                dispersion = calc(3, 'add', dispersion, term)
        else:
            data = [subtract(1, x, values[0]) for x in values] if algorithm == 15 else values
            mean = calc(1, 'div', reduction(data, linear, 0), size)
            for _ in range(1 if algorithm == 14 else 2 if algorithm == 18 else 0):
                errors = [subtract(1, x, mean) for x in data]
                mean = calc(1, 'add', mean, calc(1, 'div', reduction(errors, linear, 0), size))
            errors = [subtract(1, x, mean) for x in data]
            terms = [calc(2, 'mul', delta, x if algorithm in (16, 17) else delta)
                     for delta, x in zip(errors, data)]
            dispersion = reduction(terms, squared, 3)
            if algorithm in (15, 17):
                residual = reduction(errors, linear, 0)
                correction = calc(4, 'div', calc(2, 'mul', residual, residual), size) if algorithm == 15 else calc(2, 'mul', mean, residual)
                dispersion = subtract(4, dispersion, correction)
        population = function in ('VAR.P', 'VARP', 'VARPA', 'STDEV.P', 'STDEVP', 'STDEVPA')
        result = calc(4, 'add', dispersion, Q(0)) if function == 'DEVSQ' else calc(4, 'div', dispersion, size if population else size-1)
        return calc(4, 'sqrt', result) if function.startswith('STDEV') else result
    if algorithm >= 7:
        total = reduction(values, linear, 0)
        squares = reduction([calc(2, 'mul', x, x) for x in values], squared, 3)
        square_total = calc(2, 'mul', total, total)
        population = function in ('VAR.P', 'VARP', 'VARPA', 'STDEV.P', 'STDEVP', 'STDEVPA')
        divisor = Q(1) if function == 'DEVSQ' else (size if population else size-1)
        if algorithm == 7:
            numerator = subtract(4, squares, calc(1, 'div', square_total, size))
            result = numerator if function == 'DEVSQ' else calc(4, 'div', numerator, divisor)
        elif algorithm == 8:
            numerator = subtract(4, calc(4, 'mul', size, squares), square_total)
            denominator = size if function == 'DEVSQ' else calc(4, 'mul', size, divisor)
            result = calc(4, 'div', numerator, denominator)
        elif algorithm == 9:
            left = calc(4, 'div', squares, divisor)
            right = calc(4, 'div', square_total, calc(4, 'mul', size, divisor))
            result = subtract(4, left, right)
        elif algorithm == 10:
            left = calc(4, 'div', squares, divisor)
            mean = calc(1, 'div', total, size)
            scaled = calc(1, 'div', total, divisor)
            result = subtract(4, left, calc(2, 'mul', mean, scaled))
        elif algorithm == 11:
            mean = calc(1, 'div', total, size)
            difference = subtract(4, calc(4, 'div', squares, size), calc(2, 'mul', mean, mean))
            scale = size if function == 'DEVSQ' else calc(4, 'div', size, divisor)
            result = calc(4, 'mul', difference, scale)
        elif algorithm == 12:
            mean = calc(1, 'div', total, size)
            terms = []
            for x in values:
                difference = subtract(1, x, mean)
                product = calc(2, 'mul', difference, difference)
                terms.append(product if function == 'DEVSQ' else calc(4, 'div', product, divisor))
            result = reduction(terms, squared, 3)
        else:
            terms = []
            for i in range(len(values)):
                for j in range(i+1, len(values)):
                    difference = subtract(1, values[i], values[j])
                    terms.append(calc(2, 'mul', difference, difference))
            numerator = reduction(terms, squared, 3)
            denominator = size if function == 'DEVSQ' else calc(4, 'mul', size, divisor)
            result = calc(4, 'div', numerator, denominator)
        return calc(4, 'sqrt', result) if function.startswith('STDEV') else result
    if algorithm in (5, 6):
        data = values if linear == 0 else values[::-1]
        mean, total, dispersion = data[0], data[0], Q(0)
        for count, value in enumerate(data[1:], 2):
            k = Q(count)
            if algorithm == 5:
                difference = subtract(1, value, mean)
                increment = calc(1, 'div', difference, k)
                mean = calc(1, 'add', mean, increment)
                remainder = subtract(1, value, mean)
                term = calc(2, 'mul', difference, remainder)
            else:
                total = calc(0, 'add', total, value)
                difference = subtract(1, calc(1, 'mul', k, value), total)
                product = calc(2, 'mul', difference, difference)
                denominator = calc(4, 'mul', k, k-1)
                term = calc(4, 'div', product, denominator)
            dispersion = calc(3, 'add', dispersion, term)
    else:
        data = values
        if algorithm == 4:
            data = [subtract(1, x, values[0]) for x in values]
        total = reduction(data, linear, 0)
        mean = calc(1, 'div', total, size)
        if algorithm in (0, 1):
            differences = [subtract(1, x, mean) for x in data]
            products = [calc(2, 'mul', x, x) for x in differences]
            dispersion = reduction(products, squared, 3)
            if algorithm == 1:
                error = reduction(differences, squared, 0)
                correction = calc(4, 'div', calc(2, 'mul', error, error), size)
                dispersion = subtract(4, dispersion, correction)
        else:
            products = [calc(2, 'mul', x, x) for x in data]
            uncentered = reduction(products, squared, 3)
            if algorithm == 3:
                scaled = calc(4, 'mul', size, uncentered)
                correction = calc(2, 'mul', total, total)
                dispersion = calc(4, 'div', subtract(4, scaled, correction), size)
            else:
                correction = calc(2, 'mul', total, mean)
                dispersion = subtract(4, uncentered, correction)
    if function == 'DEVSQ':
        return dispersion
    is_population = function in ('VAR.P', 'VARP', 'VARPA', 'STDEV.P', 'STDEVP', 'STDEVPA')
    variance = calc(4, 'div', dispersion, size if is_population else size-1)
    return calc(4, 'sqrt', variance) if function.startswith('STDEV') else variance


def decode_scoped(row):
    raw = row['raw']
    require(raw['result_index'] is None, 'selected result outside grammar')
    bits = []
    for arg in raw['args']:
        if isinstance(arg, list):
            items = list(itertools.chain.from_iterable(arg)) if all(isinstance(r, list) for r in arg) else arg
        else:
            items = [arg]
        for item in items:
            require(isinstance(item, str) and item.startswith('0x') and len(item) == 18,
                    'non-numeric original input in numeric scope')
            bits.append(int(item, 16))
    require(bits == row['values_bits'], 'derived input differs from original typed input bits')
    require(2 <= len(bits) <= 32, 'length outside grammar')
    values = list(map(binary64, bits))
    require(all(x == 0 or two(-400) <= abs(x) <= two(400) for x in values), 'input outside bounded scope')
    expected = binary64(raw['expected'])
    require((expected >= 0 or raw['function'] == 'AVERAGE') and (expected == 0 or abs(expected) >= two(-1022)), 'target outside scope')
    return values


def outcome(root):
    if root is None:
        return 'error:Num'
    result = B.round(root)
    if result is None:
        return 'error:Num'
    if result == 0:
        return '0x0000000000000000'
    if abs(result) < two(-1022):
        return 'fault:NotImplemented'
    return f'0x{B.floor_coordinate(abs(result)) | ((1 << 63) if result < 0 else 0):016x}'


def verify_root(root, encoded, model):
    require(encoded.startswith('ext80-le:') if model.startswith('resident') else encoded.startswith('binary64:'),
            'incorrect root representation for store policy')
    if root is None:
        if encoded.startswith('binary64:'):
            bits = int(encoded.split(':')[1], 16)
            require((bits >> 52) & 2047 == 2047, 'negative root did not produce a nonfinite result')
        else:
            bits = bytes.fromhex(encoded.split(':')[1])
            require(int.from_bytes(bits[8:], 'little') & 32767 == 32767, 'negative root did not produce nonfinite Ext80')
    else:
        require(tagged_value(encoded) == root, 'exact full-width root mismatch')


def root_profile_key(profile, function):
    """Ignore only the label and provably dead DEVSQ expression fields.

    DEVSQ forms 0 and 12 return the identical centered sum, with no divisor.
    Their region 4 does not reach that root. All logical labels remain in the
    catalogue, first-kill coverage and execution counts.
    """
    algorithm=profile['algorithm'];models=tuple(profile['models'])
    if function=='DEVSQ' and algorithm in (0,12):
        algorithm=0;models=(*models[:4],'strict')
    return algorithm,profile['linear_order'],profile['square_order'],models


def check(scoped, directory, output):
    start = time.monotonic()
    source = Path(__file__)
    sources = {str(p): digest(p) for p in (source, source.with_name('exact.py'), source.with_name('check.py'), source.with_name('normal.py'))}
    directory = Path(directory)
    summary = json.loads((directory/'summary.json').read_text())
    require(summary['complete'] is True, 'partial search cannot certify completeness')
    function = summary['function']
    require(function in ('VAR','VAR.S','VAR.P','VARP','VARA','VARPA','STDEV','STDEV.S',
            'STDEV.P','STDEVP','STDEVA','STDEVPA','DEVSQ','AVERAGE'),'unknown variance function')
    frozen = [r for r in records(scoped) if r['raw']['function'] == function]
    require(frozen, 'empty variance observation scope')
    values = [decode_scoped(r) for r in frozen]
    ps = list(records(directory/'profiles.jsonl'))
    family = summary.get('family', 'r16')
    require(ps == list(profiles(function, family)) and len(ps) == summary['programs'], 'incomplete or modified grammar')
    keys=[root_profile_key(p,function)for p in ps]
    @lru_cache(maxsize=65536)
    def replay_root(key, at):
        algorithm,linear,squared,models=key
        profile=dict(algorithm=algorithm,linear_order=linear,square_order=squared,models=models)
        return evaluate(profile,values[at],function)
    witnesses = list(records(directory/'witnesses.jsonl'))
    require(len(witnesses) == len(frozen) == summary['observations'], 'missing observation accounting')
    seen = set()
    kill_counts = Counter()
    counts = Counter()
    for kill in records(directory/'first-kills.jsonl'):
        identity, order = kill['candidate_id'], kill['order']
        require(identity not in seen and 0 <= identity < len(ps), 'duplicate or unknown candidate')
        require(0 <= order < len(frozen), 'unknown witness order')
        require(kill['trace_index'] == frozen[order]['trace_index'], 'changed witness identity')
        require(kill['expected'] == frozen[order]['raw']['expected'], 'changed answer')
        p = ps[identity]
        for prior in range(order+1):
            root = replay_root(keys[identity], prior)
            calculated = outcome(root)
            expected = frozen[prior]['raw']['expected']
            counts['program_executions'] += 1
            if prior < order:
                require(calculated == expected, 'record is not the first contradiction')
            else:
                region = 1 if function == 'AVERAGE' else 4 if function.startswith('STDEV') else 3 if p['algorithm'] == 12 or (function == 'DEVSQ' and p['algorithm'] in (0, 5, 6)) else 4
                verify_root(root, kill['root'], p['models'][region])
                require(calculated == kill['got'] and calculated != expected, 'invalid kill')
        seen.add(identity)
        kill_counts[order] += 1
        counts['certified_first_kills'] += 1
    survivors = json.loads((directory/'survivors.json').read_text())
    require(len(set(survivors)) == len(survivors), 'duplicate survivor')
    require(seen.isdisjoint(survivors) and seen | set(survivors) == set(range(len(ps))), 'unsearched program')
    for at, row in enumerate(frozen):
        for identity in survivors:
            require(outcome(replay_root(keys[identity], at)) == row['raw']['expected'], 'invalid survivor')
            counts['program_executions'] += 1
    remaining = len(ps)
    for order, (record, row) in enumerate(zip(witnesses, frozen)):
        remaining -= kill_counts[order]
        require(record == dict(order=order, trace_index=row['trace_index'], killed=kill_counts[order], remaining=remaining),
                'incorrect witness coverage accounting')
    require(remaining == len(survivors) == summary['survivors'], 'incorrect survivor count')
    require(counts['program_executions'] == summary['executions'], 'execution count mismatch')
    require(all(digest(p) == h for p, h in sources.items()), 'checker source changed during proof')  # ubs:ignore — public integrity hashes, not secrets
    result = dict(complete=True, function=function, programs=len(ps), observations=len(frozen),
                  survivors=len(survivors), counts=counts, arithmetic_cache=arithmetic.cache_info()._asdict(),
                  distinct_exact_operations=COUNTS['distinct_exact_operations'], seconds=time.monotonic()-start,
                  root_cache=replay_root.cache_info()._asdict(),
                  source_sha256=sources,
                  input_sha256={str(p): digest(p) for p in (Path(scoped), *sorted(directory.iterdir())) if p.is_file()},
                  limits='Finite declared programs only. Coherent precision regions and original transport assumptions. '
                         'Zero signs are immaterial here: every divisor is a positive integer, roots publish zero with the product zero policy. '
                         'No alias identity, private Excel precision, extra branch, fitted cutoff or complete kernel is inferred.')
    with open(output, 'x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('scoped')
    parser.add_argument('directory')
    parser.add_argument('output')
    args = parser.parse_args()
    check(args.scoped, args.directory, args.output)
