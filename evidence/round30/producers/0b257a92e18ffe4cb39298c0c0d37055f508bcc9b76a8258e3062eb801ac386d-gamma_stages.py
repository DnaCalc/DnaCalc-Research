"""Paired public gamma calls, with one prefactor and one series stage.

Original MIT research code. Reuses the scalar recorder, ISA interpreter and
exact inverse arithmetic. Candidate labels are fixed before reading answers;
counts report contradictions and never select coefficients or a correction.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import sys

from fingerprint import (CachedService, bits, value, load_service, Native,
                         primitive_stage_graph, positive_extended_cell, Builder, sha, write)
from exact import BINARY64, EXTENDED64, Span, binary64, inverse
from x87 import load64, decode

STAGE_BINDINGS = {str(p.resolve()):sha(p) for p in (
    Path(__file__), Path(__file__).with_name('x87_programs.py'),
    Path(__file__).with_name('x87_native.c'),
    *(Path(__file__).parents[1]/'checker'/name for name in
      ('geomean.py', 'exact.py', 'check.py', 'x87.py')))}


def paired(checkout, output, cache, library, count=32, offline=False):
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    cw = 0x133f
    one = load64(int(bits(1), 16))
    exp_graph = primitive_stage_graph('EXP')

    def leaf(op, a, b):
        return native.leaf(op, a, b, cw)

    def store(a):
        return load64(native.publish(a))

    def scalar(a):
        return load64(int(bits(a), 16))

    def exp(a, stored=True):
        negative = bool(a[9] & 128)
        if negative:
            a = leaf('neg', a, a)
        result = native.evaluate(exp_graph, [], cw, extended=a)[:10]
        result = leaf('div', one, result) if negative else result
        return store(result) if stored else result

    # Adjacent windows and power-of-two edges; no current match status is read.
    shapes = (1.25, 1.5, 1.75, 2., 2.5, 3., 4., 6., 10., 19.)
    centers = (.125, .5, .7, 1., 1.1, 1.9, 2.5, 4., 8., 16.)
    cases = []
    for a in shapes:
        for center in centers:
            if center >= a:
                continue
            start = int(bits(center), 16)-count//2
            cases.extend((a, f'{start+k:016x}') for k in range(count))
    xs = sorted({x for _, x in cases})
    logs = dict(zip(xs, service.ask('LN', [[client.number(x)] for x in xs])['rows']))
    lg = dict(zip(shapes, service.ask('GAMMALN', [[client.number(a+1)] for a in shapes])['rows']))
    lg_density = dict(zip(shapes, service.ask('GAMMALN', [[client.number(a)] for a in shapes])['rows']))
    norm_gamma = dict(zip(shapes, service.ask('GAMMA', [[client.number(a+1)] for a in shapes])['rows']))
    questions = [[client.number(x), client.number(a), client.number(1.), client.boolean(flag)]
                 for a, x in cases for flag in (False, True)]
    answers = service.ask('GAMMA.DIST', questions)['rows']

    @lru_cache(maxsize=8192)
    def series(a, x, model):
        rr, term, summ = scalar(a), one, one
        terms = []
        if model in ('sum_zero', 'reverse'):
            summ = scalar(0.)
        xx = load64(int(x, 16))
        for iteration in range(1, 1000):
            rr = store(leaf('add', rr, one))
            if model == 'strict':
                t = value(f'{native.publish(term):016x}')
                z = value(f'{native.publish(rr):016x}')
                term = scalar(t*(value(x)/z))
                summ = scalar(value(f'{native.publish(summ):016x}') + value(f'{native.publish(term):016x}'))
            elif model in ('stored_quotient', 'resident_quotient', 'resident_sum', 'resident_term_sum', 'sum_zero', 'reverse'):
                quotient = leaf('div', xx, rr)
                if model in ('stored_quotient', 'sum_zero', 'reverse'):
                    quotient = store(quotient)
                term = leaf('mul', term, quotient)
                if model != 'resident_term_sum':
                    term = store(term)
                summ = leaf('add', summ, term)
                if model not in ('resident_sum', 'resident_term_sum'):
                    summ = store(summ)
                terms.append(term)
            elif model == 'multiply_then_divide':
                term = store(leaf('div', leaf('mul', term, xx), rr))
                summ = store(leaf('add', summ, term))
            else:
                raise ValueError('unknown isolated series model')
            # Fixed Cephes stopping criterion; no fitted iteration bound.
            test_sum = leaf('add', summ, one) if model in ('sum_zero', 'reverse') else summ
            if decode(term)/decode(test_sum) <= binary64(bits(2.**-53)):
                if model == 'reverse':
                    summ = scalar(0.)
                    for t in reversed(terms):
                        summ = store(leaf('add', summ, t))
                if model in ('sum_zero', 'reverse'):
                    summ = store(leaf('add', summ, one))
                return store(summ), iteration
        raise ValueError('series did not converge')

    totals, by_shape, records = Counter(), {}, []
    models = ('strict', 'stored_quotient', 'resident_quotient', 'multiply_then_divide',
              'resident_sum', 'resident_term_sum', 'sum_zero', 'reverse')
    for i, (a, x) in enumerate(cases):
        density, cdf = answers[2*i:2*i+2]
        if any(r['type'] != 'number' for r in (density, cdf, logs[x], lg[a])):
            raise ValueError('positive finite paired gamma call was not numeric')
        aa, xx = scalar(a), load64(int(x, 16))
        t = store(leaf('mul', aa, load64(int(logs[x]['bits'], 16))))
        argument = leaf('sub', leaf('sub', t, xx), load64(int(lg[a]['bits'], 16)))
        r = exp(argument)
        predicted_pdf = native.publish(leaf('div', store(leaf('mul', r, aa)), xx))
        density_consumers = dict(
            stored_exp_resident_product=native.publish(leaf('div', leaf('mul', r, aa), xx)),
            resident_exp_product=native.publish(leaf('div', leaf('mul', exp(argument, False), aa), xx)),
            ratio_first=native.publish(leaf('mul', r, store(leaf('div', aa, xx)))),
            stored_argument=native.publish(leaf('div', store(leaf('mul', exp(store(argument)), aa)), xx)))
        density_t = store(leaf('mul', scalar(a-1.), load64(int(logs[x]['bits'], 16))))
        density_argument = leaf('sub', leaf('sub', density_t, xx),
                                load64(int(lg_density[a]['bits'], 16)))
        direct_pdf = native.publish(exp(density_argument))
        divided_prefix = store(leaf('div', exp(leaf('sub', t, xx)),
                                    load64(int(norm_gamma[a]['bits'],16))))
        cell = inverse(Span(binary64(density['bits']), binary64(density['bits'])),
                       binary64(x), 'div', True, 'stored_pc64')
        # The intermediate r*a is stored before division by x.
        intermediate = BINARY64.coordinates(cell)
        r_domains = []
        for k in ([] if intermediate is None else range(intermediate[0], intermediate[1]+1)):
            y = BINARY64.value(k)
            domain = inverse(Span(y, y), binary64(bits(a)), 'mul', True, 'stored_pc64')
            coordinates = BINARY64.coordinates(domain)
            if coordinates:
                r_domains.append(coordinates)
        r_coordinates = sorted({k for lo, hi in r_domains for k in range(lo, hi+1)})
        if len(r_coordinates) > 8:
            raise ValueError('unexpectedly broad prefactor preimage')
        prefixes = [BINARY64.value(k) for k in r_coordinates]
        target = Span(binary64(cdf['bits']), binary64(cdf['bits']))
        sum_domains = [inverse(target, p, 'mul', True, 'stored_pc64') for p in prefixes]
        predictions, compatible, iterations = {}, {}, {}
        for model in models:
            summ, n = series(a, x, model)
            predictions[model] = f'{native.publish(leaf("mul", r, summ)):016x}'
            compatible[model] = any(domain.contains(decode(summ)) for domain in sum_domains)
            iterations[model] = n
        strict_sum, _ = series(a, x, 'strict')
        for label, prefactor in (('resident_exp', exp(argument, False)),
                                 ('stored_argument', exp(store(argument)))):
            predictions[label] = f'{native.publish(leaf("mul", prefactor, strict_sum)):016x}'
        predictions['divided_prefix'] = f'{native.publish(leaf("mul", divided_prefix, strict_sum)):016x}'
        row = dict(shape=bits(a), input=x, density=density, cdf=cdf,
                   public_LN=logs[x], public_GAMMALN=lg[a],
                   prefactor=f'{native.publish(r):016x}', pdf_prediction=f'{predicted_pdf:016x}',
                   direct_pdf=f'{direct_pdf:016x}', public_density_GAMMALN=lg_density[a],
                   divided_prefix=f'{native.publish(divided_prefix):016x}', public_GAMMA=norm_gamma[a],
                   density_consumers={k:f'{v:016x}' for k,v in density_consumers.items()},
                   permitted_prefactors=[bits(float(p)) for p in prefixes],
                   cdf_predictions=predictions, series_possible=compatible, iterations=iterations)
        tally = by_shape.setdefault(str(a), Counter())
        for counter in (totals, tally):
            counter['pairs'] += 1
            counter['prefactor_pdf_match'] += predicted_pdf == int(density['bits'], 16)
            counter['direct_pdf_match'] += direct_pdf == int(density['bits'], 16)
            for model, result in density_consumers.items():
                counter[model+'_pdf_match'] += result == int(density['bits'], 16)
            counter['empty_prefactor_preimage'] += not prefixes
            for model in models:
                counter[model+'_possible'] += compatible[model]
                counter[model+'_cdf_match'] += predictions[model] == cdf['bits']
            for model in ('resident_exp', 'stored_argument', 'divided_prefix'):
                counter[model+'_cdf_match'] += predictions[model] == cdf['bits']
        records.append(row)
    result = dict(format='dnacalc-gamma-paired-stages/1', complete=True,
        selection='Fixed shape grid; adjacent windows at dyadic edges and ordinary decimal sites, restricted to x<a<20.',
        counts=totals, by_shape=by_shape, records=records, receipts=service.receipts,
        source_bindings={str(p):sha(p) for p in (Path(__file__), Path(library))},
        limits='The density inverse conditions on stored (r*a)/x and a shared prefactor. Empty constraints refute that family; nonempty constraints do not identify a series. Public GAMMALN and LN are independently observed, never fitted. Four ranked textbook series graphs differ at the recurrence only; no whole-program enumeration or count-selected changes.')
    write(output, result)
    print(json.dumps(dict(counts=totals, by_shape=by_shape), indent=2), flush=True)
    return result


def power_prefactor(checkout, output, cache, library, offline=False):
    """Remove the density's final scaling exactly, then invert one EXP input.

    Shape and x are powers of two, so every proposed a/x consumer is exact.
    The resulting r is an observed binary64, not a fitted prefactor. A common
    c in EXP(stored(a*LN(x))-x-c) must satisfy every x at the same shape.
    """
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    shapes = (2., 4., 8., 16.)
    cases = [(a, 2.**e) for a in shapes for e in range(-12, int(math.log2(a)))]
    questions = [[client.number(x), client.number(a), client.number(1.), client.boolean(c)]
                 for a, x in cases for c in (False, True)]
    answers = service.ask('GAMMA.DIST', questions)['rows']
    xs = sorted({x for _, x in cases})
    logs = dict(zip(xs, service.ask('LN', [[client.number(x)] for x in xs])['rows']))
    lg = dict(zip(shapes, service.ask('GAMMALN', [[client.number(a+1)] for a in shapes])['rows']))
    g = Builder(1)
    t = g.op('mul', 0, g.op('fldl2e', store=False), store=False)
    whole = g.op('frndint', t, store=False)
    small = g.op('f2xm1', g.op('sub', t, whole, store=False), store=False)
    positive = g.op('fscale', g.op('add', g.one, small, store=False), whole, store=False)
    out = g.op('div', g.one, positive, store=False)
    graph = g.graph(out, out)
    records, groups, totals = [], {}, Counter()
    log_groups = {}
    ln2 = native.leaf('fldln2', bytes(10), bytes(10), 0x133f)
    for i, (a, x) in enumerate(cases):
        density, cdf = answers[2*i:2*i+2]
        if density['type'] != 'number' or cdf['type'] != 'number':
            raise ValueError('finite power-of-two density/CDF must be numeric')
        # Exact rescaling, independently checked as a rational.
        prefix = value(density['bits'])*x/a
        if binary64(bits(prefix)) != binary64(density['bits'])*binary64(bits(x))/binary64(bits(a)):
            raise ValueError('density prefactor extraction rounded')
        argument = positive_extended_cell(native, graph, bits(prefix),
                                          binary64(bits(-math.log(prefix))), decreasing=True)
        if argument is None:
            raise ValueError('empty resident EXP cell requires a separate proof')
        m = Span(*(decode(bytes.fromhex(argument[k])) for k in ('low', 'high')))
        t = native.leaf('mul', load64(int(bits(a), 16)), load64(int(logs[x]['bits'], 16)), 0x133f)
        t = load64(native.publish(t))
        base = decode(native.leaf('sub', t, load64(int(bits(x), 16)), 0x133f))
        domain = EXTENDED64.pullback(m).linear(1, base)
        group = groups.setdefault(str(a), dict(domain=None, contradicted=False, rows=[]))
        group['rows'].append(i)
        if not group['contradicted']:
            group['domain'] = domain if group['domain'] is None else group['domain'].meet(domain)
            if group['domain'] is None:
                group['contradicted'] = True
                group['first_contradiction'] = i
        c = binary64(lg[a]['bits'])
        base_domain = EXTENDED64.pullback(m).linear(-1, c)
        t_domain = EXTENDED64.pullback(base_domain).linear(1, binary64(bits(x)))
        log_domain = inverse(t_domain, binary64(bits(a)), 'mul', True, 'stored_pc64')
        log_group = log_groups.setdefault(bits(x), dict(domain=None, contradicted=False, rows=[]))
        log_group['rows'].append(i)
        if not log_group['contradicted']:
            log_group['domain'] = log_domain if log_group['domain'] is None else log_group['domain'].meet(log_domain)
            if log_group['domain'] is None:
                log_group['contradicted'] = True
                log_group['first_contradiction'] = i
        resident_log = native.leaf('fyl2x', ln2, load64(int(bits(x),16)), 0x133f)
        resident_t = native.leaf('mul',load64(int(bits(a),16)),resident_log,0x133f)
        resident_base = decode(native.leaf('sub',resident_t,load64(int(bits(x),16)),0x133f))
        resident_domain = EXTENDED64.pullback(m).linear(1,resident_base)
        rr, term, summ = a, 1., 1.
        while True:
            rr += 1.
            term *= x/rr
            summ += term
            if term/summ <= 2.**-53:
                break
        published = native.publish(native.leaf('mul', load64(int(bits(prefix), 16)),
                                              load64(int(bits(summ), 16)), 0x133f))
        row = dict(shape=bits(a), input=bits(x), density=density, cdf=cdf,
                   extracted_prefactor=bits(prefix), public_LN=logs[x], public_GAMMALN=lg[a],
                   exponential_magnitude_cell=argument, normalizer_domain=domain.json(),
                   public_normalizer_possible=domain.contains(binary64(lg[a]['bits'])),
                   public_normalizer_with_resident_log_product=resident_domain.contains(c),
                   resident_log_product_normalizer_domain=resident_domain.json(),
                   log_domain=log_domain.json() if log_domain else None,
                   cdf_with_extracted_prefactor=f'{published:016x}')
        totals['rows'] += 1
        totals['public_normalizer_possible'] += row['public_normalizer_possible']
        totals['public_normalizer_with_resident_log_product'] += row['public_normalizer_with_resident_log_product']
        totals['cdf_with_extracted_prefactor_match'] += published == int(cdf['bits'], 16)
        records.append(row)
    for group in list(groups.values())+list(log_groups.values()):
        domain = group.pop('domain')
        coords = BINARY64.coordinates(domain) if domain is not None else None
        group.update(domain=domain.json() if domain else None,
                     permitted_binary64_count=0 if coords is None else coords[1]-coords[0]+1,
                     binary64_values=[] if coords is None else [bits(float(BINARY64.value(k)))
                         for k in range(coords[0], min(coords[1]+1, coords[0]+8))])
    result = dict(format='dnacalc-gamma-power-prefactor/1', complete=True, counts=totals,
        groups=groups, log_groups=log_groups, records=records, graph=graph, receipts=service.receipts,
        source_bindings={str(p):sha(p) for p in (Path(__file__), Path(library))},
        limits='Power-of-two density scaling is exact. The extracted r is a conditional observed prefactor; sharing it with the CDF is separately tested. Normalizer domains condition on one stored public LN/product and the declared resident negative EXP chain, with PC64 subtraction. A nonempty common interval is a necessary constraint, never an identified coefficient or a fitted implementation.')
    write(output, result)
    print(json.dumps(dict(counts=totals, groups=groups), indent=2), flush=True)
    return result


def power_primitives(study, output, cache, library):
    """Invert one exponential primitive with the density scale removed exactly.

    The fixed textbook exponentials vary negative handling and constant
    precision. The unknown normalizer must be common to all x at a
    fixed shape. No constant inferred here is installed in an implementation.
    """
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete power-prefactor observations')
    answers, bindings, ordered = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    native = Native(library)
    graphs = {}
    # Code-shape rank is fixed before the outputs: the ordinary signed
    # textbook sequence, then the public-EXP reciprocal sequence.
    for kind in ('signed', 'reciprocal', 'signed_constant53'):
        g = Builder(1)
        exponent = 0 if kind == 'reciprocal' else g.op('neg', 0, store=False)
        constant = (g.constant('3ff71547652b82fe') if kind == 'signed_constant53'
                    else g.op('fldl2e', store=False))
        t = g.op('mul', exponent, constant, store=False)
        whole = g.op('frndint', t, store=False)
        fraction = g.op('sub', t, whole, store=False)
        small = g.op('f2xm1', fraction, store=False)
        result = g.op('fscale', g.op('add', g.one, small, store=False), whole, store=False)
        if kind == 'reciprocal':
            result = g.op('div', g.one, result, store=False)
        graphs[kind] = g.graph(result, result)
    groups, records = {}, []
    for row in doc['records']:
        a, x = map(binary64, (row['shape'], row['input']))
        density = row['density']
        original_args = [dict(type='number', bits=b) for b in
                         (row['input'], row['shape'], bits(1.))] + [dict(type='boolean', value=False, mode='cell')]
        if answers[('GAMMA.DIST', json.dumps(original_args, sort_keys=True))] != density:
            raise ValueError('power density observation changed')
        prefix = binary64(row['extracted_prefactor'])
        if prefix != binary64(density['bits'])*x/a:
            raise ValueError('power-of-two density scaling was not exact')
        t = native.leaf('mul', load64(int(row['shape'], 16)),
                        load64(int(row['public_LN']['bits'], 16)), 0x133f)
        t = load64(native.publish(t))
        base = decode(native.leaf('sub', t, load64(int(row['input'], 16)), 0x133f))
        predictions = {}
        for kind, graph in graphs.items():
            cell = positive_extended_cell(native, graph, row['extracted_prefactor'],
                                          binary64(bits(-math.log(float(prefix)))), decreasing=True)
            group = groups.setdefault((kind, row['shape']), dict(domain=None, contradicted=False,
                                                                 rows=0, public_possible=0))
            group['rows'] += 1
            if cell is None:
                domain = None
            else:
                magnitudes = Span(*(decode(bytes.fromhex(cell[k])) for k in ('low', 'high')))
                domain = EXTENDED64.pullback(magnitudes).linear(1, base)
            if not group['contradicted']:
                group['domain'] = domain if group['domain'] is None else (
                    None if domain is None else group['domain'].meet(domain))
                if group['domain'] is None:
                    group['contradicted'] = True
                    group['first_contradiction'] = len(records)
            c = binary64(row['public_GAMMALN']['bits'])
            possible = domain is not None and domain.contains(c)
            group['public_possible'] += possible
            predictions[kind] = dict(exponential_cell=cell,
                                     normalizer_domain=domain.json() if domain else None,
                                     public_normalizer_possible=possible)
        records.append(dict(shape=row['shape'], input=row['input'],
                            extracted_prefactor=row['extracted_prefactor'], primitives=predictions))
    reported = []
    for (kind, shape), group in groups.items():
        domain = group.pop('domain')
        coords = BINARY64.coordinates(domain) if domain else None
        reported.append(dict(primitive=kind, shape=shape, **group,
                             common_normalizer=domain.json() if domain else None,
                             binary64_count=0 if coords is None else coords[1]-coords[0]+1,
                             binary64_values=[] if coords is None else [bits(float(BINARY64.value(k)))
                                 for k in range(coords[0], min(coords[1]+1, coords[0]+8))]))
    result = dict(format='dnacalc-gamma-power-primitive-inverses/1', complete=True,
                  ranked_primitives=list(graphs), graphs=graphs, groups=reported, records=records,
                  original_calls=len(ordered), receipts=doc['receipts'],
                  source_bindings=frozen,
                  original_bindings=bindings,
                  limits='Original power-of-two density scaling removes the consumer exactly. '
                  'Each inverse tests one fixed native exponential with a shared but unknown normalizer. '
                  'Empty common constraints refute that primitive/producer composition; nonempty '
                  'constraints do not identify a coefficient or a private logarithm. No kernel is '
                  'chosen from match counts and no inferred constant is installed.')
    write(output, result)
    print(json.dumps(reported, indent=2), flush=True)
    return result


def shared_log_constraints(study, output, cache, library):
    """Test a shared unknown log and normalizer without fitting either.

    Powers of two make a*L exact for a binary64 L. For this input grid,
    a*L-x is also exact in PC64 whenever .5 <= abs(L) <= 16 (L(1)=0).
    Therefore each inverse gives one necessary difference constraint on
    L(x)-C(a)/a. A negative cycle refutes the whole declared family even
    after relaxing all interval endpoints to closed bounds.
    """
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete primitive inverses')
    _, bindings, _ = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    from fractions import Fraction as Q
    results = []
    for primitive in doc['ranked_primitives']:
        edges = []
        nodes = set()
        for index, row in enumerate(doc['records']):
            a, x = binary64(row['shape']), binary64(row['input'])
            cell = row['primitives'][primitive]['exponential_cell']
            if cell is None:
                raise ValueError('empty primitive cell needs its own certificate')
            magnitudes = Span(*(decode(bytes.fromhex(cell[k])) for k in ('low', 'high')))
            domain = EXTENDED64.pullback(magnitudes).linear(-1/a, x/a)
            log, normalizer = 'L:'+row['input'], 'C/a:'+row['shape']
            nodes.update((log, normalizer))
            edges.extend([dict(source=normalizer, target=log, bound=domain.high, row=index),
                          dict(source=log, target=normalizer, bound=-domain.low, row=index)])
        distance = {name:Q(0) for name in nodes}
        predecessor = {}
        changed = None
        for _ in range(len(nodes)):
            changed = None
            for i, edge in enumerate(edges):
                candidate = distance[edge['source']]+edge['bound']
                if candidate < distance[edge['target']]:
                    distance[edge['target']] = candidate
                    predecessor[edge['target']] = i
                    changed = edge['target']
            if changed is None:
                break
        cycle = []
        if changed is not None:
            node = changed
            for _ in nodes:
                node = edges[predecessor[node]]['source']
            start = node
            while True:
                edge = edges[predecessor[node]]
                cycle.append(edge)
                node = edge['source']
                if node == start:
                    break
                if len(cycle) > len(nodes):
                    raise ValueError('invalid negative cycle reconstruction')
            total = sum((e['bound'] for e in cycle), Q(0))
            if total >= 0:
                raise ValueError('cycle does not prove a contradiction')
        else:
            total = Q(0)
            if any(distance[e['target']] > distance[e['source']]+e['bound'] for e in edges):
                raise ValueError('feasible witness violates a necessary constraint')
        serial = lambda q:[str(q.numerator), str(q.denominator)]
        results.append(dict(primitive=primitive, continuous_constraints_feasible=changed is None,
                            negative_cycle=[dict(e, bound=serial(e['bound'])) for e in cycle],
                            cycle_bound_sum=serial(total),
                            edges=[dict(e, bound=serial(e['bound'])) for e in edges],
                            feasible_potentials=None if changed is not None else
                                {k:serial(v) for k,v in distance.items()}))
    result = dict(format='dnacalc-gamma-shared-log-constraints/1', complete=True,
                  results=results, receipts=doc['receipts'], source_bindings=frozen,
                  limits='Necessary constraints only: one stored binary64 private L(x), '
                  '.5<=abs(L)<=16 away from x=1, L(1)=0, one normalizer C(a), '
                  'stored a*L, resident differences and one declared exponential, with '
                  'the exact a/x density consumer. Continuous variables and closed '
                  'endpoints relax the model. A negative cycle is a refutation even '
                  'of that relaxation. A feasible witness identifies no private value '
                  'and supplies no implementation or count-selected correction.')
    write(output, result)
    print(json.dumps([dict(primitive=r['primitive'], feasible=r['continuous_constraints_feasible'],
                           cycle_rows=[e['row'] for e in r['negative_cycle']],
                           cycle_bound_sum=r['cycle_bound_sum']) for r in results]), flush=True)
    return result


def split_prefactor(checkout, study, output, cache, library, offline=False):
    """Test the literal power-times-exponential normalizer on dyadic inputs."""
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete power-prefactor study')
    _, bindings, _ = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    xs = sorted({r['input'] for r in doc['records']})
    shapes = sorted({r['shape'] for r in doc['records']})
    # Independent public leaves; original density/CDF answers are reused.
    exps = dict(zip(xs, service.ask('EXP', [[client.number(bits(-value(x)))] for x in xs])['rows']))
    gammas = dict(zip(shapes, service.ask('GAMMA', [[client.number(value(a)+1.)] for a in shapes])['rows']))
    groups = {}
    records = []
    for row in doc['records']:
        a, x = binary64(row['shape']), binary64(row['input'])
        if a.denominator != 1:
            raise ValueError('power experiment needs an integer shape')
        e, gamma = exps[row['input']], gammas[row['shape']]
        if e['type'] != 'number' or gamma['type'] != 'number':
            raise ValueError('positive finite leaf must be numeric')
        product = x**a.numerator*binary64(e['bits'])
        # On this grid the power and its multiplication by EXP(-x) are exact.
        if BINARY64.round(product) != product:
            raise ValueError('dyadic numerator scaling rounded')
        target = binary64(row['extracted_prefactor'])
        domain = inverse(Span(target, target), product, 'div', False, 'stored_pc64')
        group = groups.setdefault(row['shape'], dict(domain=None, contradicted=False, rows=0,
                                                     public_gamma_possible=0))
        group['rows'] += 1
        if not group['contradicted']:
            group['domain'] = domain if group['domain'] is None else (
                None if domain is None else group['domain'].meet(domain))
            if group['domain'] is None:
                group['contradicted'] = True
                group['first_contradiction'] = len(records)
        possible = domain is not None and domain.contains(binary64(gamma['bits']))
        group['public_gamma_possible'] += possible
        records.append(dict(shape=row['shape'], input=row['input'], public_exp=e,
                            public_gamma=gamma, density=row['density'],
                            extracted_prefactor=row['extracted_prefactor'],
                            normalizer_domain=domain.json() if domain else None,
                            public_gamma_possible=possible))
    reported = []
    for shape, group in groups.items():
        domain = group.pop('domain')
        coords = BINARY64.coordinates(domain) if domain else None
        reported.append(dict(shape=shape, **group, common_normalizer=domain.json() if domain else None,
                             binary64_count=0 if coords is None else coords[1]-coords[0]+1,
                             binary64_values=[] if coords is None else [bits(float(BINARY64.value(k)))
                                 for k in range(coords[0], min(coords[1]+1, coords[0]+8))]))
    result = dict(format='dnacalc-gamma-split-prefactor/1', complete=True,
                  groups=reported, records=records, receipts=doc['receipts']+service.receipts,
                  source_bindings=frozen,
                  limits='Exact powers of two remove the power and numerator-product rounding. '
                  'Original density outputs and new independent public EXP/GAMMA leaves constrain '
                  'one unknown shared Gamma(a+1) divisor. A common interval is necessary only; '
                  'no normalizer or algorithm is chosen from match counts or installed.')
    write(output, result)
    print(json.dumps(reported, indent=2), flush=True)
    return result


def poisson_pairs(checkout, study, output, cache, library, offline=False):
    """Check a possible public Poisson dependency before changing gamma helpers."""
    from geomean import verified_study_answers
    source = Path(study)
    doc = json.loads(source.read_text())
    if doc.get('complete') is not True:
        raise ValueError('incomplete gamma observations')
    _, bindings, _ = verified_study_answers(doc, Path(cache))
    frozen = {**STAGE_BINDINGS, **bindings,
              str(source.resolve()):sha(source), str(Path(library).resolve()):sha(library)}
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    questions = [[client.number(value(r['shape'])-1.), client.number(r['input']), client.boolean(flag)]
                 for r in doc['records'] for flag in (False, True)]
    answers = service.ask('POISSON.DIST', questions)['rows']
    counts, records = Counter(), []
    for index, row in enumerate(doc['records']):
        pmf, cdf = answers[2*index:2*index+2]
        same = pmf == row['density']
        complement = cdf['type'] == 'number' and row['cdf'] == dict(type='number', bits=bits(1.-value(cdf['bits'])))
        counts['pairs'] += 1
        counts['pdf_equals_public_poisson'] += same
        counts['cdf_equals_public_poisson_complement'] += complement
        records.append(dict(shape=row['shape'], input=row['input'], gamma_density=row['density'],
                            gamma_cdf=row['cdf'], poisson_pmf=pmf, poisson_cdf=cdf,
                            density_equal=same, cdf_complement_equal=complement))
    result = dict(format='dnacalc-gamma-poisson-public-pairs/1', complete=True, counts=counts,
                  records=records, receipts=doc['receipts']+service.receipts, source_bindings=frozen,
                  limits='Original gamma answers, independently measured Poisson calls at integer '
                  'shape minus one. An equality is a functional observation, not proof of a '
                  'private call or of the noninteger route. Complement roundoff is explicit; '
                  'no correction, selector or implementation is inferred from these counts.')
    write(output, result)
    print(json.dumps(counts), flush=True)
    return result


def gammaln_store_edges(checkout, output, cache, library, offline=False, count=1024):
    """Observe the mean-error jump at each final polynomial's binade.

    A final binary64 sum rounds the constant coefficient on a coarser
    lattice when the polynomial doubles. An unstored final PC64 sum has
    eleven more bits there. Edges come from the existing coefficient table,
    before any new Excel answer. Coefficients are never changed in this study.
    """
    import re
    import struct
    try:
        import mpmath as mp
    except ImportError:
        sys.path.insert(0, '/home/ubuntu/claude-work/gamma/pylib')
        import mpmath as mp
    mp.mp.prec = 200
    path = Path(checkout)/'crates/calc-func/src/ported/excel_numeric/gammaln.rs'
    frozen = {**STAGE_BINDINGS, str(path.resolve()):sha(path),
              str(Path(library).resolve()):sha(library), str(Path(mp.__file__).resolve()):sha(mp.__file__)}
    text = path.read_text()
    tables = {}
    for name in ('P', 'Q'):
        block = text.split('const B3_'+name+':', 1)[1].split('];', 1)[0]
        tables[name] = [v.replace('_', '').lower() for v in re.findall(r'0x([0-9A-Fa-f_]+)', block)]
        if len(tables[name]) != 9:
            raise ValueError('changed B3 coefficient table shape')
    def polynomial(name, x):
        result = mp.mpf(0)
        for b in reversed(tables[name]):
            result = result*x+mp.mpf(value(b))
        return result
    edges = []
    for name, level in (('P', 4096), ('P', 8192), ('Q', 8192), ('Q', 16384)):
        low, high = mp.mpf(2), mp.mpf(3)
        for _ in range(180):
            middle = (low+high)/2
            if polynomial(name, middle) < level:
                low = middle
            else:
                high = middle
        center = float((low+high)/2)
        windows = []
        for delta in (2.**-16, 2.**-20):
            for sign in (-1, 1):
                first = int(bits(center+sign*delta), 16)-count//2
                windows.append(dict(center=bits(center+sign*delta), delta=bits(delta), sign=sign,
                                    inputs=[f'{first+k:016x}' for k in range(count)]))
        edges.append(dict(polynomial=name, level=level, estimated_edge=bits(center), windows=windows))
    graphs = {}
    for numerator_stored, denominator_stored in ((False,False), (True,False), (False,True)):
        g = Builder(1)
        def op(kind, a, b, stored=True):
            index = g.op(kind, a, b, store=False)
            if stored:
                g.nodes[index]['site'] = 0
                g.sites = 1
            return index
        def chain(name, stored):
            v = g.constant(tables[name][-1])
            for i in range(7, -1, -1):
                v = op('mul', v, 0)
                v = op('add', v, g.constant(tables[name][i]), stored=i != 0 or stored)
            return v
        p, q = chain('P', numerator_stored), chain('Q', denominator_stored)
        ratio = op('div', p, q)
        difference = op('sub', 0, g.constant(bits(2.)), stored=False)
        result = op('mul', difference, ratio)
        label = 'P'+('stored' if numerator_stored else 'resident')+'_Q'+('stored' if denominator_stored else 'resident')
        graphs[label] = g.graph(result, result)
    manifest = dict(format='dnacalc-gammaln-store-edge-questions/1', complete=True,
                    tables=tables, edges=edges, graphs=graphs, source_bindings=frozen,
                    selection='Existing polynomial crossings of four fixed powers of two; '
                    'two fixed symmetric distances and consecutive-double windows. No Excel output selects a site.')
    write(str(output)+'.questions.json', manifest)
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    records = []
    for edge in edges:
        for window in edge['windows']:
            xs = window['inputs']
            answers = service.ask('GAMMALN', [[client.number(x)] for x in xs])['rows']
            center = mp.mpf(value(window['center']))
            f0, f1, f2 = mp.loggamma(center), mp.digamma(center), mp.polygamma(1, center)
            sums = dict(Excel=mp.mpf(0), **{k:mp.mpf(0) for k in graphs})
            squares = dict(Excel=mp.mpf(0), **{k:mp.mpf(0) for k in graphs})
            for x, answer in zip(xs, answers):
                if answer['type'] != 'number':
                    raise ValueError('GAMMALN edge answer must be numeric')
                xx = mp.mpf(value(x)); dx = xx-center
                truth = f0+f1*dx+f2*dx*dx/2
                unit = mp.mpf(math.ulp(value(answer['bits'])))
                predictions = {k:native.publish(native.evaluate(graph, [int(x,16)], 0x133f, mask=1)[:10])
                               for k,graph in graphs.items()}
                observed = dict(Excel=int(answer['bits'],16), **predictions)
                for name, raw in observed.items():
                    residual = (mp.mpf(struct.unpack('>d',struct.pack('>Q',raw))[0])-truth)/unit
                    sums[name] += residual
                    squares[name] += residual*residual
            n = len(xs)
            means = {k:float(v/n) for k,v in sums.items()}
            sem = {k:float(mp.sqrt(max(mp.mpf(0), squares[k]/n-(sums[k]/n)**2)/n)) for k in sums}
            records.append(dict(polynomial=edge['polynomial'], level=edge['level'],
                                center=window['center'], delta=window['delta'], sign=window['sign'],
                                first_input=xs[0], last_input=xs[-1], calls=n,
                                mean_error_in_result_ulps=means, descriptive_standard_error=sem))
    result = dict(format='dnacalc-gammaln-final-sum-store-edges/1', complete=True,
                  records=records, receipts=service.receipts, source_bindings=frozen,
                  questions_sha256=sha(str(output)+'.questions.json'),
                  limits='Stage-binade observations only. Means use a 200-bit mathematical '
                  'log-gamma and a quadratic expansion in each consecutive window; '
                  'the third-derivative remainder is below 1e-20 result ULP on this grid. '
                  'Standard errors describe within-window residual variation, not independent '
                  'random sampling. Existing coefficients remain approximate and unchanged; '
                  'no table, store or threshold is chosen by exact-match counts.')
    write(output, result)
    print(json.dumps(records, indent=2), flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout')
    p.add_argument('output')
    p.add_argument('--cache', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--count', type=int, default=32)
    p.add_argument('--offline', action='store_true')
    p.add_argument('--power-prefactor', action='store_true')
    p.add_argument('--power-primitives')
    p.add_argument('--shared-log-constraints')
    p.add_argument('--split-prefactor')
    p.add_argument('--poisson-pairs')
    p.add_argument('--gammaln-store-edges', action='store_true')
    args = p.parse_args()
    if not 2 <= args.count <= (4096 if args.gammaln_store_edges else 256):
        raise ValueError('small adjacent windows only')
    if args.gammaln_store_edges:
        gammaln_store_edges(args.checkout, args.output, args.cache, args.library, args.offline, args.count)
    elif args.poisson_pairs:
        poisson_pairs(args.checkout, args.poisson_pairs, args.output, args.cache, args.library, args.offline)
    elif args.split_prefactor:
        split_prefactor(args.checkout, args.split_prefactor, args.output, args.cache, args.library, args.offline)
    elif args.shared_log_constraints:
        shared_log_constraints(args.shared_log_constraints, args.output, args.cache, args.library)
    elif args.power_primitives:
        power_primitives(args.power_primitives, args.output, args.cache, args.library)
    elif args.power_prefactor:
        power_prefactor(args.checkout, args.output, args.cache, args.library, args.offline)
    else:
        paired(args.checkout, args.output, args.cache, args.library, args.count, args.offline)
