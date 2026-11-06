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


def power_primitives(study, output, cache, library, direct_density=False):
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
        target = density['bits'] if direct_density else row['extracted_prefactor']
        prefix = binary64(target)
        if not direct_density and prefix != binary64(density['bits'])*x/a:
            raise ValueError('power-of-two density scaling was not exact')
        coefficient=bits(float(a-1)) if direct_density else row['shape']
        t = native.leaf('mul', load64(int(coefficient, 16)),
                        load64(int(row['public_LN']['bits'], 16)), 0x133f)
        t = load64(native.publish(t))
        base = decode(native.leaf('sub', t, load64(int(row['input'], 16)), 0x133f))
        predictions = {}
        for kind, graph in graphs.items():
            cell = positive_extended_cell(native, graph, target,
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
            possible = None if direct_density else domain is not None and domain.contains(c)
            group['public_possible'] += possible is True
            predictions[kind] = dict(exponential_cell=cell,
                                     normalizer_domain=domain.json() if domain else None,
                                     public_normalizer_possible=possible)
        records.append(dict(shape=row['shape'], input=row['input'],
                            density=density, exponential_target=target,
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
    result = dict(format=('dnacalc-gamma-density-primitive-inverses/1' if direct_density else
                          'dnacalc-gamma-power-primitive-inverses/1'), complete=True,
                  exponent_coefficient='a-1' if direct_density else 'a',
                  ranked_primitives=list(graphs), graphs=graphs, groups=reported, records=records,
                  original_calls=len(ordered), receipts=doc['receipts'],
                  source_bindings=frozen,
                  original_bindings=bindings,
                  limits=('The observed density directly targets EXP(stored((a-1)*public_LN(x))-x-C(a)). '
                          'The shape-a public GAMMALN is not observed here; no normalizer value is compared. '
                          if direct_density else 'Original power-of-two density scaling removes the consumer exactly. ')+
                  'Each inverse tests one fixed native exponential with a shared but unknown normalizer. '
                  'Empty common constraints refute that primitive/producer composition; nonempty '
                  'constraints do not identify a coefficient or a private logarithm. No kernel is '
                  'chosen from match counts and no inferred constant is installed.')
    write(output, result)
    print(json.dumps(reported, indent=2), flush=True)
    return result


def shared_log_constraints(study, output, cache, library, scaled_log=False):
    """Test a shared unknown log and normalizer without fitting either.

    Powers of two make a*L exact for a binary64 L. For this input grid,
    a*L-x is also exact in PC64 whenever .5 <= abs(L) <= 16 (L(1)=0).
    Therefore each inverse gives one necessary difference constraint on
    L(z)-C(a)/a, with z=x or the explicitly selected exact x/a. A negative cycle refutes the whole declared family even
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
            log_input=bits(float(x/a)) if scaled_log else row['input']
            if scaled_log and binary64(log_input)!=x/a:raise ValueError('scaled log argument rounded')
            log, normalizer = 'L:'+log_input, 'C/a:'+row['shape']
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
                  logarithm_argument='x/a' if scaled_log else 'x',
                  results=results, receipts=doc['receipts'], source_bindings=frozen,
                  limits='Necessary constraints only: one stored binary64 private L at the stated argument, '
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


def unit_shape_exponential(checkout,bank,manifest,output,cache,library,offline=False):
    """Use shape-one density to expose an independent exponential wrapper.

    Existing original density inputs choose centers only. Their answers are
    retained, never re-asked; both successor inputs get paired new calls.
    """
    import gzip
    source=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz'
    old=[]
    for line in gzip.open(source,'rt'):
        raw=json.loads(line);args=raw['args']
        if len(args)!=1 or len(args[0])!=4 or any(a['kind']!='number' for a in args[0]):continue
        x,a,b,c=(v['value'] for v in args[0])
        if a==b==1. and c==0. and 0.<x<700.:old.append(raw)
    if not old:raise ValueError('no original unit-shape densities')
    primitive=Path(manifest);doc=json.loads(primitive.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete exponential graphs')
    graphs=doc['graphs'];known={r['args'][0][0]['bits_hex'].removeprefix('0x') for r in old}
    inputs=sorted({f'{int(x,16)+offset:016x}' for x in known for offset in (1,2)}-known)
    if len(inputs)>1500:raise ValueError('bounded adjacent-input study only')
    native=Native(library);predictions={x:{name:f'{native.publish(native.evaluate(graph,[int(x,16)],0x133f)[:10]):016x}'
                                        for name,graph in graphs.items()} for x in inputs}
    frozen={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(primitive.resolve()):sha(primitive),
            str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-unit-shape-questions/1',complete=True,
        original_centers=old,inputs=inputs,predictions=predictions,graphs=graphs,source_bindings=frozen,
        selection='First two successors of every original positive finite unit-scale shape-one density below700. All original density arguments are excluded before new calls; no current match result selects inputs.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    densities=service.ask('GAMMA.DIST',[[client.number(x),client.number(1.),client.number(1.),client.number(0.)] for x in inputs])['rows']
    exps=service.ask('EXP',[[client.number(bits(-value(x)))] for x in inputs])['rows']
    counts=Counter();records=[]
    for x,d,e in zip(inputs,densities,exps):
        if d['type']!='number' or e['type']!='number':raise ValueError('nonnumeric bounded density/EXP')
        counts['paired_rows']+=1;counts['same_public_EXP']+=d==e
        for name,prediction in predictions[x].items():
            counts['density/'+name]+=d==dict(type='number',bits=prediction)
            counts['EXP/'+name]+=e==dict(type='number',bits=prediction)
        records.append(dict(input=x,density=d,public_EXP=e,predictions=predictions[x]))
    result=dict(format='dnacalc-gamma-unit-shape-exponential/1',complete=True,counts=counts,records=records,
        original_centers=old,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Paired original current-build calls at new adjacent inputs. Shape1/unit scale '
        'mathematically isolates exp(-x); this does not assume other shapes call the same '
        'private exponential or identify a resident-argument ABI. The three fixed native '
        'graphs are comparisons, never selected by their number of matches. Historical '
        'center answers retain their separate environment and are not new coverage.')
    write(output,result);print(json.dumps(counts),flush=True);return result


def unit_shape_prefactor(checkout,study,output,cache,library,offline=False):
    """Test the literal log/prefactor/divide route without cancelling x."""
    from geomean import verified_study_answers
    import itertools
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete unit-shape pairs')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    questions=Path(str(source)+'.questions.json')
    if sha(questions)!=doc['questions_sha256']:raise ValueError('changed original question manifest')
    graphs=json.loads(questions.read_text())['graphs'];inputs=[r['input'] for r in doc['records']]
    profiles=[dict(logarithm=log,argument_stored=arg,prefix_stored=prefix,division=division,exponential=exp)
              for log,arg,prefix,division,exp in itertools.product(
                  ('public_LN','resident_LN'),(True,False),(True,False),('stored_pc64','pc53'),graphs)]
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-unit-prefactor-questions/1',complete=True,
        inputs=inputs,profiles=profiles,source_bindings=STAGE_BINDINGS,
        selection='EXP(LN(x)-x)/x, three fixed exponentials and one-at-a-time publications. '
        'No Gamma1=exp(-x) algebraic cancellation, inferred coefficient or fitted normalizer.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    logs=service.ask('LN',[[client.number(x)] for x in inputs])['rows']
    native=Native(library);ln2=native.leaf('fldln2',bytes(10),bytes(10),0x133f)
    counts=[Counter() for _ in profiles];records=[]
    for row,log in zip(doc['records'],logs):
        if log['type']!='number':raise ValueError('nonnumeric public logarithm')
        xx=load64(int(row['input'],16));ll=load64(int(log['bits'],16))
        resident=native.leaf('fyl2x',ln2,xx,0x133f);predictions=[]
        for index,p in enumerate(profiles):
            magnitude=native.leaf('sub',xx,ll if p['logarithm']=='public_LN' else resident,0x133f)
            if p['argument_stored']:magnitude=load64(native.publish(magnitude))
            prefix=native.evaluate(graphs[p['exponential']],[],0x133f,extended=magnitude)[:10]
            if p['prefix_stored']:prefix=load64(native.publish(prefix))
            raw=native.leaf('div',prefix,xx,0x133f if p['division']=='stored_pc64' else 0x123f)
            predicted=f'{native.publish(raw):016x}';predictions.append(predicted)
            counts[index]['matched' if row['density']==dict(type='number',bits=predicted) else 'mismatched']+=1
        records.append(dict(input=row['input'],density=row['density'],public_EXP=row['public_EXP'],public_LN=log,predictions=predictions))
    frozen={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source),str(questions.resolve()):sha(questions),
            str(Path(library).resolve()):sha(library)}
    result=dict(format='dnacalc-gamma-unit-prefactor-stores/1',complete=True,profiles=profiles,counts=counts,
        records=records,receipts=doc['receipts']+service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='All48 literal profiles are fixed before new independent LN leaves. Original '
        'density and EXP answers remain unchanged. A full agreement establishes only this '
        'sampled functional composition, not uniqueness, other shapes or an unmeasured '
        'normalizer. A failed profile is never replaced by its nearest answer or chosen by counts.')
    write(output,result);print(json.dumps(dict(profiles=len(profiles),complete_survivors=[(i,p) for i,(p,c) in enumerate(zip(profiles,counts)) if not c['mismatched']])));return result


def unit_shape_constants(study,output,cache,library):
    """Constrain one common normalization or exponent shift, without installing it."""
    from geomean import verified_study_answers
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete shape-one pairs')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    questions=Path(str(source)+'.questions.json')
    if sha(questions)!=doc['questions_sha256']:raise ValueError('changed original primitive questions')
    graphs=json.loads(questions.read_text())['graphs'];native=Native(library)
    groups={name:dict(domain=None,contradicted=False,rows=0) for name in
            ('divide_stored_pc64','multiply_stored_pc64',*(('shift/'+k) for k in graphs))}
    records=[]
    for index,row in enumerate(doc['records']):
        target=binary64(row['density']['bits']);e=binary64(row['public_EXP']['bits']);x=binary64(row['input'])
        domains={name:inverse(Span(target,target),e,op,left,'stored_pc64') for name,op,left in
                 [('divide_stored_pc64','div',False),('multiply_stored_pc64','mul',True)]}
        cells={}
        for kind,graph in graphs.items():
            cell=positive_extended_cell(native,graph,row['density']['bits'],binary64(bits(-math.log(float(target)))),decreasing=True,max_radius_bits=56)
            cells[kind]=cell
            if cell is None:domain=None
            else:
                magnitudes=Span(*(decode(bytes.fromhex(cell[k])) for k in ('low','high')))
                domain=EXTENDED64.pullback(magnitudes).linear(1,-x)
            domains['shift/'+kind]=domain
        for name,domain in domains.items():
            group=groups[name];group['rows']+=1
            if not group['contradicted']:
                group['domain']=domain if group['domain'] is None else (None if domain is None else group['domain'].meet(domain))
                if group['domain'] is None:group['contradicted']=True;group['first_contradiction']=index
        records.append(dict(input=row['input'],density=row['density'],public_EXP=row['public_EXP'],
            domains={k:v.json() if v is not None else None for k,v in domains.items()},exponential_cells=cells))
    reported={}
    for name,group in groups.items():
        domain=group.pop('domain');coords=BINARY64.coordinates(domain) if domain is not None else None
        reported[name]=dict(group,common_interval=domain.json() if domain is not None else None,
            binary64_count=0 if coords is None else coords[1]-coords[0]+1)
    result=dict(format='dnacalc-gamma-unit-shape-normalizer-constraints/1',complete=True,groups=reported,records=records,
        source_bindings={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source),str(questions.resolve()):sha(questions),str(Path(library).resolve()):sha(library)},
        receipts=doc['receipts'],
        limits='One common multiplicative/divisive normalization of the measured public EXP '
        'or one common additive shift before a declared native exponential. Intersections '
        'are necessary only; no constant is chosen or installed. An empty cell rules out '
        'that complete declared scalar normalization, not a log/power route or a private '
        'exponential. Shape-one density is mathematically exp(-x), but its implementation '
        'is not assumed to be an exponential-only wrapper.')
    write(output,result);print(json.dumps(reported),flush=True);return result


def unit_shape_consumers(study,output,cache):
    """Invert only the final density consumer, with an arbitrary prefix."""
    from geomean import verified_study_answers
    from exact import operation
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete original unit-shape study')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    profiles=[dict(kind=kind,model=model) for kind in ('divide_by_x','multiply_reciprocal_x')
              for model in ('strict','stored_pc64')]
    counts=[Counter() for _ in profiles];records=[]
    for row in doc['records']:
        x=binary64(row['input']);y=binary64(row['density']['bits']);constraints=[]
        for index,p in enumerate(profiles):
            if p['kind']=='divide_by_x':constant=x;op='div'
            else:constant=operation(1,x,'div',p['model']);op='mul'
            domain=inverse(Span(y,y),constant,op,True,p['model'])
            coords=BINARY64.coordinates(domain) if domain is not None else None
            resident=EXTENDED64.coordinates(domain) if domain is not None else None
            n=0 if coords is None else coords[1]-coords[0]+1
            tally=counts[index];tally['rows']+=1;tally['empty_binary64_prefix']+=n==0
            tally['empty_resident64_prefix']+=resident is None
            neighbors=[]
            if n==0 and domain is not None:
                ks=[BINARY64.floor_coordinate(domain.low),BINARY64.ceiling_coordinate(domain.high)]
                for k in ks:
                    prefix=BINARY64.value(k);result=operation(prefix,constant,op,p['model'])
                    neighbors.append(dict(prefix=bits(float(prefix)),calculated=bits(float(result))))
                if not binary64(neighbors[0]['calculated'])<y<binary64(neighbors[1]['calculated']):
                    raise ValueError('binary64 prefix gap not certified by its actual neighboring results')
            constraints.append(dict(profile=index,required_prefix=domain.json() if domain is not None else None,
                binary64_prefix_count=n,resident64_prefix_count=0 if resident is None else resident[1]-resident[0]+1,
                adjacent_binary64_outputs=neighbors))
        records.append(dict(input=row['input'],density=row['density'],constraints=constraints))
    result=dict(format='dnacalc-gamma-unit-final-consumer-preimages/1',complete=True,profiles=profiles,counts=counts,records=records,
        receipts=doc['receipts'],source_bindings={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source)},
        limits='Any positive prefix, without a log, exponential, coefficient or normalizer producer. '
        'A binary64 gap excludes that one final stored-prefix consumer; its exact neighboring '
        'calculated outputs certify the gap. A resident preimage is necessary only, not an '
        'ABI or private algorithm identification. No prefix is fitted or installed and '
        'these feasibility counts are never product pass counts.')
    write(output,result);print(json.dumps(list(zip(profiles,counts))),flush=True);return result


def scaled_density_consumers(checkout,study,bank,output,cache,offline=False):
    """Hold shape and exact x/beta fixed while changing the final scale."""
    from geomean import verified_study_answers
    from exact import operation
    import gzip
    source=Path(study);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete power-density anchors')
    _,bindings,_=verified_study_answers(doc,Path(cache))
    profiles=[dict(kind=kind,model=model,prefix_role='post_alpha_product' if kind=='divide' else 'unscaled_prefactor')
              for kind in ('divide','stored_ratio_multiply') for model in ('strict','stored_pc64')]+[
                  dict(kind='resident_ratio_multiply',model='stored_pc64',prefix_role='unscaled_prefactor')]
    rows=[];groups=[]
    for g,anchor in enumerate(doc['records']):
        a,y=binary64(anchor['shape']),binary64(anchor['input'])
        prefix=binary64(anchor['density']['bits'])*y/a
        if BINARY64.round(prefix)!=prefix:raise ValueError('anchor rescaling was not exact')
        groups.append(dict(shape=anchor['shape'],ratio=anchor['input'],density=anchor['density'],prefix_publication=bits(float(prefix))))
        for beta in (3.,5.,7.,10.):
            x=y*int(beta)
            if binary64(bits(float(x)))!=x:raise ValueError('scaled input rounded')
            args=[bits(float(x)),anchor['shape'],bits(beta),bits(0.)]
            rows.append(dict(group=g,args=args))
    original=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz';lookup={};wanted={tuple(r['args']) for r in rows}
    for line in gzip.open(original,'rt'):
        raw=json.loads(line);aa=raw['args']
        if len(aa)==1 and len(aa[0])==4 and all(a['kind']=='number' for a in aa[0]):
            key=tuple(a['bits_hex'].removeprefix('0x').lower() for a in aa[0])
            if key in wanted:lookup[key]=raw
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-scaled-consumer-questions/1',complete=True,
        rows=rows,groups=groups,profiles=profiles,original_counterparts=list(lookup.values()),
        source_bindings={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(original.resolve()):sha(original)},
        selection='Four fixed non-dyadic scales for each prior power-of-two input/shape. x/beta is exact and unchanged; original matching bank questions are reused, never re-asked.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    missing=[r for r in rows if tuple(r['args']) not in lookup]
    fresh=service.ask('GAMMA.DIST',[[client.number(x) for x in r['args']] for r in missing])['rows'] if missing else []
    answers={tuple(r['args']):a for r,a in zip(missing,fresh)}
    for key,raw in lookup.items():
        if raw['outcome']['kind']!='number':raise ValueError('nonnumeric historical scaled density')
        answers[key]=dict(type='number',bits=raw['outcome']['bits_hex'].removeprefix('0x'))
    accum=[]
    for g in groups:
        d,y=binary64(g['density']['bits']),binary64(g['ratio']);r=binary64(g['prefix_publication'])
        accum.append([inverse(Span(d,d),y,'div',True,p['model']) if p['kind']=='divide'
                      else BINARY64.pullback(Span(r,r)) for p in profiles])
    dead=[[v is None for v in group] for group in accum];records=[]
    for row in rows:
        answer=answers[tuple(row['args'])]
        if answer['type']!='number':raise ValueError('nonnumeric scaled-density answer')
        x,a=map(binary64,row['args'][:2]);y=binary64(answer['bits']);domains=[]
        for index,p in enumerate(profiles):
            if p['kind']=='divide':constant=x;op='div'
            else:
                constant=(EXTENDED64.round(a/x) if p['kind']=='resident_ratio_multiply'
                          else operation(a,x,'div',p['model']));op='mul'
            domain=inverse(Span(y,y),constant,op,True,p['model']);domains.append(domain.json() if domain is not None else None)
            gi=row['group']
            if not dead[gi][index]:
                accum[gi][index]=None if domain is None else accum[gi][index].meet(domain)
                if accum[gi][index] is None:dead[gi][index]=True
        records.append(dict(row,answer=answer,original_bank=lookup.get(tuple(row['args'])),required_prefixes=domains))
    counts=[Counter() for _ in profiles]
    for gi,g in enumerate(groups):
        g['prefix_constraints']=[]
        for i,p in enumerate(profiles):
            domain=accum[gi][i];b=BINARY64.coordinates(domain) if domain is not None else None;e=EXTENDED64.coordinates(domain) if domain is not None else None
            counts[i]['groups']+=1;counts[i]['binary64_nonempty']+=b is not None;counts[i]['resident64_nonempty']+=e is not None
            g['prefix_constraints'].append(dict(required_shared_prefix=domain.json() if domain is not None else None,
                binary64_count=0 if b is None else b[1]-b[0]+1,resident64_count=0 if e is None else e[1]-e[0]+1))
    result=dict(format='dnacalc-gamma-scaled-final-consumer-preimages/1',complete=True,profiles=profiles,counts=counts,groups=groups,records=records,
        receipts=doc['receipts']+service.receipts,source_bindings={**STAGE_BINDINGS,**bindings,str(source.resolve()):sha(source),str(original.resolve()):sha(original)},
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='One arbitrary shared prefix per fixed shape/exact x/beta. For division it is '
        'the post-alpha product S, constrained directly by old D=S/y; for multiplication it '
        'is unscaled R, whose publication is fixed by the exact alpha/y power-of-two anchor. '
        'This allows every possible rounded alpha product even for fractional shapes. '
        'Both binary64 and resident64 possibilities '
        'are counted without choosing a value. Original bank counterparts retain their own '
        'environment. A contradiction refutes only that shared-prefix final-consumer model; '
        'a nonempty intersection identifies neither a prefix producer nor a complete gamma '
        'algorithm. No store, coefficient, threshold or per-row helper is installed.')
    write(output,result);print(json.dumps(dict(new_calls=len(fresh),original_bank_calls=len(lookup),profiles=list(zip(profiles,counts)))),flush=True);return result


def density_scale_holdout(checkout,bank,output,cache,offline=False):
    """Fresh fractional shapes with alpha/x an exact power of two at the anchor."""
    import gzip
    shapes=(.375,.625,.875,1.125,1.375,1.625,2.25,2.75,3.5,5.5,7.5,12.5,24.,48.)
    rows=[[bits(a*2.**e),bits(a),bits(1.),bits(0.)] for a in shapes for e in (-4,-1,0,1,2)]
    original=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz';wanted={tuple(r) for r in rows};known={}
    for line in gzip.open(original,'rt'):
        raw=json.loads(line);aa=raw['args']
        if len(aa)==1 and len(aa[0])==4 and all(a['kind']=='number' for a in aa[0]):
            key=tuple(a['bits_hex'].removeprefix('0x').lower() for a in aa[0])
            if key in wanted:known[key]=raw
    anchors=Path(str(output)+'.anchors.json')
    write(str(anchors)+'.questions.json',dict(format='dnacalc-gamma-scale-holdout-questions/1',complete=True,
        rows=rows,shapes=shapes,original_bank_counterparts=list(known.values()),source_bindings=STAGE_BINDINGS,
        selection='Fourteen new fixed fractional/integer shapes, five exact alpha/x powers of two. All historical matches are reused before asking new anchors; every anchor then gets all four predeclared scales.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    missing=[r for r in rows if tuple(r) not in known]
    fresh=service.ask('GAMMA.DIST',[[client.number(x) for x in r] for r in missing])['rows'] if missing else []
    answers={tuple(r):a for r,a in zip(missing,fresh)}
    for key,raw in known.items():
        if raw['outcome']['kind']!='number':raise ValueError('nonnumeric original anchor')
        answers[key]=dict(type='number',bits=raw['outcome']['bits_hex'].removeprefix('0x'))
    records=[]
    for r in rows:
        a=answers[tuple(r)]
        if a['type']!='number' or value(a['bits'])<=0:raise ValueError('positive normal holdout anchor required')
        records.append(dict(shape=r[1],input=r[0],density=a,original_bank=known.get(tuple(r))))
    write(anchors,dict(format='dnacalc-gamma-scale-holdout-anchors/1',complete=True,records=records,
        receipts=service.receipts,source_bindings={**STAGE_BINDINGS,str(original.resolve()):sha(original)},
        questions_sha256=sha(str(anchors)+'.questions.json'),
        limits='Original positive density anchors only; exact power-of-two scaling constrains '
        'the prefix publication without a logarithm, exponential or gamma approximation. '
        'No response selects a shape, scale, coefficient or store programme.'))
    return scaled_density_consumers(checkout,anchors,bank,output,cache,offline)


def density_exponents(checkout,bank,output,cache,offline=False):
    """Fixed normal-input, ratio-overflow and zero controls before a consumer change."""
    import gzip
    shapes=(.375,.5,1.,1.25,2.,4.,16.,64.)
    xs=[0.,2.**-1022,2.**-1000,2.**-500,2.**-100,2.**-10,1.,2.**100,2.**500,value('7fefffffffffffff')]
    rows=[[bits(x),bits(a),bits(b),bits(0.)] for a in shapes for x in xs for b in (2.**-1022,1.,value('7fefffffffffffff'))]
    for a in (4.,16.,64.):
        center=int(bits(a/value('7fefffffffffffff')),16)
        for k in range(-4,5):
            x=f'{center+k:016x}'
            if value(x)>=2.**-1022:
                rows.extend([[x,bits(a),bits(b),bits(0.)] for b in (2.**-1022,1.,value('7fefffffffffffff'))])
    rows=list(dict.fromkeys(map(tuple,rows)));wanted=set(rows);known={}
    original=Path(bank)/'cache/oracle/build-20131/GAMMA.DIST.jsonl.gz'
    for line in gzip.open(original,'rt'):
        raw=json.loads(line);aa=raw['args']
        if len(aa)==1 and len(aa[0])==4 and all(a['kind']=='number' for a in aa[0]):
            key=tuple(a['bits_hex'].removeprefix('0x').lower() for a in aa[0])
            if key in wanted:known[key]=raw
    write(str(output)+'.questions.json',dict(format='dnacalc-gamma-density-exponent-questions/1',complete=True,
        rows=rows,original_bank_counterparts=list(known.values()),source_bindings=STAGE_BINDINGS,
        selection='Fixed positive shapes, zero and normal exponent edges, extreme scales and adjacent normal doubles around alpha/MAX. No subnormal input or output-dependent cutoff.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    missing=[r for r in rows if r not in known]
    fresh=service.ask('GAMMA.DIST',[[client.number(x) for x in r] for r in missing])['rows'] if missing else []
    records=[dict(args=list(r),answer=a,original_bank=None) for r,a in zip(missing,fresh)]
    for r,raw in known.items():
        o=raw['outcome'];a=(dict(type='number',bits=o['bits_hex'].removeprefix('0x')) if o['kind']=='number'
                          else dict(type='error',code=o['code']))
        records.append(dict(args=list(r),answer=a,original_bank=raw))
    result=dict(format='dnacalc-gamma-density-exponents/1',complete=True,records=records,receipts=service.receipts,
        source_bindings={**STAGE_BINDINGS,str(original.resolve()):sha(original)},questions_sha256=sha(str(output)+'.questions.json'),
        limits='Original numeric/error readings at predeclared admitted input extremes. No prediction '
        'is an Excel answer; exponent/ratio observations justify no unmeasured numerical cutoff.')
    write(output,result);print(json.dumps(dict(new_calls=len(fresh),original_calls=len(known),kinds=Counter(r['answer']['type'] for r in records))),flush=True);return result


def classify_consumer_losses(prefixes,output):
    """Classify every verified old match without inventing a private prefix."""
    from exact import operation
    source=Path(prefixes);doc=json.loads(source.read_text())
    if doc.get('complete') is not True:raise ValueError('incomplete registered prefix replay')
    records=[];counts=Counter()
    for row in doc['records']:
        if row['classification']!='prefix_replay_verified':raise ValueError('unresolved typed loss requires a separate route audit')
        old=row['original'];want=old['expected'];r=value(row['prefix']);factor=value(row['stored_factor'])
        if not want.startswith('0x') or not math.isfinite(r) or not math.isfinite(factor) or factor==0:
            classification='exception_or_scale_domain';constraints=[];alternatives=[]
        else:
            target=binary64(want);rr=binary64(row['prefix'])
            parse=lambda a:binary64(bits(float(a))) if isinstance(a,bool) else binary64(a)
            x,a=map(parse,old['args'][:2]);constraints=[];alternatives=[];any_binary64=False
            for factor_model in ('strict','stored_pc64'):
                f=operation(a,x,'div',factor_model)
                if factor_model=='stored_pc64' and f!=binary64(row['stored_factor']):raise ValueError('independent factor and registered replay differ')
                for multiply_model in ('strict','stored_pc64'):
                    predicted=operation(rr,f,'mul',multiply_model)
                    domain=inverse(Span(target,target),f,'mul',True,multiply_model)
                    coords=BINARY64.coordinates(domain) if domain is not None else None
                    any_binary64|=coords is not None
                    possible=domain is not None and domain.contains(rr)
                    if possible!=(predicted==target):raise ValueError('prefix inverse and actual arithmetic disagree')
                    if possible:alternatives.append(dict(factor_word=factor_model,multiply_word=multiply_model))
                    constraints.append(dict(factor_word=factor_model,multiply_word=multiply_model,
                        factor=bits(float(f)),prediction=bits(float(predicted)),
                        required_prefix=domain.json() if domain is not None else None,
                        binary64_prefix_count=0 if coords is None else coords[1]-coords[0]+1,
                        current_prefix_possible=possible))
            if any(p['factor_word']==p['multiply_word']=='stored_pc64' for p in alternatives):
                raise ValueError('listed loss still matches through the current consumer')
            classification=('arithmetic_words_remain_open' if alternatives else
                'current_prefix_outside_all_declared_word_preimages' if any_binary64 else
                'prefix_return_storage_remains_open')
        counts[classification]+=1;records.append(dict(row,loss_class=classification,word_alternatives=alternatives,prefix_constraints=constraints))
    result=dict(format='dnacalc-gamma-consumer-loss-classification/1',complete=True,counts=counts,records=records,
        source_bindings={**STAGE_BINDINGS,str(source.resolve()):sha(source)},
        limits='Each old matched return and new registered return is verified by the parent '
        'cleared-prefix replay. Necessary preimages allow both quotient/multiply words. '
        'Word alternatives are explicitly unresolved, never silently attributed to a prefix. '
        'Otherwise the retained prefix lies outside every declared preimage; its logarithm, '
        'exponential, normalizer or return-storage cause is not uniquely identified. Original '
        'answers stay unchanged and every loss stays a failure; no correction or selector is installed.')
    write(output,result);print(json.dumps(dict(rows=len(records),counts=counts)),flush=True);return result


def gammaln_near_stores(checkout,manifest,output,cache,library,offline=False,count=4096):
    """Resolve the existing fixed polynomial crossings with adjacent inputs."""
    import struct
    try:
        import mpmath as mp
    except ImportError:
        sys.path.insert(0,'/home/ubuntu/claude-work/gamma/pylib')
        import mpmath as mp
    mp.mp.prec=200
    source=Path(manifest);prior=json.loads(source.read_text())
    if prior.get('complete') is not True:raise ValueError('incomplete earlier edge manifest')
    frozen={**STAGE_BINDINGS,str(source.resolve()):sha(source),str(Path(library).resolve()):sha(library),
            str(Path(mp.__file__).resolve()):sha(mp.__file__)}
    windows=[];questions=[]
    for edge in prior['edges']:
        center=int(edge['estimated_edge'],16);xs=[f'{center+i-count//2:016x}' for i in range(count)]
        windows.append(dict(polynomial=edge['polynomial'],level=edge['level'],center=edge['estimated_edge'],inputs=xs))
        questions.extend([[dict(type='number',bits=x)] for x in xs])
    write(str(output)+'.questions.json',dict(format='dnacalc-gammaln-near-store-questions/1',complete=True,
        windows=windows,source_bindings=frozen,selection='Four earlier polynomial crossings, each centered on4096 consecutive binary64 inputs. No new answer selects a center or a model.'))
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    answers=service.ask('GAMMALN',questions)['rows'];native=Native(library);records=[]
    graphs=prior['graphs'];remainder=mp.mpf(0)
    for wi,window in enumerate(windows):
        center=mp.mpf(value(window['center']));f0=mp.loggamma(center);f1=mp.digamma(center);f2=mp.polygamma(1,center)
        samples=[]
        for index,(x,a) in enumerate(zip(window['inputs'],answers[wi*count:(wi+1)*count])):
            if a['type']!='number':raise ValueError('nonnumeric positive gamma observation')
            dx=mp.mpf(value(x))-center;truth=f0+f1*dx+f2*dx*dx/2
            unit=mp.mpf(math.ulp(value(a['bits'])));remainder=max(remainder,abs(dx)**3/(6*unit))
            predictions={name:f'{native.publish(native.evaluate(graph,[int(x,16)],0x133f,mask=1)[:10]):016x}' for name,graph in graphs.items()}
            ys=dict(Excel=a['bits'],**predictions)
            residuals={name:float((mp.mpf(value(raw))-truth)/unit) for name,raw in ys.items()}
            samples.append(dict(offset=index-count//2,input=x,excel=a,predictions=predictions,error_in_result_ulps=residuals))
        blocks=[]
        for start in range(0,count,128):
            part=samples[start:start+128]
            blocks.append(dict(first_offset=part[0]['offset'],last_offset=part[-1]['offset'],
                means={name:sum(p['error_in_result_ulps'][name] for p in part)/len(part) for name in ys}))
        records.append(dict(polynomial=window['polynomial'],level=window['level'],center=window['center'],samples=samples,blocks=blocks))
    result=dict(format='dnacalc-gammaln-near-store-windows/1',complete=True,records=records,
        source_bindings=frozen,receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        quadratic_remainder_bound_in_ulps=str(remainder),
        limits='Exact adjacent-input observations at previously chosen centers; block means are descriptive, '
        'not independent random samples or a fitted change point. The old coefficient table and '
        'all three predeclared store graphs remain frozen. On2<=x<=3, abs(psi2(x))<1 bounds '
        'the stated Taylor remainder. Neither a table nor a store change follows from a pattern alone.')
    write(output,result);print(json.dumps(dict(windows=len(records),calls=len(answers),remainder=str(remainder))),flush=True);return result


def stirling_stage_constraints(checkout, output, cache, library, offline=False, count=64, consumers=False):
    """Invert the last GAMMA product without choosing Stirling coefficients.

    Independent public LN and EXP leaves remove the corresponding producer
    assumptions. Each declared stage graph admits an arbitrary binary64 w(x).
    Common-window and monotone constraints are additional, explicitly
    conditional tests, not fitted values or a replacement for the application.
    """
    if not 2 <= count <= 256:
        raise ValueError('small adjacent windows required')
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    native = Native(library)
    centers = (10., 12., 16., 20., 24., 30., 32., 36., 48., 64., 80., 100., 128., 160.)
    prefix_profiles = ['resident_product_resident_sqrt', 'stored_product_resident_sqrt',
                       'resident_product_stored_sqrt', 'stored_product_stored_sqrt']
    forms = ('constant_then_correction', 'correction_then_constant', 'stored_constant_correction',
             'resident_constant_correction', 'resident_constant_body', 'resident_square_body') if consumers else ('constant_then_correction',)
    profiles = [dict(prefix=p,finish=f) for p in prefix_profiles for f in forms]
    groups = [dict(center=bits(c), inputs=[f'{int(bits(c),16)+i:016x}' for i in range(-count//2,count//2+1) if i])
              for c in centers]
    # Stay in the high band even at its first exact boundary.
    groups[0]['inputs'] = [f'{int(bits(10.),16)+i:016x}' for i in range(1,count+1)]
    xs = [x for group in groups for x in group['inputs']]
    frozen = {**STAGE_BINDINGS, str(Path(library).resolve()):sha(library)}
    write(str(output)+'.questions.json', dict(format='dnacalc-gamma-stirling-stage-questions/1',
        complete=True, groups=groups, profiles=profiles, source_bindings=frozen,
        selection='Fourteen fixed integer and exponent boundaries, consecutive inputs in increasing order. '
                  'Four literal product/square-root store graphs precede every new answer. '
                  'All logarithms and exponentials are independently measured public calls.'))
    logs = dict(zip(xs,service.ask('LN',[[client.number(x)] for x in xs])['rows']))
    ems = dict(zip(xs,service.ask('EXP',[[client.number(-value(x))] for x in xs])['rows']))
    answers = dict(zip(xs,service.ask('GAMMA',[[client.number(x)] for x in xs])['rows']))
    def op(kind, a, b):
        return native.leaf(kind, a, b, 0x133f)
    def scalar(x):
        return load64(int(bits(x),16))
    def stored(a):
        return load64(native.publish(a))
    root_arguments = {}
    for x in xs:
        half = .5*value(x)
        root_arguments[x] = f'{native.publish(op("mul",scalar(half-math.floor(half)),load64(int(logs[x]["bits"],16)))):016x}'
    roots_unique = sorted(set(root_arguments.values()))
    roots = dict(zip(roots_unique,service.ask('EXP',[[client.number(x)] for x in roots_unique])['rows']))
    def power(x, n):
        product, square = None, load64(int(x,16))
        while n:
            if n & 1:
                product = square if product is None else stored(op('mul',product,square))
            n >>= 1
            if n:
                square = stored(op('mul',square,square))
        return product if product is not None else scalar(1.)
    records = []
    counts = [Counter() for _ in profiles]
    for group in groups:
        common = [None for _ in profiles]
        common_dead = [False for _ in profiles]
        monotone_upper = [None for _ in profiles]
        monotone_dead = [False for _ in profiles]
        group_records = []
        for x in group['inputs']:
            target = answers[x]
            if any(a['type'] != 'number' for a in (target,logs[x],ems[x],roots[root_arguments[x]])):
                raise ValueError('high-band leaves must be finite numbers')
            v = stored(op('mul',power(x,math.floor(.5*value(x))),load64(int(roots[root_arguments[x]]['bits'],16))))
            square_root = op('sqrt',load64(int(x,16)),scalar(0.))
            product = op('mul',v,load64(int(ems[x]['bits'],16)))
            domains = []
            for index, profile in enumerate(profiles):
                label, finish = profile['prefix'], profile['finish']
                pp = stored(product) if label.startswith('stored_product') else product
                qq = stored(square_root) if label.endswith('stored_sqrt') else square_root
                u = stored(op('div',pp,qq))
                square = stored(op('mul',u,v))
                constant = load64(0x40040d931ff62706)
                body = stored(op('mul',square,constant))
                interval = Span(binary64(target['bits']),binary64(target['bits']))
                if finish in ('constant_then_correction','resident_constant_body','resident_square_body'):
                    if finish == 'resident_constant_body':
                        body = op('mul',square,constant)
                    elif finish == 'resident_square_body':
                        body = stored(op('mul',op('mul',u,v),constant))
                    domain = inverse(interval,decode(body),'mul',True,'stored_pc64')
                elif finish == 'correction_then_constant':
                    interval = inverse(interval,decode(constant),'mul',True,'stored_pc64')
                    domain = inverse(interval,decode(square),'mul',True,'stored_pc64') if interval else None
                elif finish in ('stored_constant_correction','resident_constant_correction'):
                    interval = inverse(interval,decode(square),'mul',True,'stored_pc64')
                    model = 'stored_pc64' if finish.startswith('stored_') else 'resident_pc64'
                    domain = inverse(interval,decode(constant),'mul',True,model) if interval else None
                else:
                    raise ValueError('unknown fixed final association')
                coordinates = BINARY64.coordinates(domain) if domain else None
                counts[index]['rows'] += 1
                counts[index]['empty_binary64_w'] += coordinates is None
                if not common_dead[index]:
                    common[index] = domain if common[index] is None else common[index].meet(domain) if domain else None
                    common_dead[index] = common[index] is None
                if not monotone_dead[index]:
                    if coordinates is None:
                        monotone_dead[index] = True
                    else:
                        low, high = coordinates
                        upper = high if monotone_upper[index] is None else min(monotone_upper[index],high)
                        monotone_dead[index] = upper < low
                        monotone_upper[index] = upper
                domains.append(dict(body=f'{native.publish(body):016x}',resident_body=body.hex(),square=f'{native.publish(square):016x}',
                    required_w=domain.json() if domain else None,
                    binary64_count=0 if coordinates is None else coordinates[1]-coordinates[0]+1,
                    binary64_values=[] if coordinates is None else [bits(float(BINARY64.value(k)))
                        for k in range(coordinates[0],min(coordinates[1]+1,coordinates[0]+4))]))
            group_records.append(dict(input=x,answer=target,public_log=logs[x],public_negative_exp=ems[x],
                root_argument=root_arguments[x],public_root=roots[root_arguments[x]],
                powered=f'{native.publish(v):016x}',resident_sqrt=square_root.hex(),profiles=domains))
        results = []
        for index in range(len(profiles)):
            coordinates = BINARY64.coordinates(common[index]) if common[index] else None
            counts[index]['groups'] += 1
            counts[index]['common_w_possible'] += coordinates is not None
            counts[index]['monotone_w_possible'] += not monotone_dead[index]
            results.append(dict(common_w=common[index].json() if common[index] else None,
                common_binary64_count=0 if coordinates is None else coordinates[1]-coordinates[0]+1,
                monotone_w_possible=not monotone_dead[index]))
        records.append(dict(center=group['center'],samples=group_records,profiles=results))
    result = dict(format='dnacalc-gamma-stirling-stage-constraints/1',complete=True,
        profiles=profiles,counts=counts,records=records,receipts=service.receipts,source_bindings=frozen,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Exact final-product preimages with an arbitrary binary64 correction per input. '
               'Empty cells refute the declared preceding body and final product. '
               'Common or nonincreasing corrections across a window are stronger conditional hypotheses; '
               'a failed one does not refute an arbitrary nonmonotonic computed correction. '
               'The source integer-power route and stored public LN/EXP leaves stay fixed. '
               'No correction, coefficient or stage is selected by counts or installed.')
    write(output,result)
    print(json.dumps(dict(profiles=profiles,counts=counts)),flush=True)
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
    p.add_argument('--density-primitives')
    p.add_argument('--shared-log-constraints')
    p.add_argument('--scaled-log-constraints')
    p.add_argument('--split-prefactor')
    p.add_argument('--poisson-pairs')
    p.add_argument('--gammaln-store-edges', action='store_true')
    p.add_argument('--gammaln-near-stores',metavar='MANIFEST')
    p.add_argument('--unit-shape-exp',nargs=2,metavar=('BANK','GRAPHS'))
    p.add_argument('--unit-prefactor',metavar='STUDY')
    p.add_argument('--unit-constants',metavar='STUDY')
    p.add_argument('--unit-consumers',metavar='STUDY')
    p.add_argument('--scaled-consumers',nargs=2,metavar=('STUDY','BANK'))
    p.add_argument('--scale-holdout',metavar='BANK')
    p.add_argument('--density-exponents',metavar='BANK')
    p.add_argument('--classify-losses',metavar='PREFIX_REPLAY')
    p.add_argument('--stirling-stages',action='store_true')
    p.add_argument('--stirling-consumers',action='store_true')
    args = p.parse_args()
    if not 2 <= args.count <= (4096 if args.gammaln_store_edges or args.gammaln_near_stores else 256):
        raise ValueError('small adjacent windows only')
    if args.stirling_stages or args.stirling_consumers:
        stirling_stage_constraints(args.checkout,args.output,args.cache,args.library,args.offline,args.count,args.stirling_consumers)
    elif args.classify_losses:
        classify_consumer_losses(args.classify_losses,args.output)
    elif args.density_exponents:
        density_exponents(args.checkout,args.density_exponents,args.output,args.cache,args.offline)
    elif args.scale_holdout:
        density_scale_holdout(args.checkout,args.scale_holdout,args.output,args.cache,args.offline)
    elif args.scaled_consumers:
        scaled_density_consumers(args.checkout,args.scaled_consumers[0],args.scaled_consumers[1],args.output,args.cache,args.offline)
    elif args.unit_consumers:
        unit_shape_consumers(args.unit_consumers,args.output,args.cache)
    elif args.unit_constants:
        unit_shape_constants(args.unit_constants,args.output,args.cache,args.library)
    elif args.unit_prefactor:
        unit_shape_prefactor(args.checkout,args.unit_prefactor,args.output,args.cache,args.library,args.offline)
    elif args.unit_shape_exp:
        unit_shape_exponential(args.checkout,args.unit_shape_exp[0],args.unit_shape_exp[1],args.output,args.cache,args.library,args.offline)
    elif args.gammaln_near_stores:
        gammaln_near_stores(args.checkout,args.gammaln_near_stores,args.output,args.cache,args.library,args.offline,args.count)
    elif args.gammaln_store_edges:
        gammaln_store_edges(args.checkout, args.output, args.cache, args.library, args.offline, args.count)
    elif args.poisson_pairs:
        poisson_pairs(args.checkout, args.poisson_pairs, args.output, args.cache, args.library, args.offline)
    elif args.split_prefactor:
        split_prefactor(args.checkout, args.split_prefactor, args.output, args.cache, args.library, args.offline)
    elif args.scaled_log_constraints:
        shared_log_constraints(args.scaled_log_constraints,args.output,args.cache,args.library,scaled_log=True)
    elif args.shared_log_constraints:
        shared_log_constraints(args.shared_log_constraints, args.output, args.cache, args.library)
    elif args.density_primitives:
        power_primitives(args.density_primitives,args.output,args.cache,args.library,direct_density=True)
    elif args.power_primitives:
        power_primitives(args.power_primitives, args.output, args.cache, args.library)
    elif args.power_prefactor:
        power_prefactor(args.checkout, args.output, args.cache, args.library, args.offline)
    else:
        paired(args.checkout, args.output, args.cache, args.library, args.count, args.offline)
