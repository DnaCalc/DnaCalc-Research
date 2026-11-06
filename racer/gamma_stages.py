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


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('checkout')
    p.add_argument('output')
    p.add_argument('--cache', required=True)
    p.add_argument('--library', required=True)
    p.add_argument('--count', type=int, default=32)
    p.add_argument('--offline', action='store_true')
    p.add_argument('--power-prefactor', action='store_true')
    args = p.parse_args()
    if not 2 <= args.count <= 256:
        raise ValueError('small adjacent windows only')
    if args.power_prefactor:
        power_prefactor(args.checkout, args.output, args.cache, args.library, args.offline)
    else:
        paired(args.checkout, args.output, args.cache, args.library, args.count, args.offline)
