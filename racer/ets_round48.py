"""R48 ETS structural experiments. Original MIT, implementer, tier 1 only.

Sources: original Excel20527 questions, IEEE754, Hyndman's published
mathematics and Microsoft's 2014 Power View article (links in design note).
No restricted code, disassembly, or excluded checkout is read.
"""
import argparse
from collections import Counter
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path
import random

import numpy as np

from . import ets, ets_round47 as prior, round41
from .fingerprint import bits, value

OUT = Path('/home/ubuntu/codex-runs/round48')


def save(name, obj):
    with (OUT / (name + '.json')).open('x') as stream:
        json.dump(obj, stream, separators=(',', ':'), allow_nan=False)


def freeze(name, obj):
    source = Path(__file__).read_bytes()
    save(name + '-frozen', dict(obj, source=source.decode(),
                              source_sha256=hashlib.sha256(source).hexdigest()))


def service():
    round41.OUT = OUT
    return round41.service(priority='interactive')


def assess(predictions, answers):
    failures = [dict(index=i, prediction=p, answer=a) for i, (p, a) in
                enumerate(zip(predictions, answers)) if p != a.get('bits')]
    if len(predictions) != len(answers):
        raise ValueError('incomplete answers')
    return dict(selections=len(answers), failures=failures)


def project(pair):
    a = min(.998, max(.002, pair[0]))
    return a, min(a - .001, max(.001, pair[1]))


def nelder_mead(score, start, step=.1, iterations=160):
    """Textbook two-dimensional reflection/expansion/contraction/shrink.

    Feasibility is enforced at every stored vertex. Coefficients 1,2,1/2,1/2;
    fixed iteration budget avoids an unidentified floating tolerance.
    """
    points = [project(start), project((start[0]+step, start[1])),
              project((start[0], start[1]+step))]
    pairs = [(score(p), p) for p in points]
    for _ in range(iterations):
        pairs.sort(key=lambda p: p[0])
        best, middle, worst = [p[1] for p in pairs]
        center = tuple((a+b)/2 for a, b in zip(best, middle))
        reflected = project(tuple(2*c-w for c, w in zip(center, worst)))
        sr = score(reflected)
        if sr < pairs[0][0]:
            expanded = project(tuple(c+2*(r-c) for c, r in zip(center, reflected)))
            se = score(expanded)
            pairs[-1] = (se, expanded) if se < sr else (sr, reflected)
        elif sr < pairs[1][0]:
            pairs[-1] = sr, reflected
        else:
            target = reflected if sr < pairs[2][0] else worst
            contracted = project(tuple((c+t)/2 for c, t in zip(center, target)))
            sc = score(contracted)
            if sc < min(sr, pairs[2][0]):
                pairs[-1] = sc, contracted
            else:
                points = [project(tuple((a+b)/2 for a, b in zip(best, p)))
                          for p in (middle, worst)]
                pairs = [pairs[0]] + [(score(p), p) for p in points]
    return min(pairs, key=lambda p: p[0])[1]


def score(y, pair, window='hybrid', error='additive'):
    """Profiled Gaussian additive or multiplicative innovation likelihood.

    For a fixed residual window, additive likelihood is monotone in MSE.
    Multiplicative likelihood adds the observation Jacobian; it is a distinct
    model and explicitly conditional on the additive state recurrence here.
    """
    a, b = pair
    n = len(y)
    start = n - (n if n <= 10 else min(int(n*.4), 10))
    l, t, seasons = y[0], ets.initial_trend(y), [0., 0.]
    residuals = []; logs = []
    for j in range(1, n):
        p = (l+t)+seasons[j % 2]
        e = y[j]-p
        if window != 'tail' or j >= start:
            if error == 'multiplicative':
                if p <= 0: return math.inf
                residuals.append((e/p)**2); logs.append(math.log(p))
            else:
                residuals.append(e*e)
        if window != 'full' and start and j >= start:
            l += t
        else:
            l += t+a*e; t += b*e; seasons[j % 2] += 2.**-52*e
    mse = sum(residuals)/len(residuals)
    if error == 'additive': return mse
    return len(residuals)*math.log(mse)+2*sum(logs) if mse else -math.inf


def refine(score_fn, start):
    """Coarse winner followed by four fixed halved coordinate meshes."""
    best = start
    for width in [.1, .05, .025, .0125]:
        candidates = [project((best[0]+dx*width, best[1]+dy*width))
                      for dx, dy in itertools.product([-1, 0, 1], repeat=2)]
        best = min(candidates, key=score_fn)
    return best


def search_shapes():
    rng = random.Random(480101)
    cases = []
    for i in range(96):
        n = [8, 11, 16, 24, 31, 32][i % 6]
        # Positive offset permits the alternative multiplicative likelihood.
        y = [256.+4*j+rng.randrange(-128, 129)/8 for j in range(n)]
        if i % 4 == 1: y[n//2] += 128.
        elif i % 4 == 2: y = [v+(64. if j >= n//2 else 0.) for j, v in enumerate(y)]
        elif i % 4 == 3: y = [v+16*math.sin(2*math.pi*j/5) for j, v in enumerate(y)]
        cases.append(dict(y=y, season=0, family=i % 4))
    pools = dict(sparse=prior.POOL,
                 cartesian=sorted(set(project((a, b)) for a, b in
                     itertools.product([.1, .25, .5, .75, .9], repeat=2))))
    predictions = {}; raw = {}
    native = prior.Native()
    for window, error, pool_name, method in itertools.product(
            ['hybrid', 'tail', 'full'], ['additive', 'multiplicative'],
            ['sparse', 'cartesian'], ['grid', 'refine', 'simplex']):
        key = ':'.join([window, error, pool_name, method]); pairs = []; pred = []
        for z in cases:
            y = z['y']; fn = lambda p: score(y, p, window, error)
            chosen = min(pools[pool_name], key=fn)
            if method == 'refine': chosen = refine(fn, chosen)
            elif method == 'simplex': chosen = nelder_mead(fn, chosen)
            rounded = tuple(float(format(v, '.6g')) for v in chosen)
            f = native.fit(y, [rounded]); pairs.append(list(chosen))
            pred.extend([bits(f[1]), bits(f[2])])
        predictions[key] = pred; raw[key] = pairs
    distinctions = {}
    for name, left, right in [
        ('grid_vs_simplex', 'hybrid:additive:sparse:grid', 'hybrid:additive:sparse:simplex'),
        ('additive_vs_multiplicative', 'hybrid:additive:sparse:grid', 'hybrid:multiplicative:sparse:grid'),
        ('training_vs_holdout', 'full:additive:sparse:grid', 'tail:additive:sparse:grid')]:
        distinctions[name] = [i for i in range(len(cases)) if
                             predictions[left][2*i:2*i+2] != predictions[right][2*i:2*i+2]]
    freeze('search-shapes', dict(cases=cases, predictions=predictions, raw=raw,
        distinctions=distinctions, pools=pools,
        identity='For the same additive innovations and window: -2 log L = n log(SSE/n) + constant. No real-valued separating history exists.',
        rank=['finite sparse grid with held-out transitions', 'coarse grid then refinement',
              'bounded two-dimensional Nelder-Mead', 'multiplicative error likelihood'],
        limits='These complete finite trees include their initial state and score. Failure refutes the tree, not every optimizer in that family.'))
    c, s = service(); d = prior.ask(c, s, cases, 'FORECAST.ETS.STAT', [1, 2])
    save('search-shapes-answers', d)
    results = {k: assess(p, d['answers']) for k, p in predictions.items()}
    save('search-shapes-result', dict(programs=results,
        survivors=[k for k, r in results.items() if not r['failures']]))
    print('search trees', len(results), 'survivors', [k for k, r in results.items() if not r['failures']],
          'separating histories', {k: len(v) for k, v in distinctions.items()}, flush=True)


def joined_graph(n, form=0, extend=True, initial_store=False, initial_precision=0x133f):
    """Splice existing initializer and recurrence graphs without a f64 boundary.

    Arithmetic stores of the recurrence stay fixed to the established strict
    program, except for the measured PC64 final addition. Only initializer
    stores vary. This directly tests a resident initializer quotient.
    """
    from .ets_loop import graph, slope_graph, TAGS
    init = slope_graph(n, form, extend=extend)
    loop = graph(n, form=3)
    nodes = copy.deepcopy(init['nodes'])
    trend = init['outputs'][0]
    if initial_store:
        nodes.append(dict(op='add', a=trend, b=n+1, site=8))
        trend = len(nodes)-1
    mapping = {i: i for i in range(n)}
    for i, node in enumerate(loop['nodes']):
        if node['op'] == 'input':
            if node['a'] == n+1: mapping[i] = trend
            elif node['a'] < n: mapping[i] = node['a']
            else:
                new = dict(node); new['a'] = n+node['a']-(n+1) if node['a'] > n+1 else n
                # y..., horizon, alpha, beta (initializer replaces trend input).
                nodes.append(new); mapping[i] = len(nodes)-1
        else:
            new = dict(node)
            if node['op'] != 'const':
                new['a'] = mapping[node['a']]; new['b'] = mapping[node['b']]
                tag = node.get('tag')
                new.pop('site', None)
                # Every arithmetic expression computes at PC53 with no wider
                # resident intermediates; publication alone uses PC64.
                new['phase'] = 1 if tag == 'publication' else 0
                if tag != 'publication': new['site'] = 9
            nodes.append(new); mapping[i] = len(nodes)-1
    # Publication always stays PC64. Initializer precision changes independently.
    for node in nodes[:len(init['nodes'])]: node['phase'] = int(initial_precision == 0x133f)
    return dict(nodes=nodes, outputs=[mapping[i] for i in loop['outputs']], sites=10)


def short_resident(centered=False):
    from .x87_programs import Native, build_library
    # Retained random short observations excluding resets in the fixed prior
    # candidate. Reported parameters remain conditional inputs: this selection
    # does not certify whether Excel itself reset a private raw beta.
    d = ets.read(prior.OUT/'short-map-frozen.json')
    actual = ets.read(prior.OUT/'short-map-forecast.json')['answers']
    st = ets.read(prior.OUT/'short-map-stat.json')['answers']
    old = prior.Native(); cases = []
    for i, z in enumerate(d['cases']):
        f = old.fit(z['y'])
        if f[7] or any(bits(f[1+j]) != st[3*i+j].get('bits') for j in range(2)): continue
        cases.append(dict(z, original_index=i, a=value(st[3*i]['bits']), b=value(st[3*i+1]['bits'])))
    lib = OUT/'x87-ets48.so'; build_library(lib); native = Native(lib)
    programs = []; survivors = []
    prefix = 'short-centered-v2' if centered else 'short-resident-v2'
    freeze(prefix, dict(cases=cases, forms=list(range(3,5) if centered else range(3)),
        initializer_source=Path(__file__).with_name('ets_loop.py').read_text(),
        correction='Preserve the original result. V2 varies initializer PC53/PC64 independently of the measured final PC64 publication. The original PC53 half varied both, so it was not the stated fixed-publication family.',
        conditional='Reported fitted parameters and the fixed prior candidate reset decision; not a certification of private raw beta.'))
    for form, initial_store, precision in itertools.product(range(3,5) if centered else range(3), [False, True], [0x123f, 0x133f]):
        graphs = {n: joined_graph(n, form, initial_store=initial_store, initial_precision=precision) for n in range(3, 8)}
        for mask in range(256):
            first = None; failures = 0
            for n in range(3, 8):
                subset = [z for z in cases if len(z['y']) == n]
                if not subset: continue
                inputs = [[int(bits(v), 16) for v in z['y']+[h, z['a'], z['b']]]
                          for z in subset for h in [0., 1., 2., 5.]]
                out = native.grid(graphs[n], inputs, cw=0x123f, inner_cw=0x133f,
                                  mask=mask | (1<<9) | ((1<<8) if initial_store else 0))
                for j, z in enumerate(subset):
                    for k, h in enumerate([0, 1, 2, 5]):
                        prediction = f'{int(out[4*j+k, (n-1+h)%2]):016x}'
                        answer = actual[4*z['original_index']+k]
                        if prediction != answer.get('bits'):
                            failures += 1
                            if first is None: first = dict(case=z['original_index'], horizon=h, prediction=prediction, answer=answer)
            result = dict(form=form, initial_store=initial_store, precision=precision, mask=mask,
                          contradictions=failures, first=first)
            programs.append(result)
            if not failures: survivors.append(result)
    save(prefix+'-result', dict(cases=cases, programs=programs, survivors=survivors,
        limits='Retained discovery only; fixed postfit-reset candidate chooses the reset-free subset. No product behavior is inferred from selected compatibility counts.'))
    print('resident initializer programs', len(programs), 'survivors', len(survivors), survivors[:3], flush=True)


def detection_shapes():
    from .ets_round40 import padded_history
    from .ets_seasonal import detection_candidate, moving_initial
    rng = random.Random(480102); cases = []
    for n in [5, 6, 7, 8, 11, 12, 16, 20, 24, 31, 32]:
        for k in range(20):
            p = 2+k%8; phase = k%p; drift = [0., .125, 16., 65536.][k%4]
            pattern = [rng.randrange(-16, 17)/4 for _ in range(p)]
            y = [pattern[(j+phase)%p]+drift*j for j in range(n)]
            if k >= 10: y = [v+rng.randrange(-8, 9)/16 for v in y]
            ex = [-300, 0, 300][k%3]
            cases.append(dict(y=[math.ldexp(v, ex) for v in y], p=p, season=1, phase=phase, kind='periodic' if k<10 else 'noise'))
    predictions = {}
    for pre, norm, order, cutoff in itertools.product(
            ['linear', 'difference'], ['biased', 'pearson'], ['first', 'global'], [.5, 1.96]):
        key = f'acf:{pre}:{norm}:{order}:{cutoff}'; pred = []
        for z in cases:
            y = padded_history(z['y']) if len(z['y']) < 7 else z['y']
            limit = .5 if cutoff == .5 else cutoff/math.sqrt(len(y))
            p, _ = detection_candidate(y, pre, norm, order, 0, limit)
            pred.append(bits(float(p if p <= 2*len(z['y'])//5 else 0)))
        predictions[key] = pred
    # Spectrum supplies bounded candidate intervals; correlation chooses inside
    # them. These explicit simplified trees test the architecture, not every
    # heuristic in the Power View description.
    for detrend, spectrum_size, rank in itertools.product([True, False], ['n', 'power2'], ['power', 'pearson']):
        key = f'spectrum:{detrend}:{spectrum_size}:{rank}'; pred = []
        for z in cases:
            y = padded_history(z['y']) if len(z['y']) < 7 else z['y']
            y = np.asarray(y)/max(abs(v) for v in y)
            n = len(y); t = ets.initial_trend(y)
            x = y-np.mean(y)-(np.arange(n)-(n-1)/2)*t if detrend else y-np.mean(y)
            size = n if spectrum_size == 'n' else 1 << (n-1).bit_length()
            powers = abs(np.fft.rfft(x, n=size))**2
            top = sorted(range(1, len(powers)), key=lambda j: powers[j], reverse=True)[:3]
            candidates = []
            _, corr = detection_candidate(list(y), 'linear' if detrend else 'mean', 'pearson', 'global', 0, 0.)
            for frequency in top:
                lo = max(2, math.ceil(size/(frequency+1))); hi = min(2*len(z['y'])//5, math.floor(size/max(.5, frequency-1)))
                for p in range(lo, hi+1):
                    if p < len(corr): candidates.append((powers[frequency] if rank == 'power' else corr[p], p))
            p = max(candidates, default=(0, 0), key=lambda r:r[0])[1]
            pred.append(bits(float(p)))
        predictions[key] = pred
    # Direct linear phase regression is an error-minimizing detector without
    # inheriting the still unidentified parameter search/seasonal initializer.
    for scope, penalty in itertools.product(['full', 'holdout'], ['none', 'aic', 'bic']):
        key = f'phase-error:{scope}:{penalty}'; pred = []
        for z in cases:
            y = np.asarray(z['y'])/max(abs(v) for v in z['y']); n=len(y)
            start = n if scope == 'full' else max(3, n-min(int(n*.4), 10))
            candidates = []
            for p in [1]+list(range(2, 2*n//5+1)):
                matrix = np.column_stack([np.arange(n)]+[(np.arange(n)%p == q).astype(float) for q in range(p)])
                coef, _, _, _ = np.linalg.lstsq(matrix[:start], y[:start], rcond=None)
                errors = y-matrix@coef
                residual = errors if scope == 'full' else errors[start:]
                mse = float(np.mean(residual**2)); k=p+1
                cost = mse if penalty == 'none' else len(residual)*math.log(max(mse, 1e-300))+(2*k if penalty == 'aic' else k*math.log(n))
                candidates.append((cost, p))
            p = min(candidates)[1]; pred.append(bits(float(0 if p == 1 else p)))
        predictions[key] = pred
    freeze('detection-shapes', dict(cases=cases, predictions=predictions,
        rank=['spectrum intervals with correlation ranking', 'ACF with a fixed statistical cutoff',
              'period minimizing fitted or held-out phase error'],
        limits='Power View supplies a method-family hypothesis, not an Excel worksheet specification. Explicit finite simplified trees only.'))
    c, s = service(); d = prior.ask(c, s, cases, 'FORECAST.ETS.SEASONALITY', [0])
    save('detection-shapes-answers', d)
    programs = {k: assess(p, d['answers']) for k, p in predictions.items()}
    save('detection-shapes-result', dict(programs=programs,
        survivors=[k for k, r in programs.items() if not r['failures']]))
    print('detection trees', len(programs), 'survivors', [k for k, r in programs.items() if not r['failures']], flush=True)


def short_seasons():
    from .ets_seasonal import initial, fit
    rng = random.Random(480301); cases = []; forecasts = []; stats = []; confidence = []
    for n in [5, 6, 7, 8, 9, 10]:
        for k in range(64):
            exponent = [-380, -100, 0, 100, 380][k % 5]
            offset = rng.randrange(-1024, 1025)/4
            amplitude = rng.choice([-1, 1])*rng.randrange(1, 129)/8
            drift = rng.randrange(-32, 33)/16
            y = [math.ldexp(offset+(-amplitude if j % 2 else amplitude)+drift*j, exponent) for j in range(n)]
            l, t, se, e, _ = fit(y, 2, .751, .001, .001, initial(y, 2, 'firsttwo', 'first_detrend'))
            if any(e): raise ValueError('nonzero constructed innovation')
            cases.append(dict(y=y, season=1, n=n, kind='exact_alternating_affine'))
            forecasts.extend(bits((l+h*t)+se[(n-1+int(h)) % 2]) for h in [0., 1., 2., 5.])
            stats.extend(bits(v) for v in [.751, .001, .001, 0., 0., 0., 0.])
            confidence.extend([bits(0.)]*2)
    freeze('short-seasons', dict(cases=cases, forecast_predictions=forecasts,
        statistic_predictions=stats, confidence_predictions=confidence,
        period_predictions=[bits(2.)]*len(cases),
        method='First-two-cycle exact trend and first-season centered states. All transitions are zero before questions. Lengths straddle the three-cycle initializer boundary.'))
    c, s = service(); results = {}
    for name, function, variants, pred in [
            ('period', 'FORECAST.ETS.SEASONALITY', [0], [bits(2.)]*len(cases)),
            ('forecast', 'FORECAST.ETS', [0., 1., 2., 5.], forecasts),
            ('stat', 'FORECAST.ETS.STAT', list(range(1, 8)), stats),
            ('confidence', 'FORECAST.ETS.CONFINT', [(1., .95), (5., .95)], confidence)]:
        d = prior.ask(c, s, cases, function, variants); save('short-seasons-'+name, d)
        results[name] = assess(pred, d['answers'])
    save('short-seasons-result', results)
    print('short seasonal exact', {k: (v['selections'], len(v['failures']),
         dict(Counter(cases[f['index']//(1 if k=='period' else 4 if k=='forecast' else 7 if k=='stat' else 2)]['n'] for f in v['failures'])))
         for k, v in results.items()}, flush=True)


def alternating_proof():
    from .ets_seasonal import initial, fit
    rng = random.Random(480302); cases=[]; predictions={k:[] for k in ['forecast','stat','confidence','period']}
    for n in range(5, 33):
        bases=[]
        for k in range(16):
            exponent=rng.randrange(-370, 371); offset=rng.randrange(-1024,1025)/8
            amplitude=rng.choice([-1,1])*rng.randrange(1,129)/16; drift=rng.randrange(-32,33)/16
            bases.append([math.ldexp(offset+(-amplitude if j%2 else amplitude)+drift*j,exponent) for j in range(n)])
        # Actual input-magnitude endpoints and adjacent coarse-lattice values.
        bases.extend([[math.ldexp(float(j%2),exponent) for j in range(n)] for exponent in [-400,400]])
        for y in bases:
            for season in [1,2]:
                l,t,se,e,_=fit(y,2,.751,.001,.001,initial(y,2,'firsttwo','first_detrend'))
                if any(e): raise ValueError('proof innovations not zero')
                cases.append(dict(y=y,season=season,n=n,kind='fresh_exact_alternating'))
                predictions['forecast'].extend(bits((l+h*t)+se[(n-1+int(h))%2]) for h in [0.,1.,2.,5.])
                predictions['stat'].extend(bits(v) for v in [.751,.001,.001,0.,0.,0.,0.])
                predictions['confidence'].extend([bits(0.)]*2);predictions['period'].append(bits(2.))
    # Explicit period2 constant/affine controls test the zero-amplitude boundary.
    for n in range(5,33):
        for drift in [0.,.25,-.5]:
            y=[16.+drift*j for j in range(n)]
            l,t,se,e,_=fit(y,2,.751,.001,.001,initial(y,2,'firsttwo','first_detrend'))
            cases.append(dict(y=y,season=2,n=n,kind='zero_amplitude_boundary'))
            predictions['forecast'].extend(bits((l+h*t)+se[(n-1+int(h))%2])for h in[0.,1.,2.,5.])
            predictions['stat'].extend(bits(v)for v in[.751,.001,.001,0.,0.,0.,0.])
            predictions['confidence'].extend([bits(0.)]*2);predictions['period'].append(bits(0.))
    freeze('alternating-proof',dict(cases=cases,predictions=predictions,
        method='Independent seed; every length5..32, automatic/explicit2 pairs, input magnitude endpoints, explicit zero-amplitude controls. Exact zero innovations; no reported parameters condition predictions.'))
    c,s=service();results={}
    for name,fun,variants in [('forecast','FORECAST.ETS',[0.,1.,2.,5.]),
                             ('stat','FORECAST.ETS.STAT',list(range(1,8))),
                             ('confidence','FORECAST.ETS.CONFINT',[(1.,.95),(5.,.95)]),
                             ('period','FORECAST.ETS.SEASONALITY',[0])]:
        d=prior.ask(c,s,cases,fun,variants);save('alternating-proof-'+name,d)
        results[name]=assess(predictions[name],d['answers'])
    save('alternating-proof-result',results)
    print('alternating fresh proof',{k:(v['selections'],len(v['failures']))for k,v in results.items()},flush=True)


def seasonal_programs():
    from .ets_seasonal import initial, moving_initial, fit
    rng = random.Random(480202); cases = []
    for i in range(96):
        p = [2, 3, 4][i % 3]; n = rng.choice(list(range(max(8, 3*p), 33)))
        exponent = [-300, 0, 300][i % 3]
        pattern = [rng.randrange(-32, 33)/4 for _ in range(p)]
        y = [256.+128*j+pattern[j % p] for j in range(n)]
        if i % 2: y = [v+rng.randrange(-8, 9)/32 for v in y]
        else: y[rng.randrange(n)] += 1/32
        cases.append(dict(y=[math.ldexp(v, exponent) for v in y], season=p, kind='positive_trend_innovations'))
    freeze('seasonal-programs', dict(cases=cases,
        initializers=['first-season / first-two-cycle trend', 'moving2p centered', 'moving2p raw adjustment', 'rolling2p centered'],
        trend_programs=['raw OLS', 'centered covariance', 'resident PC64 centered covariance'],
        conditional='STAT parameters; a strong positive trend separates postfit trend reset. Forecasts frozen after parameter answers and before forecast answers.'))
    c, s = service(); parameters = prior.ask(c, s, cases, 'FORECAST.ETS.STAT', [1, 2, 3]); save('seasonal-programs-stat', parameters)
    programs = {}
    for initializer, trend_program in itertools.product(['firsttwo', 'moving', 'moving_raw', 'rolling'], ['raw', 'centered', 'centered64']):
        pred = []
        for i, z in enumerate(cases):
            y, p = z['y'], z['season']; a, b, g = [value(parameters['answers'][3*i+j]['bits']) for j in range(3)]
            init = initial(y, p, 'firsttwo', 'first_detrend') if initializer == 'firsttwo' else moving_initial(
                y, p, ma_order='subtract_add' if initializer=='rolling' else 'direct',
                adjustment='raw' if initializer=='moving_raw' else 'centered')
            if trend_program != 'raw':
                dtype = np.longdouble if trend_program=='centered64' else float
                adjusted = [dtype(v)-dtype(init[2][j % p]) for j, v in enumerate(y)]
                mean = sum(adjusted, dtype(0))/len(y); mx=dtype((len(y)-1)/2)
                numerator = sum(((dtype(j)-mx)*(v-mean) for j, v in enumerate(adjusted)), dtype(0))
                denominator = sum(((dtype(j)-mx)**2 for j in range(len(y))), dtype(0))
                init = init[0], float(numerator/denominator), init[2]
            l, t, se, _, _ = fit(y, p, a, b, g, init)
            pred.extend(bits(float(np.longdouble(l+h*t)+np.longdouble(se[(len(y)-1+int(h)) % p]))) for h in [0., 1., 2., 5.])
        programs[initializer+':'+trend_program] = pred
    freeze('seasonal-programs-conditional', dict(cases=cases, parameter_answers=parameters,
                                                predictions=programs, recurrence='Existing measured strict innovations; PC64 final addition.'))
    d = prior.ask(c, s, cases, 'FORECAST.ETS', [0., 1., 2., 5.]); save('seasonal-programs-forecast', d)
    results = {k: assess(p, d['answers']) for k, p in programs.items()}
    save('seasonal-programs-result', dict(programs=results, survivors=[k for k, r in results.items() if not r['failures']]))
    print('seasonal initializers', len(results), 'survivors', [k for k, r in results.items() if not r['failures']], flush=True)


def cycle_pairs():
    rng = random.Random(480401); cases=[]
    for i in range(128):
        exponent = [-400, -100, 0, 100, 400][i % 5]
        if i < 64: y = [math.ldexp(rng.randrange(-256, 257)/8, exponent) for _ in range(4)]
        else:
            center = int(bits(math.ldexp(rng.choice([1., 1/3]), exponent)), 16)
            y = [value(f'{center+rng.randrange(-16, 17):016x}') for _ in range(4)]
        cases.extend([dict(y=y, season=season, pair=i) for season in [0, 2]])
    freeze('cycle-pairs', dict(cases=cases,
        prediction='Every consecutive period0/period2 pair agrees exactly on every selection; no conditional numerical formula is assumed.',
        design='Four observations, dyadic and one-bit-stepped histories, including the input magnitude edges.'))
    c, s = service(); results={}
    for name, fun, variants in [('forecast', 'FORECAST.ETS', [0., 1., 2., 5.]),
                               ('stat', 'FORECAST.ETS.STAT', list(range(1, 9))),
                               ('confidence', 'FORECAST.ETS.CONFINT', [(1., .5), (1., .95), (5., .95)])]:
        d = prior.ask(c, s, cases, fun, variants); save('cycle-pairs-'+name, d)
        stride=len(variants); failures=[]
        for i in range(128):
            left=d['answers'][2*i*stride:(2*i+1)*stride]; right=d['answers'][(2*i+1)*stride:(2*i+2)*stride]
            if left != right: failures.append(dict(pair=i, left=left, right=right))
        results[name]=dict(selections=len(d['answers']), paired_histories=128, failures=failures)
    save('cycle-pairs-result', results)
    print('cycle pairs', {k:len(v['failures']) for k,v in results.items()}, flush=True)


def exact_period_map():
    """One fixed extension of exact repetition detection; stop on refutation."""
    rng=random.Random(480103);cases=[];pred=[]
    for p in range(3,9):
        for i in range(48):
            n=p*rng.randrange(3,32//p+1);exponent=[-300,0,300][i%3]
            pattern=[rng.randrange(-64,65)/8 for _ in range(p)]
            trend=[0.,.125,16.,65536.][i%4];offset=rng.randrange(-512,513)/4
            y=[math.ldexp(offset+pattern[j%p]+trend*j,exponent)for j in range(n)]
            period=0
            for q in range(2,min(8,n//3)+1):
                delta=y[q]-y[0]
                if all(y[j+q]-y[j]==delta for j in range(n-q)):
                    period=q;break
            cases.append(dict(y=y,season=1,constructed_period=p,kind='exact_affine_repetition'))
            pred.append(bits(float(period)))
    freeze('exact-period-map',dict(cases=cases,predictions=pred,
        method='Smallest exact repeating affine phase program with at least three complete cycles. This is one fixed larger-period extension of the existing exact p2 region, not a fitted cutoff or peak correction.'))
    c,s=service();d=prior.ask(c,s,cases,'FORECAST.ETS.SEASONALITY',[0]);save('exact-period-map-answers',d)
    result=assess(pred,d['answers']);save('exact-period-map-result',result)
    print('exact repetition detector',result['selections'],'contradictions',len(result['failures']),flush=True)


if __name__ == '__main__':
    steps = {'search': search_shapes, 'detection': detection_shapes, 'short-resident': short_resident,
             'short-centered': lambda: short_resident(True), 'short-seasons': short_seasons,
             'seasonal': seasonal_programs, 'cycles': cycle_pairs, 'alternating-proof': alternating_proof,
             'exact-period-map': exact_period_map}
    parser = argparse.ArgumentParser(); parser.add_argument('step', choices=steps)
    steps[parser.parse_args().step]()
