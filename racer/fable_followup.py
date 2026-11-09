"""Fable-1 follow-up, round49. Original MIT, implementer, tier1 only.

Sources: retained Excel observations, IEEE754 and elementary mathematics.
No restricted implementation source or excluded checkout is read. Uses the
existing audited cache and retains predictions before fresh requests.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path
import random
from datetime import date

import numpy as np
import mpmath as mp

from . import round41
from .fingerprint import bits, value, sha

OUT = Path('/home/ubuntu/codex-runs/round49')
INPUT = OUT
CACHE = None
ROOT = Path(__file__).resolve().parents[1]


def save(name, obj):
    with (OUT / (name + '.json')).open('x') as stream:
        json.dump(obj, stream, separators=(',', ':'), allow_nan=False)


def freeze(name, obj):
    save(name + '-frozen', dict(obj, source_sha256=sha(__file__),
                               source=Path(__file__).read_text()))


def service(offline=False):
    round41.OUT = OUT
    c,s=round41.service(offline=offline, priority='interactive')
    if CACHE is not None:s.directory=Path(CACHE)
    return c,s


def windows(centers, count):
    return [dict(center=bits(c), start=f'{int(bits(c),16)-count//2:016x}', count=count)
            for c in centers]


def window_inputs(ws):
    return [f'{int(w["start"],16)+i:016x}' for w in ws for i in range(w['count'])]


def gamma_map(offline=False):
    # Forty uniform centers plus exact joins and fixed candidate binade edges.
    mp.mp.dps = 65
    shifts = [5., 7., 607/128, 9.]
    centers = list(np.linspace(.70001, 2.99999, 40)) + [1., 1.5, 2.]
    edges = []
    for g in shifts:
        for power in [1., 2., 4., 8.]:
            def term(x):
                z = x+g+.5
                return z-(x+.5)*math.log(z)-power
            if term(.7)*term(3.) < 0:
                lo, hi = .7, 3.
                for _ in range(60):
                    mid=(lo+hi)/2
                    if term(mid)>0: lo=mid
                    else: hi=mid
                edge=(lo+hi)/2
                edges.append(dict(g=g,power=power,x=edge))
                centers += [edge-.00001,edge+.00001]
    ws=windows(sorted(set(float(c) for c in centers)),256)
    xs=window_inputs(ws)
    freeze('gamma-map',dict(windows=ws,edges=edges,
        hypothesis='Lanczos temporary binade crossings predict an error-sd doubling; no coefficients fitted.'))
    c,s=service(offline)
    ans=s.ask_rows('GAMMALN',[[c.number(x)] for x in xs])
    save('gamma-map-answers',dict(windows=ws,inputs=xs,answers=ans,receipts=s.receipts))
    summaries=[];offset=0
    for w in ws:
        raw=xs[offset:offset+w['count']];ys=ans[offset:offset+w['count']]
        if any(a.get('type')!='number' for a in ys):
            raise ValueError('unread/non-numeric core observation')
        errors=[float((mp.mpf(value(a['bits']))-mp.loggamma(mp.mpf(value(x))))*2**52)
                for x,a in zip(raw,ys)]
        summaries.append(dict(center=value(w['center']),sd=float(np.std(errors)),
                              mean=float(np.mean(errors)),low=min(errors),high=max(errors)))
        offset+=w['count']
    save('gamma-map-result',dict(units='2^-52 absolute',windows=summaries,edges=edges))
    print(json.dumps(summaries),flush=True)


def gamma_reduction(offline=False):
    rng=random.Random(49011)
    xs=[rng.uniform(3.,8.) for _ in range(1024)]
    xs += [value(f'{int(bits(b),16)+i:016x}') for b in [3.,4.,5.,6.,7.,8.]
           for i in range(-32,33) if 3.<=value(f'{int(bits(b),16)+i:016x}')<8.]
    lows=[value(f'{int(bits(.7),16)+i:016x}') for i in range(-128,129)]
    lows += [math.ldexp(rng.uniform(1.,2.),e) for e in range(-1021,0,7)]
    freeze('gamma-reduction',dict(inputs=[bits(x) for x in xs],lower=[bits(x) for x in lows],
        hypotheses=['ascending stored product; stored PC64 final sum','plain binary64 final sum',
                    'literal .7 join; stored PC64 final subtraction','plain subtraction']))
    c,s=service(offline)
    def ask(fn,x):
        out=s.ask_rows(fn,[[c.number(v)] for v in x])
        if any(a.get('type')!='number' for a in out): raise ValueError('non-numeric leaf')
        return [value(a['bits']) for a in out]
    ps=[]
    for x in xs:
        k=int(x)-2;p=x-k
        for j in range(1,k): p=float(np.longdouble(p)*np.longdouble(x-k+j))
        ps.append(p)
    target=ask('GAMMALN',xs);core=ask('GAMMALN',[x-(int(x)-2) for x in xs]);logs=ask('LN',ps)
    up=ask('GAMMALN',[x+1 for x in lows]);lower=ask('GAMMALN',lows);ll=ask('LN',lows)
    rows=[dict(x=bits(x),answer=bits(a),core=bits(b),log=bits(l),
               plain=bits(b+l),stored=bits(float(np.longdouble(b)+np.longdouble(l))))
          for x,a,b,l in zip(xs,target,core,logs)]
    lr=[dict(x=bits(x),answer=bits(a),upper=bits(u),log=bits(l),
             plain=bits(u-l),stored=bits(float(np.longdouble(u)-np.longdouble(l))))
        for x,a,u,l in zip(lows,lower,up,ll)]
    save('gamma-reduction-result',dict(rows=rows,lower=lr,receipts=s.receipts))
    print('upper',len(rows),{k:sum(r[k]==r['answer'] for r in rows) for k in ['plain','stored']})
    for flag in [True,False]:
        rs=[r for r in lr if (value(r['x'])<.7)==flag]
        print('below .7',flag,len(rs),{k:sum(r[k]==r['answer'] for r in rs) for k in ['plain','stored']})


def gamma_reference_audit():
    path=ROOT/'evidence/fable-1/e1c_gammaln_core_error.json'
    doc=json.loads(path.read_text());mp.mp.dps=80
    rows=[]
    for r in doc['rows']:
        x=mp.mpf(value(r['x']));truth=mp.loggamma(x);a=mp.mpf(value(r['gammaln']))
        rows.append(dict(x=r['x'],answer=r['gammaln'],truth_rn=bits(float(truth)),
                         error=float(a-truth),reference_error=float(mp.mpf(value(r['true']))-truth)))
    summary=[]
    for lo,hi in [(.7,1.),(1.,1.5),(1.5,2.),(2.,2.5),(2.5,3.)]:
        rs=[r for r in rows if lo<=value(r['x'])<hi]
        summary.append(dict(region=[lo,hi],n=len(rs),
            max_error=max(abs(r['error']) for r in rs),
            old_reference_max_error=max(abs(r['reference_error']) for r in rs)))
    result=dict(reference=doc['truth'],source=str(path),sha256=sha(path),rows=rows,summary=summary)
    save('gamma-reference-audit',result)
    print(json.dumps(summary),flush=True)


def lattice(start, answers):
    """The retained R43 phase-coherence estimator, in absolute units 2^-53.

    The largest coherent period below 2.5 sd avoids selecting a shorter
    harmonic. Low-coherence windows stay diagnostics, not exact constraints.
    """
    xs=np.array([value(a) for a in answers],dtype=float)
    differences=np.diff(xs/2.**-53)
    sd=float(np.std(differences));u=2. if xs[0]>=1. else 1.
    qs=np.arange(1.5,40.,.01)*u
    scores=np.abs(np.exp(2j*np.pi*np.outer(1/qs,differences)).mean(axis=1))
    ok=qs<2.5*sd
    if not ok.any(): return dict(sd=sd,period=None,coherence=0.)
    maximum=float(scores[ok].max())
    candidates=np.flatnonzero(ok & (scores>=.6*maximum))
    i=int(candidates[-1]);lo=max(0,i-30);hi=min(len(qs),i+30)
    j=lo+int(np.argmax(scores[lo:hi]));period=float(qs[j])
    fine=np.linspace(max(1.5*u,period-.015*u),period+.015*u,121)
    fs=np.abs(np.exp(2j*np.pi*np.outer(1/fine,differences)).mean(axis=1))
    j=int(np.argmax(fs))
    return dict(sd=sd,period=float(fine[j]),coherence=float(fs[j]))


def normal_fit():
    path=INPUT/'retained-normal-windows.json'
    doc=json.loads(path.read_text());records=[]
    for w in doc['windows']:
        if any(v is None for v in w['rows']): raise ValueError('missing retained answer')
        estimate=lattice(w['start'],w['rows'])
        records.append(dict(t=w['t'],s=w['t']**2,**estimate))
    selected=[r for r in records if r['t']>=.52 and r['coherence']>=.35]
    ss=np.array([r['s'] for r in selected]);ll=np.array([r['period'] for r in selected])
    # Test both the proposed absolute spacing and the spacing relative to x.
    mp.mp.dps=50
    yy=np.array([float(mp.sqrt(2)*mp.erfinv(r['t'])) for r in selected])
    models=[]
    for label,ordinates in [('inverse_absolute',1/ll),('inverse_relative',yy/ll)]:
        line=np.polyfit(ss,ordinates,1);fit=np.polyval(line,ss)
        root=float(-line[1]/line[0]);prediction=[]
        for point in [.46,.48,.495]:
            norm=1. if label=='inverse_absolute' else float(mp.sqrt(2)*mp.erfinv(mp.sqrt(point)))
            prediction.append(norm/float(np.polyval(line,point)))
        models.append(dict(label=label,line=line.tolist(),zero=root,
            relative_rms=float(np.sqrt(np.mean(((fit-ordinates)/ordinates)**2))),
            predicted_period=prediction))
    freeze('normal-denominator',dict(retained=str(path),sha256=sha(path),records=records,
        selected=len(selected),models=models,probe_s=[.46,.48,.495],count=2048,
        controls_t=[.75],limits='Lattice estimates have phase noise; denominator scale and C are not identified.'))
    print(json.dumps(models),flush=True)


def normal_probe(offline=False):
    frozen=json.loads((INPUT/'normal-denominator-frozen.json').read_text())
    centers=[.5+math.sqrt(s)/2 for s in frozen['probe_s']]+[.875]
    ws=windows(centers,2048);xs=window_inputs(ws)
    freeze('normal-probe',dict(windows=ws,predictions=frozen['models'],
        note='Frozen offline denominator fits; fourth window t=.75 is piece B.'))
    c,s=service(offline)
    ans=s.ask_rows('NORM.S.INV',[[c.number(x)] for x in xs])
    results=[]
    for i,w in enumerate(ws):
        rows=ans[i*2048:(i+1)*2048]
        if any(a.get('type')!='number' for a in rows): raise ValueError('non-numeric normal observation')
        results.append(dict(center=value(w['center']),**lattice(w['start'],[a['bits'] for a in rows])))
    save('normal-probe-result',dict(windows=ws,inputs=xs,answers=ans,results=results,receipts=s.receipts))
    print(json.dumps(results),flush=True)


def bond_date(months):
    year=2020+months//12;month=1+months%12
    return float((date(year,month,1)-date(1899,12,30)).days)


def yield_scan(offline=False):
    # Exact whole coupon counts. Equal issue/settlement odd controls below
    # are retained as typed NUM results; odd_scan supplies the valid scan.
    designs=[(2,1,25.),(2,1,50.),(2,1,80.),(4,1,6.25),(8,1,100/256),(4,2,6.25)]
    cases=[]
    for n,f,price in designs:
        for exponent in range(-500,501):
            args=[bond_date(0),bond_date(12*n//f),0.,math.ldexp(price,exponent),
                  math.ldexp(100.,exponent),float(f),0.]
            odd=[args[0],args[1],args[0],bond_date(12//f),*args[2:]]
            cases.append(dict(n=n,f=f,price=price,exponent=exponent,
                              regular=[bits(v) for v in args],odd=[bits(v) for v in odd]))
    freeze('yield-scan',dict(cases=cases,seed_bounds='finite limit containing face100, not fitted',
        inverse='z_next=((n+1)z-r*z^(n+1))/n, z=1+y/f, r=price/redemption',
        warning='Inverse is cubic for n=2; up to three real preimages, at most two positive.'))
    c,s=service(offline)
    regular=s.ask_rows('YIELD',[[c.number(x) for x in z['regular']] for z in cases])
    odd=s.ask_rows('ODDFYIELD',[[c.number(x) for x in z['odd']] for z in cases])
    summaries=[]
    for j,(n,f,price) in enumerate(designs):
        rr=regular[j*1001:(j+1)*1001];oo=odd[j*1001:(j+1)*1001]
        runs=[]
        for a,b in zip(rr,oo):
            runs.append(a==b)
        summaries.append(dict(n=n,f=f,price=price,regular_small=rr[0],regular_large=rr[-1],
            odd_small=oo[0],odd_large=oo[-1],same=sum(runs),selections=1001,
            regular_numeric=sum(a.get('type')=='number' for a in rr),
            odd_numeric=sum(a.get('type')=='number' for a in oo),
            odd_distinct=len({json.dumps(a,sort_keys=True) for a in oo})))
    save('yield-scan-result',dict(cases=cases,regular=regular,odd=odd,summary=summaries,receipts=s.receipts))
    print(json.dumps(summaries),flush=True)


def inverse_newton_interval(interval, ratio, n):
    """Real preimages on all three monotone branches of even-n Newton.

    Decimal arithmetic only locates the conditional analytic map; endpoints
    are widened in the caller. This does not bound a private price kernel.
    """
    r=mp.mpf(ratio);root=r**(-mp.mpf(1)/n)
    def fun(z): return ((n+1)*z-r*z**(n+1))/n
    def positive_outer(w):
        # F maps [root,infinity) onto (-infinity,root].
        if w>root: raise ValueError('outside positive outer branch')
        lo=root;hi=2*root+abs(w)+1
        for _ in range(64):
            if fun(hi)<=w: break
            hi*=2
        else: raise ArithmeticError('failed to bracket inverse Newton branch')
        for _ in range(100):
            mid=(lo+hi)/2
            if fun(mid)>w:lo=mid
            else:hi=mid
        return lo,hi
    def central(w):
        lo=-root;hi=root
        for _ in range(100):
            mid=(lo+hi)/2
            if fun(mid)<w:lo=mid
            else:hi=mid
        return lo,hi
    low,high=interval;out=[]
    if low<=root:
        a=positive_outer(min(high,root))[0];b=positive_outer(low)[1];out.append((a,b))
    if high>=-root:
        a=-positive_outer(-high)[1];b=-positive_outer(min(-low,root))[0];out.append((a,b))
    a,b=max(low,-root),min(high,root)
    if a<=b:out.append((central(a)[0],central(b)[1]))
    return out


def merge_intervals(intervals):
    merged=[]
    for a,b in sorted(intervals):
        if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
        else:merged.append((a,b))
    return merged


def yield_backward():
    doc=json.loads((INPUT/'yield-scan-result.json').read_text());mp.mp.dps=60
    freeze('yield-backward-v3',dict(max_depth=24,output_padding_ulps=64,
        inverse='All three real monotone branches, endpoint bisection brackets; no restriction to positive predecessors.'))
    records=[]
    for q in doc['summary']:
        n,f,price=q['n'],q['f'],q['price'];ratio=mp.mpf(price)/100
        for side in ['small','large']:
            answer=q['regular_'+side]
            if answer.get('type')!='number':continue
            y=value(answer['bits']);z=1+mp.mpf(y)/f
            # Conditional envelope includes output rounding and 64 output ulps.
            # It is deliberately not presented as an unknown-kernel error bound.
            pad=mp.mpf(math.ulp(y))*64/f
            intervals=[(z-pad,z+pad)];steps=[]
            probes={'two_f_over_n':mp.mpf(2*f)/n,'f_over_n':mp.mpf(f)/n,
                    'zero':mp.mpf(0),'point_one':mp.mpf('0.1'),
                    'negative_two_f_over_n':-mp.mpf(2*f)/n,
                    'cash_limit':mp.mpf(2*f)/n*(1/ratio-1)}
            for depth in range(25):
                admitted=[name for name,seed in probes.items()
                          if any(a<=1+seed/f<=b for a,b in intervals)]
                clipped=[(max(a,mp.mpf(0)),min(b,mp.mpf(3))) for a,b in intervals if b>=0 and a<=3]
                steps.append(dict(depth=depth,intervals=len(intervals),seed_names=admitted,
                    bounded_z_0_3=[[str(a),str(b)] for a,b in clipped]))
                if depth<24:
                    intervals=merge_intervals([v for interval in intervals
                        for v in inverse_newton_interval(interval,ratio,n)])
            records.append(dict(n=n,f=f,price=price,side=side,answer=answer,steps=steps))
            print(n,f,price,side,[(z['depth'],z['intervals'],z['seed_names']) for z in steps],flush=True)
    save('yield-backward-v3-result',dict(records=records,source=Path(__file__).read_text(),
        limits='Conditional exact analytic Newton in yield; 64-ULP output envelope is a sensitivity control, not a bound on private price or slope arithmetic.'))


def odd_scan(offline=False):
    doc=json.loads((INPUT/'yield-scan-result.json').read_text());cases=[]
    for q in doc['cases']:
        row=q['odd'].copy();row[2]=bits(bond_date(-12//q['f']))
        cases.append(dict(**{k:q[k] for k in ['n','f','price','exponent']},args=row))
    freeze('odd-scan-v2',dict(cases=cases,
        design='Strict issue < settlement < firstcoupon. Whole remaining periods; zero coupon removes the long first coupon payment. Initial equal-issue requests were all NUM and are retained as invalid controls.'))
    c,s=service(offline)
    answers=s.ask_rows('ODDFYIELD',[[c.number(x) for x in q['args']] for q in cases])
    results=[]
    for j,q in enumerate(doc['summary']):
        rr=doc['regular'][j*1001:(j+1)*1001];oo=answers[j*1001:(j+1)*1001]
        results.append(dict(n=q['n'],f=q['f'],price=q['price'],small=oo[0],large=oo[-1],
            numeric=sum(a.get('type')=='number' for a in oo),same_as_yield=sum(a==b for a,b in zip(rr,oo)),
            distinct=len({json.dumps(a,sort_keys=True) for a in oo})))
    save('odd-scan-v2-result',dict(cases=cases,answers=answers,summary=results,receipts=s.receipts))
    print(json.dumps(results),flush=True)


def regular_first_control(offline=False):
    cases=[]
    for f in [1,2,4]:
        for exponent in range(-500,501,10):
            # First period is regular; settlement is strictly inside it.
            settle=bond_date(max(1,6//f));maturity=bond_date(24)
            rate=math.ldexp(.05,exponent);price=math.ldexp(80.,exponent)
            redemption=math.ldexp(100.,exponent)
            regular=[settle,maturity,rate,price,redemption,float(f),0.]
            odd=[settle,maturity,bond_date(0),bond_date(12//f),*regular[2:]]
            cases.append(dict(f=f,exponent=exponent,regular=[bits(v) for v in regular],
                              odd=[bits(v) for v in odd]))
    freeze('regular-first-control',dict(cases=cases,
        hypothesis='With a regular first period and scaled coupon, YIELD and ODDFYIELD have the same price function. Equality of their answers tests shared solver behaviour.'))
    c,s=service(offline)
    rr=s.ask_rows('YIELD',[[c.number(v) for v in q['regular']] for q in cases],chunk_size=100)
    oo=s.ask_rows('ODDFYIELD',[[c.number(v) for v in q['odd']] for q in cases],chunk_size=100)
    save('regular-first-control-result',dict(cases=cases,regular=rr,odd=oo,
        equal=sum(a==b for a,b in zip(rr,oo)),
        regular_numeric=sum(a.get('type')=='number' for a in rr),
        odd_numeric=sum(a.get('type')=='number' for a in oo),receipts=s.receipts))
    print('regular first',len(cases),'equal',sum(a==b for a,b in zip(rr,oo)),flush=True)


def odd_calendar_controls(offline=False):
    rows=[];labels=[]
    for f in [1,2,4]:
        for rate in [0.,.05]:
            for offset in [-1.,0.,1.]:
                rows.append([bond_date(max(1,6//f)),bond_date(24),bond_date(0)+offset,
                             bond_date(12//f),rate,80.,100.,float(f),0.])
                labels.append(dict(f=f,rate=rate,issue_offset_days=offset))
    rows.append([39763.,44256.,39736.,39873.,.0575,84.5,100.,2.,0.])
    labels.append(dict(control='ld1 corpus bond'))
    freeze('odd-calendar',dict(rows=[[bits(v) for v in q] for q in rows],labels=labels))
    c,s=service(offline)
    ans=s.ask_rows('ODDFYIELD',[[c.number(v) for v in q] for q in rows])
    price_rows=[row[:5]+[.1]+row[6:] for row in rows]
    prices=s.ask_rows('ODDFPRICE',[[c.number(v) for v in q] for q in price_rows])
    save('odd-calendar-result',dict(rows=rows,labels=labels,answers=ans,prices=prices,receipts=s.receipts))
    print(json.dumps(list(zip(labels,ans,prices))),flush=True)


def corpus_function_counts():
    """Attribute existing verdicts with the existing corpus formula labeller."""
    from .corpus_cases import workbook_directory, formula_label
    summary=json.loads((INPUT/'summary.json').read_text())
    counts=defaultdict(Counter);total=Counter()
    for workbook in summary['workbooks']:
        total['awaiting_capture']+=workbook['counts']['awaiting_capture']
        directory=workbook_directory(INPUT,workbook,summary.get('version',1))
        captures=sorted(directory.glob('*.capture.json'))
        runs=[(p,directory/p.name.removesuffix('.capture.json')) for p in captures]
        if not runs and workbook['capture'].startswith('committed capture '):
            runs=[(Path(workbook['capture'].removeprefix('committed capture ')),directory)]
        if not runs:
            if workbook['counts']['awaiting_capture']!=workbook['counts']['formula_cells']:
                raise ValueError('unaccounted workbook')
        for capture_path,run in runs:
            cap=json.loads(capture_path.read_text())
            key=lambda row:(row['sheet'],row['row'],row['column'])
            formulas={key(row):row for row in cap['formulas']}
            triage=json.loads((run/'triage.json').read_text())
            for verdict in triage['cells']:
                k=key(verdict)
                label=formula_label(formulas.get(k,{}).get('formula',''),formulas,k)
                result=verdict['verdict']
                if result in ['matched','agreed_rejection']:
                    result=('circular_rule' if verdict.get('circular_rule') else
                            'external_links' if verdict.get('external_links') else 'matched')
                counts[label][result]+=1;total[result]+=1
    headline=summary['totals']['counts']
    if total['matched']!=headline['pass'] or sum(total.values())!=headline['formula_cells']:
        raise ValueError(f'function attribution does not reconcile with headline: {dict(total)}')
    save('corpus-function-counts',dict(totals=dict(total),functions={k:dict(v) for k,v in counts.items()},
        source_sha256=sha(INPUT/'summary.json'),
        rule='Existing formula_label, one formula cell per label, existing verdicts; external/circular/pending cells remain failures. Nested calls are not multiplied.'))
    print('corpus attribution',dict(total),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    commands={'gamma-map':gamma_map,'gamma-reduction':gamma_reduction,
              'normal-probe':normal_probe,'yield-scan':yield_scan,'odd-scan':odd_scan,
              'regular-first-control':regular_first_control,'odd-calendar':odd_calendar_controls}
    local={'gamma-reference-audit':gamma_reference_audit,'normal-fit':normal_fit,
           'yield-backward':yield_backward,'corpus-counts':corpus_function_counts}
    parser.add_argument('command',choices=list(commands)+list(local))
    parser.add_argument('--offline',action='store_true')
    parser.add_argument('--out',type=Path,default=OUT)
    parser.add_argument('--inputs',type=Path,default=INPUT)
    parser.add_argument('--cache',type=Path)
    args=parser.parse_args()
    OUT=args.out;INPUT=args.inputs;CACHE=args.cache
    OUT.mkdir(parents=True,exist_ok=True)
    if args.command in commands: commands[args.command](args.offline)
    else: local[args.command]()
