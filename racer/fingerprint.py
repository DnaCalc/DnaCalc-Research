"""Consecutive-double Excel fingerprints, using the existing scalar service.

Original MIT tooling. Inputs, typed answers and worker provenance are retained
in small gzip files. Spectral peaks are observations, not a count of hidden
stores: only a discriminated, fully matching computation tree can supply that.
Probe the mapping, single-bit steps and exponent edges to find the shape. Recover
its exact constants by bisection or fitting all known samples, then prove on fresh
answers. Search for inputs separating adjacent constant values; otherwise report
them as undetermined. Growing degrees/pieces or failing fresh answers calls the
shape into question. Rank trees by plausible hand-written 1990s code, not match
percentage; exact answer comparisons remain exact.
"""
import argparse
from collections import Counter
import cmath
from fractions import Fraction as Q
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'checker'))
from exact import BINARY64, EXTENDED64
from reference import elementary_bounds
try:
    from .x87_programs import Builder,Native,build_library
except ImportError:
    from x87_programs import Builder,Native,build_library

PRODUCER_BINDINGS={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in
                   (Path(__file__),ROOT/'checker/reference.py',ROOT/'checker/exact.py')}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bits(value):
    return struct.pack('>d', float(value)).hex()


def value(raw):
    return struct.unpack('>d', bytes.fromhex(raw.removeprefix('0x')))[0]


def write(path, document):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(document,dict)and str(document.get('format','')).startswith('dnacalc-'):
        bindings=dict(PRODUCER_BINDINGS)
        for source,digest in document.get('source_bindings',{}).items():
            if source in bindings and bindings[source]!=digest:
                raise ValueError('conflicting analysis source binding: '+source)
            bindings[source]=digest
        if bindings!={p:sha(p)for p in bindings}:
            raise ValueError('analysis source changed; raw replies remain retained, regenerate this analysis')
        document=dict(document,producer_bindings=bindings,
                      analysis_model='Codex session; exact model identifier unavailable',
                      sources='Tier1 original Excel scalar outputs and IEEE754; original MIT arithmetic/DFT and installed MPFR interface. No tier3/4 implementation source read.')
    temporary=path.with_suffix(path.suffix+'.'+uuid.uuid4().hex+'.tmp')
    with temporary.open('x')as stream:
        stream.write(json.dumps(document, indent=2, allow_nan=False)+'\n')
    os.replace(temporary,path)


def load_service(checkout):
    sys.path.insert(0, str(Path(checkout).resolve() / 'tools/excel-service'))
    import client
    import policy
    return client, policy


def typed(argument, client):
    if isinstance(argument, bool):
        return client.boolean(argument)
    if isinstance(argument, (float, int)):
        return client.number(float(argument))
    if isinstance(argument, str) and re.fullmatch(r'0x[0-9a-f]{16}', argument):
        return client.number(argument)
    if isinstance(argument, str) and not re.search(r'[\\/:\[\]!|\x00-\x1f]', argument):
        return client.text(argument)
    raise ValueError('non-scalar or service-inadmissible seed argument')


# Fixed ordinary examples and canonical textbook boundaries, chosen before
# reading any R26 outputs. Other seeds use the existing pins, without reading
# their current predictions or sorting on whether DNA Calc matches Excel.
SEEDS = {
    'EXP': ([0.7], 0), 'LN': ([1.3], 0), 'LOG': ([1.3, 10.0], 0),
    'LOG10': ([1.3], 0), 'POWER': ([1.3, 1.0/3.0], 0),
    'SQRT': ([1.3], 0), 'GEOMEAN': ([1.3, 2.1, 3.7], 0),
    'BINOM.DIST': ([3, 20, 0.3, False], 2),
    'BINOMDIST': ([3, 20, 0.3, False], 2),
    'BINOM.DIST.RANGE': ([20, 0.3, 2, 5], 1),
    'POISSON.DIST': ([3, 2.3, False], 1), 'POISSON': ([3, 2.3, False], 1),
    'NEGBINOM.DIST': ([3, 5, 0.3, False], 2),
    'NEGBINOMDIST': ([3, 5, 0.3], 2),
    'BESSELJ': ([2.3, 0], 0), 'BESSELY': ([2.3, 0], 0),
    'BETA.DIST': ([0.3, 2, 3, True], 0), 'BETADIST': ([0.3, 2, 3], 0),
    'GAMMA.DIST': ([0.7, 2, 1, True], 0), 'GAMMADIST': ([0.7, 2, 1, True], 0),
    'GAMMA': ([1.3], 0), 'GAMMALN': ([1.3], 0), 'GAMMALN.PRECISE': ([1.3], 0),
    'NORM.S.DIST': ([0.7, True], 0), 'NORM.DIST': ([0.7, 0, 1, True], 0),
    'NORMDIST': ([0.7, 0, 1, True], 0), 'ERF': ([0.7], 0),
    'ERF.PRECISE': ([0.7], 0), 'ERFC': ([0.7], 0), 'ERFC.PRECISE': ([0.7], 0),
    'GAUSS': ([0.7], 0), 'LOGNORM.DIST': ([1.3, 0, 1, True], 0),
    'LOGNORMDIST': ([1.3, 0, 1], 0), 'KURT': ([1.3, 2.1, 4.7, 6.2, 8.1], 0),
    'SKEW': ([1.3, 2.1, 4.7, 6.2, 8.1], 0),
    'TEXT': ([1.3, '0.0000000000000000'], 0),
    'FIXED': ([1.3, 14, True], 0), 'JIS': (['abc123'], None),
    'TEXTSPLIT': (['a,b', ','], None),
    'CONCAT': ([1.3, 'abc'], 0), 'CONCATENATE': ([1.3, 'abc'], 0),
    'SWITCH': ([1.3, 1, 2, 3], 0), 'COMPLEX': ([1.3, 2.1], 0),
    'IMDIV': ([1.3, 2.1], 0),
    'MEDIAN': ([1.3, 2.1, 3.7], 0), 'RRI': ([10, 1.3, 2.1], 1),
    'STANDARDIZE': ([1.3, 0.7, 2.1], 0),
    'HYPGEOM.DIST': ([3, 10, 12, 30, False], 0),
    'HYPGEOMDIST': ([3, 10, 12, 30], 0),
    'T.TEST': ([1.3, 2.1, 2, 1], 0), 'TTEST': ([1.3, 2.1, 2, 1], 0),
    'F.TEST': ([1.3, 2.1], 0), 'Z.TEST': ([1.3, 0.7, 2.1], 1),
    'ZTEST': ([1.3, 0.7, 2.1], 1), 'CHISQ.TEST': ([1.3, 2.1], 0),
    'IRR': ([1.3, 0.1], 0),
    'NPV': ([0.13, -100, 35, 40, 50], 0),
    'CUMIPMT': ([0.013, 12, 1000, 1, 6, 0], 0),
    'ODDFYIELD': ([43831, 45657, 43500, 44197, 0.05, 97, 100, 2, 0], 5),
    'CONFIDENCE': ([0.05, 1.3, 20], 1),
    'CONFIDENCE.NORM': ([0.05, 1.3, 20], 1),
    'CONFIDENCE.T': ([0.05, 1.3, 20], 1),
    'GAMMAINV': ([0.3, 2, 1], 0), 'GAMMA.INV': ([0.3, 2, 1], 0),
    'LOGINV': ([0.3, 0, 1], 0), 'LOGNORM.INV': ([0.3, 0, 1], 0),
    'BETA.INV': ([0.3, 2, 3], 0), 'BETAINV': ([0.3, 2, 3], 0),
    'F.INV': ([0.3, 4, 8], 0), 'F.INV.RT': ([0.3, 4, 8], 0),
    'T.DIST': ([0.7, 8, True], 0), 'T.DIST.2T': ([0.7, 8], 0),
    'T.INV.2T': ([0.3, 8], 0),
    'QUOTIENT': ([1.3,.7],0), 'EUROCONVERT': ([1.3,'DEM','EUR'],0),
}
for _name in ('SINH', 'COSH', 'TANH', 'COTH', 'CSCH', 'SECH', 'IMCOSH',
              'IMCOT', 'IMCSCH', 'IMLOG10', 'IMLOG2', 'IMSECH', 'IMSINH', 'IMTAN',
              'ABS', 'DEGREES', 'RADIANS', 'INT', 'N', 'NOT', 'EVEN', 'ODD',
              'ISEVEN', 'ISODD', 'SIGN'):
    SEEDS[_name] = ([1.3], 0)
for _name in ('ROUND', 'ROUNDDOWN', 'ROUNDUP', 'TRUNC'):
    SEEDS[_name] = ([1.375, 2], 0)
for _name in ('CEILING', 'CEILING.PRECISE', 'CEILING.MATH', 'ISO.CEILING',
              'FLOOR', 'FLOOR.PRECISE', 'FLOOR.MATH'):
    SEEDS[_name] = ([1.3, 0.1], 0)


def seed(function, checkout, client):
    if function in SEEDS:
        args, index = SEEDS[function]
        return [typed(a, client) for a in args], index, 'fixed textbook example'
    path = Path(checkout) / 'crates/calc-func/tests/witnesses' / (function + '.json')
    if not path.exists():
        return None, None, 'no scalar seed in the existing pin inventory'
    for row in json.loads(path.read_text())['witnesses']:
        try:
            args = [typed(a, client) for a in row['args']]
        except ValueError:
            continue
        indices = [i for i, a in enumerate(args) if a['type'] == 'number'
                   and 2**-20 < value(a['bits']) < 2**20]
        if indices and row['expected_bits'].startswith('0x'):
            return args, indices[0], 'first ordinary numeric pin: ' + row['id']
    return None, None, 'pins require arrays, references, or lack an ordinary positive numeric coordinate'


def fft(values):
    """Radix-two DFT; no external FFT dependency or resampling of input doubles."""
    n = len(values)
    if n == 1:
        return [complex(values[0])]
    if n == 0 or n & (n-1):
        raise ValueError('FFT length must be a power of two')
    a, b = fft(values[::2]), fft(values[1::2])
    turns = [cmath.exp(-2j*math.pi*k/n)*b[k] for k in range(n//2)]
    return [a[k]+turns[k] for k in range(n//2)] + [a[k]-turns[k] for k in range(n//2)]


def spectrum(answers):
    if not answers or any(a['type'] != 'number' or not math.isfinite(value(a['bits'])) for a in answers):
        return dict(status='non-numeric window', lines=[], resolved_stored_stages=None)
    exact = [Q(value(a['bits'])) for a in answers]
    unit = Q(math.ulp(value(answers[len(answers)//2]['bits'])))
    # Subtract in rational arithmetic before converting to double, so low bits
    # survive a large output offset. Least-squares affine detrending is fixed.
    y = [float((v-exact[0])/unit) for v in exact]
    n = len(y)
    xmean, ymean = (n-1)/2, sum(y)/n
    denominator = n*(n*n-1)/12
    slope = sum((i-xmean)*(v-ymean) for i, v in enumerate(y))/denominator
    residual = [v-ymean-slope*(i-xmean) for i, v in enumerate(y)]
    rms = math.sqrt(sum(v*v for v in residual)/n)
    if rms == 0:
        return dict(status='constant or exactly affine', rms_ulp=0, slope_ulp=slope,
                    lines=[], resolved_stored_stages=None)
    # Hann window suppresses leakage from the finite interval. A local maximum
    # above 16 times the median non-DC power is a reported line, never a store.
    window = [v*(0.5-0.5*math.cos(2*math.pi*i/(n-1))) for i, v in enumerate(residual)]
    powers = [abs(v)**2 for v in fft(window)][:n//2+1]
    median = sorted(powers[2:])[len(powers[2:])//2]
    peaks = [k for k in range(1, len(powers))
             if powers[k] > powers[k-1] and (k==len(powers)-1 or powers[k] >= powers[k+1])
             and powers[k] > 16*median]
    peaks.sort(key=lambda k: -powers[k])
    lines = [dict(bin=k, cycles_per_input=k/n, period_inputs=n/k,
                  power_fraction=powers[k]/sum(powers[1:])) for k in peaks[:8]]
    return dict(status='measured', rms_ulp=rms, slope_ulp=slope,
                median_power=median, lines=lines, detected_lines=len(peaks),
                resolved_stored_stages=None,
                limitation='Harmonics, aliasing and unobserved/redundant stores prevent converting peak count into store count.')


def window_residual(answers):
    """Affine residual in output ulps, preserving all binary64 low bits.

    Optional NumPy accelerates phase studies. This is an observation of the
    public output lattice, with no mathematical Gamma approximation subtracted.
    A window must stay inside one finite, nonzero output binade.
    """
    import numpy as np
    if any(a.get('type') != 'number' for a in answers):
        raise ValueError('numeric phase window required')
    ys = [value(a['bits']) for a in answers]
    if not ys or not all(math.isfinite(y) and y != 0 and (y > 0) == (ys[0] > 0)
                         and math.ulp(y) == math.ulp(ys[0]) for y in ys):
        raise ValueError('phase window crosses an output binade')
    # Integer lattice differences are exact even when the offset is large.
    origin = int(answers[0]['bits'], 16)
    y = np.array([int(a['bits'], 16)-origin for a in answers], dtype=float)
    if ys[0] < 0:
        y = -y
    k = np.arange(len(y))-(len(y)-1)/2
    y -= y.mean()
    y -= k * (k @ y)/(k @ k)
    return y


class PhaseTrack:
    """Read one spectral phase across a small, contiguous interval.

    The polynomial is a local phase model, in cycles, not a function kernel.
    Its constant is arbitrary and its linear term has integer/sign alias
    ambiguities. Decimal keeps fractional cycles when the phase reaches 2^52.
    The caller must inspect amplitude, innovations and residuals; a fit does
    not establish a stage, a store, or a coefficient table.
    """
    def __init__(self, name, origin, frequency, slope, degree=7):
        from decimal import Decimal as D, localcontext
        if degree < 2:
            raise ValueError('phase model needs room for the initial curvature')
        self.name, self.origin, self.degree = name, D(origin), degree
        with localcontext() as ctx:
            ctx.prec = 70
            dx = D(math.ulp(origin))
            self.coefficients = [D(0) for _ in range(degree+1)]
            self.coefficients[1] = D(frequency)/dx
            self.coefficients[2] = D(slope)/(2*dx)
        self.observations, self.misses = [], []
        self.active_degree = 1
        self.reference_phase = self.reference_amplitude = None

    def evaluate(self, x, derivative=0):
        from decimal import Decimal as D, localcontext
        with localcontext() as ctx:
            ctx.prec = 70
            w = D(x)-self.origin
            y = D(0)
            for j in range(len(self.coefficients)-1, derivative-1, -1):
                factor = math.prod(range(j-derivative+1, j+1))
                y = y*w+self.coefficients[j]*factor
            return y

    def demodulate(self, xs, residual):
        """Measure a window without changing the phase model or calibration."""
        import numpy as np
        from decimal import Decimal as D, localcontext
        center = len(xs)//2
        x = xs[center]
        theta = self.evaluate(x)
        frequency = float(self.evaluate(x, 1)*D(xs[center+1]-x))
        k = np.arange(len(xs))-center
        w = np.hanning(len(xs))
        with localcontext() as ctx:
            ctx.prec = 70
            phase = float(theta % 1)+(frequency % 1)*k
        z = np.sum(residual*w*np.exp(-2j*np.pi*phase))*2/w.sum()
        amplitude = float(abs(z))
        observed_phase = float(np.angle(z)/(2*np.pi))
        return x, theta, amplitude, observed_phase

    def observe(self, xs, residual):
        import numpy as np
        from decimal import Decimal as D, localcontext
        x, theta, amplitude, observed_phase = self.demodulate(xs,residual)
        if self.reference_phase is None:
            self.reference_phase, self.reference_amplitude = observed_phase, amplitude
        if amplitude < max(.025, .35*self.reference_amplitude):
            self.misses.append(dict(x=bits(x), amplitude=amplitude))
            return dict(name=self.name, status='weak', amplitude=amplitude)
        innovation = (observed_phase-self.reference_phase+.5) % 1-.5
        with localcontext() as ctx:
            ctx.prec = 70
            self.observations.append((x, theta+D(innovation), amplitude))
        # Do not fit high derivatives to a sub-resolution interval. In this
        # diagnostic continuation a unit-sized smooth Taylor coefficient must
        # span .2 cycles before it becomes a fitted coordinate. Short windows'
        # phase noise otherwise invents huge curvature and slips whole cycles.
        # This is an explicit continuation prior, never a coefficient proof.
        radius = max(abs(x-float(self.origin)) for x, _, _ in self.observations)
        dx = math.ulp(float(self.origin))
        ceiling = max([1]+[j for j in range(2,self.degree+1) if radius**j/dx > .2])
        degree = min(self.active_degree, (len(self.observations)-1)//2)
        errors = self.fit(degree)
        while degree < ceiling and len(self.observations) >= 2*(degree+2):
            saved = self.coefficients.copy()
            candidate = self.fit(degree+1)
            change = abs(float(self.coefficients[degree+1]-saved[degree+1]))
            if change <= 6*self.fit_standard_errors[degree+1]:
                self.coefficients = saved
                break
            degree += 1
            errors = candidate
        self.active_degree = max(1,degree)
        return dict(name=self.name, status='read', amplitude=amplitude,
                    innovation_cycles=innovation, degree=self.active_degree,
                    rms_cycles=float(np.std(errors)), max_cycles=float(np.max(np.abs(errors))))

    def fit(self, degree):
        import numpy as np
        from decimal import Decimal as D, localcontext
        if degree < 1:
            return np.zeros(len(self.observations))
        with localcontext() as ctx:
            ctx.prec = 70
            xs = [D(x)-self.origin for x, _, _ in self.observations]
            design = np.array([[float(x**j) if j else 1. for j in range(degree+1)] for x in xs])
            scales = np.max(np.abs(design), axis=0)
            scales[scales == 0] = 1
            weights = np.array([a for _, _, a in self.observations])
            errors = np.array([float(y-self.evaluate(x)) for x, y, _ in self.observations])
            correction, _, _, _ = np.linalg.lstsq(design/scales*weights[:, None], errors*weights, rcond=None)
            for j, c in enumerate(correction/scales):
                self.coefficients[j] += D(float(c))
            residual = np.array([float(y-self.evaluate(x)) for x, y, _ in self.observations])
            # A fixed weighted noise floor prevents exact-looking early fits
            # from admitting an unmeasurable high derivative. These errors are
            # diagnostic estimates, not rigorous enclosures of hidden stages.
            sigma = max(.005, float(np.linalg.norm(residual*weights))/math.sqrt(max(1,len(residual)-degree-1)))
            matrix = design/scales*weights[:, None]
            covariance = np.linalg.pinv(matrix.T @ matrix)
            self.fit_standard_errors = sigma*np.sqrt(np.maximum(np.diag(covariance),0))/scales
            return residual

    def record(self):
        return dict(name=self.name, origin=bits(float(self.origin)),
                    coefficients_cycles=[str(c) for c in self.coefficients],
                    active_degree=self.active_degree,
                    reference_phase=self.reference_phase,reference_amplitude=self.reference_amplitude,
                    observations=[dict(x=bits(x), cycles=str(y), amplitude=a) for x, y, a in self.observations],
                    misses=self.misses,
                    limits='Phase modulo an arbitrary constant and integer/sign slope aliases. Continuation assumes smooth Taylor coefficients of ordinary scale: each added degree spans .2 cycles and its correction exceeds six estimated standard errors (weighted noise floor .005 cycles). These are diagnostic estimates, not rigorous stage enclosures. No stage scale, constant, coefficient table or literal store is identified.')


def accept(reply, function, rows):
    # Readbacks contain actual cell values. The request's literal/cell mode is
    # transport metadata, not a field on a read-back worksheet value.
    expected=[[{k:v for k,v in a.items() if k!='mode'} for a in row] for row in rows]
    if reply.get('version')!=1 or not re.fullmatch(r'[0-9a-f]{32}',str(reply.get('id',''))):
        raise ValueError('invalid cached Excel reply identity/version')
    if not reply.get('ok') or reply.get('inputs') != expected or len(reply.get('rows', [])) != len(rows):
        raise ValueError('failed, incomplete, reordered or changed Excel operands: ' + function)
    if reply.get('automation_security_before') != reply.get('automation_security_after'):
        raise ValueError('worker changed AutomationSecurity')
    if not all(k in reply for k in ('automation_security_before','automation_security_after',
                                   'excel_build', 'machine', 'timestamp', 'request_sha256',
                                   'session_receipt_sha256')):
        raise ValueError('missing worker provenance')


def publication(bounds):
    """Certify both ideal publication candidates across a directed enclosure."""
    if bounds is None:
        return None
    lo, hi = bounds
    if lo>hi:
        raise ValueError('reference bounds reversed')
    a,b = BINARY64.round(lo),BINARY64.round(hi)
    c,d = EXTENDED64.round(lo),EXTENDED64.round(hi)
    if a is None or b is None or c is None or d is None:
        return None
    c,d=BINARY64.round(c),BINARY64.round(d)
    if a!=b or c!=d:
        raise ValueError('reference precision cannot certify publication')
    return bits(float(a)),bits(float(c))


def rounding_evidence(function, rows, answers):
    counts=Counter(); examples=[]
    for row,answer in zip(rows,answers):
        if answer['type']!='number':
            counts['non_numeric']+=1;continue
        candidates=publication(elementary_bounds(function,row))
        if candidates is None:
            counts['no_reference']+=1;continue
        single,double=candidates;got=answer['bits']
        if single==double:
            key='indistinguishable_match' if got==single else 'neither'
        else:
            key='single_only' if got==single else 'double_only' if got==double else 'neither'
            if len(examples)<6:
                examples.append(dict(arguments=row,answer=answer,single=single,double=double))
        counts[key]+=1
    if counts['no_reference']+counts['non_numeric']:
        label='unresolved: no complete scalar mathematical reference'
    elif counts['neither']:
        label='ideal single/double leaf rejected on this window'
    elif counts['single_only'] and counts['double_only']:
        label='mixed ideal single/double evidence; route unresolved'
    elif counts['single_only']:
        label='single-round mathematical leaf consistent'
    elif counts['double_only']:
        label='double-round mathematical leaf consistent'
    else:
        label='single/double indistinguishable'
    return dict(label=label,counts=counts,discriminants=examples,
                limitation='Ideal real-function publication only. Agreement does not identify the inner algorithm or its stores; any neither result refutes this complete leaf model.')


class CachedRefusal(RuntimeError):
    def __init__(self,response):
        self.response=response
        super().__init__(response.get('error','recorded service refusal'))


class CachedService:
    def __init__(self, client, directory, offline=False):
        self.client = client
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.receipts = []
        self.refusals = []

    def ask(self, function, rows, *, result_index=None):
        request = dict(function=function, rows=rows)
        if result_index is not None:
            request['result_index'] = list(result_index)
        identity = hashlib.sha256(json.dumps(request, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        path = self.directory / (identity + '.json.gz')
        pending_path=self.directory/(identity+'.pending.json')
        # Client timeouts retain queue files. A later invocation first imports
        # its own delayed successful reply, without resubmitting that question.
        if not path.exists() and pending_path.exists():
            pending=json.loads(pending_path.read_text())
            if pending['request']!=request:raise ValueError('pending request identity changed')
            for attempt in pending['attempts']:
                answer=self.client.queue/'out'/(attempt['id']+'.json') if self.client is not None else None
                if answer is None or not answer.exists():continue
                reply=json.loads(answer.read_text())
                if reply.get('id')!=attempt['id'] or reply.get('version')!=1:
                    raise ValueError('delayed reply identity changed')
                if not reply.get('ok'):continue
                if attempt.get('queue_request_sha256')!=reply.get('request_sha256'):
                    raise ValueError('delayed reply request hash changed')
                accept(reply,function,rows)
                with gzip.open(path,'xt',compresslevel=9)as stream:
                    json.dump(dict(request=request,reply=reply),stream,separators=(',',':'))
                break
        if path.exists():
            with gzip.open(path, 'rt') as stream:
                doc = json.load(stream)
            if doc['request'] != request:
                raise ValueError('cached request changed')
            reply = doc['reply']
        else:
            if self.offline:
                refused=sorted(self.directory.glob(identity+'.refusal-*.json.gz'))
                if refused:
                    selected=refused[-1]
                    with gzip.open(selected,'rt')as stream:doc=json.load(stream)
                    if doc['request']!=request or doc['reply'].get('ok')is not False:
                        raise ValueError('recorded refusal identity changed')
                    self.refusals.append(dict(function=function,file=selected.name,sha256=sha(selected),calls_requested=len(rows)))
                    raise CachedRefusal(doc['reply'])
                raise FileNotFoundError('unanswered service request ' + identity)
            try:
                reply = self.client.ask(function, rows, **(
                    {'result_index': result_index} if result_index is not None else {}))
            except TimeoutError as error:
                match=re.search(r'Excel question ([0-9a-f]{32}) timed out',str(error))
                if match is None:raise
                request_id=match.group(1)
                queued=self.client.queue/'in'/(request_id+'.json')
                attempt=dict(id=request_id,timed_out_utc=time.time(),
                             queue_request_sha256=sha(queued)if queued.exists()else None)
                pending=json.loads(pending_path.read_text())if pending_path.exists()else dict(request=request,attempts=[])
                pending['attempts'].append(attempt);write(pending_path,pending)
                raise
            except RuntimeError as error:
                refusal=getattr(error,'response',None)
                if isinstance(refusal,dict) and refusal.get('ok') is False:
                    refused=self.directory/(identity+'.refusal-'+str(refusal.get('id','unknown'))+'.json.gz')
                    with gzip.open(refused,'xt',compresslevel=9)as stream:
                        json.dump(dict(request=request,reply=refusal),stream,separators=(',',':'))
                    self.refusals.append(dict(function=function,file=refused.name,sha256=sha(refused),calls_requested=len(rows)))
                raise
            with gzip.open(path, 'xt', compresslevel=9) as stream:
                json.dump(dict(request=request, reply=reply), stream, separators=(',', ':'))
        accept(reply, function, rows)
        self.receipts.append(dict(function=function, file=path.name, sha256=sha(path), calls=len(rows),
                                  excel_build=reply['excel_build'], machine=reply['machine'],
                                  timestamp=reply['timestamp'], request_sha256=reply['request_sha256']))
        return reply

    def ask_rows(self, function, rows, chunk_size=10000, *, result_index=None):
        """Keep each service request and receipt intact within its call limit."""
        if not 1<=chunk_size<=20000:raise ValueError('invalid scalar request chunk size')
        return [answer for start in range(0,len(rows),chunk_size)
                for answer in self.ask(function,rows[start:start+chunk_size],
                                       result_index=result_index)['rows']]


def scan(checkout, output, cache, count=512, offline=False, only=None):
    client, policy = load_service(checkout)
    checkout, output = Path(checkout).resolve(), Path(output)
    problem = checkout / 'docs/generated/problem-functions.json'
    source_paths=[Path(__file__),ROOT/'checker/reference.py',ROOT/'checker/exact.py',problem,
                  checkout/'tools/excel-service/client.py',checkout/'tools/excel-service/policy.py']
    source_bindings={str(p):sha(p)for p in source_paths}
    inventory = json.loads(problem.read_text())
    names = [r['function'] for r in inventory['functions']]
    # Explicit primitive controls supplement, but do not replace, list members.
    names += [f for f in ('EXP', 'LN', 'LOG', 'LOG10', 'POWER', 'SQRT') if f not in names]
    if only:
        wanted = set(only.split(','))
        names = [f for f in names if f in wanted]
    allowed = policy.catalogue_names(checkout)
    service = CachedService(client.Client(timeout=900), cache, offline)
    results = []
    for function in names:
        record = dict(function=function, problem_list_member=function in {r['function'] for r in inventory['functions']},
                      rounding='unresolved', stored_stages=None, band_edges=[], windows=[])
        if function not in allowed:
            record.update(status='workbook-only', reason='outside the scalar service catalogue; no service restriction changed')
        else:
            args, index, origin = seed(function, checkout, client)
            record['seed_origin'] = origin
            if args is None or index is None:
                record.update(status='non-numeric-coordinate', reason=origin)
            else:
                # Two separated windows per scalar coordinate; retain the
                # exact seed arguments even where integer truncation is flat.
                center = value(args[index]['bits'])
                window_count=max(count,4096) if function in ('EXP','LN','LOG','LOG10','POWER','SQRT') else count
                starts = [int(bits(center), 16)-window_count//2,
                          int(bits(center*1.25), 16)-window_count//2]
                rows = []
                for first in starts:
                    for b in range(first, first+window_count):
                        row = [dict(a) for a in args]
                        row[index] = client.number(f'{b:016x}')
                        rows.append(row)
                record.update(arguments=args,varying_argument=index,requested_calls=len(rows),
                              requested_windows=[dict(start=f'{first:016x}',count=window_count)for first in starts])
                try:
                    reply = service.ask(function, rows)
                except RuntimeError as e:
                    if not isinstance(getattr(e,'response',None),dict) or e.response.get('ok')is not False:
                        raise
                    # An actual refusal is retained as such, never credited as
                    # a measured window. Do not silently replace its operands.
                    record.update(status='service-refused', refusal=e.response)
                else:
                    record.update(status='observed', arguments=args, varying_argument=index)
                    for i, first in enumerate(starts):
                        answers = reply['rows'][i*window_count:(i+1)*window_count]
                        types = Counter(a['type'] for a in answers)
                        record['windows'].append(dict(start=f'{first:016x}', end=f'{first+window_count-1:016x}',
                            count=window_count, answer_types=types, distinct_answers=len({json.dumps(a, sort_keys=True) for a in answers}),
                            spectrum=spectrum(answers)))
                    record['rounding_evidence']=rounding_evidence(function,rows,reply['rows'])
                    record['rounding']=record['rounding_evidence']['label']
                    if any(w['answer_types'].get('number', 0) != window_count for w in record['windows']):
                        record['rounding'] = 'not identifiable on non-numeric/error output'
                    elif all(w['distinct_answers'] == 1 for w in record['windows']):
                        record['rounding'] = 'unresolved: flat quantized windows'
                print(function, record['status'], flush=True)
        results.append(record)
        # A restart resumes the immutable cache, while the partial report
        # explicitly remains incomplete if any request has not been visited.
        write(output.with_suffix('.partial.json'), dict(complete=False, functions=results, receipts=service.receipts))
    if source_bindings!={str(p):sha(p)for p in source_paths}:
        raise ValueError('fingerprint producer changed during run; cache remains valid, regenerate report')
    edge_path=Path(cache).resolve().parent/'edges.json'
    bands_path=Path(cache).resolve().parent/'bessel-bands.json'
    log_path=Path(cache).resolve().parent/'log-consumers.json'
    auxiliary={}
    if edge_path.exists():
        edge_doc=json.loads(edge_path.read_text());auxiliary[str(edge_path)]=sha(edge_path)
        for record in results:
            record['band_edges']=[e for e in edge_doc.get('edges',[])if e['function']==record['function']and e['status']=='adjacent observed edge']
    if bands_path.exists():
        band_doc=json.loads(bands_path.read_text());auxiliary[str(bands_path)]=sha(bands_path)
        for record in results:
            record['observed_approximation_boundaries']=[b for b in band_doc.get('bands',[])if b['function']==record['function']]
    if log_path.exists():
        log_doc=json.loads(log_path.read_text())
        if log_doc.get('complete'):
            auxiliary[str(log_path)]=sha(log_path)
            for record in results:
                if record['function']!='LOG':continue
                record['conditional_consumer_study']=log_doc
                division=Counter()
                for row in log_doc['records']:division.update(row['division_rounding'])
                record['consumer_division_rounding']=division
                coherent=all(any(p['shape']=='divide'and p['model']=='stored_pc64'and p['consistent']
                                 for p in row['profiles'])for row in log_doc['records'])
                if coherent and division['double_only']and not(division['single_only']or division['neither']):
                    record['rounding']='double division in the public-LN quotient model; ideal leaf rejected'
    primitive_path=ROOT/'evidence/round27/primitives.json'
    if primitive_path.is_file():
        primitive_doc=json.loads(primitive_path.read_text())
        if primitive_doc.get('complete'):
            auxiliary[str(primitive_path)]=sha(primitive_path)
            for record in results:
                record['broader_primitive_windows']=[r for r in primitive_doc['leaves']
                    if r['function']==record['function']]
    stage_path=ROOT/'evidence/round28/primitives-final-independent.json'
    if stage_path.is_file():
        stage_paths=[stage_path,ROOT/'evidence/round28/root-normalization-independent.json',
                     ROOT/'evidence/round28/hypergeom-final-independent.json',
                     ROOT/'evidence/round28/poisson-final-independent.json',
                     ROOT/'evidence/round28/root-public-stages-independent.json',
                     ROOT/'evidence/round28/root-store-contrasts-public-exp.json',
                     ROOT/'evidence/round28/poisson-exp-arguments.json',
                     ROOT/'evidence/round28/poisson-log-difference.json']
        stages=[]
        for i,p in enumerate(stage_paths):
            archived=p.with_suffix(p.suffix+'.gz')
            if archived.is_file():p=archived;stage_paths[i]=p
            with (gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as stream:
                stages.append(json.load(stream))
        if not all(d.get('complete')is True for d in stages):
            raise ValueError('incomplete independently checked R28 stages')
        auxiliary.update({str(p):sha(p)for p in stage_paths})
        primitive,normalization,hypergeom,poisson,paired,contrasts,poisson_exp,poisson_log=stages
        for record in results:
            name=record['function']
            if name in ('EXP','LN'):
                constants=(primitive['refined_log_constants']if name=='LN'else
                    next(f['first_pc64_candidates']for f in primitive['constant_families']if f['function']==name))
                record['isolated_stage']=dict(constant_pc64_little_endian=constants,
                    model=('FYL2X with resident ln(2)'if name=='LN'else
                           'resident multiply by log2(e), F2XM1 reduction and FSCALE'),
                    scope='Unique coefficient within the declared native instruction model on the measured windows; no private helper identification.')
            elif name=='GEOMEAN':
                record['isolated_stage']=dict(normalization='e=ilogb(P); q=trunc(e/n); M=P*2^(-q*n); root(M,n)*2^q',
                    original_calls=normalization['original_calls'],normalized_controls=normalization['normalized_calls'],
                    convention_counts=normalization['counts'],
                    bounded_core='Store LN(M), store 1/n from PC64, store their product from PC64, call the measured public EXP; singleton identity.',
                    paired_counts=paired['counts'],singleton_counts=paired['singleton_counts'],
                    store_contrasts=contrasts['counts'],
                    scope='Measured positive normal stored products and unit factors. The sampled normalized functional composition is identified; literal private implementation internals and the product underflow/overflow fallback remain open. Earlier shared-log contradictions concern the unreduced direct factorization only.')
            elif name in ('HYPGEOM.DIST','HYPGEOMDIST'):
                record['isolated_stage']=dict(counts=hypergeom['counts'],scope=hypergeom['limits'])
            elif name in ('POISSON.DIST','POISSON'):
                record['isolated_stage']=dict(counts=poisson['counts'],exp_argument_counts=poisson_exp['counts'],
                    log_difference_counts=poisson_log['counts'],
                    scope='Conditional public EXP preimages and one subtraction isolate necessary stage constraints. Confirmed adjacent public EXP gaps exclude a stored binary64 exponential input within that primitive. Neither the private PMF nor its log/exp implementation is identified.')
    result = dict(format='dnacalc-numeric-fingerprints/1', complete=True,inventory_complete=not bool(only),
        coverage=dict(problem_functions=len(inventory['functions']), reported=len(results),
                      statuses=Counter(r['status'] for r in results)),
        generator_sha256=sha(__file__), problem_list_sha256=sha(problem),
        source_bindings=source_bindings,
        auxiliary_bindings=auxiliary,
        refused_requests=service.refusals,
        sources='Original Excel scalar service outputs; IEEE754 lattice; original MIT rational arithmetic/DFT. No third-party implementation read.',
        candidate_order='Textbook product loops, Horner, plain recurrences, raw-moment variance, atanh series, one guard; fixed structural order, never accuracy rank.',
        limitations=['Exact single/double-rounding identification requires discriminating computation-tree predictions; output spectra alone cannot supply it.',
                    'Detected spectral lines are not identified stores. Absence of a line cannot exclude a stored stage.',
                    'Two local windows do not establish global route boundaries. Each unmeasured edge remains empty.',
                    'The EXP double-rounding label applies to its listed local discriminants. R27 broader windows contain74 ideal-leaf contradictions among2048 inputs; the registered x87 chain matches every measured input. An ideal mathematical EXP leaf is not a global identification.',
                    'Separate-cell scalar observations do not establish range, array, reference or workbook behavior.',
                    'All problem-list entries are reported, including workbook-only, text/discrete and unidentifiable cases.'],
        cache_files=service.receipts, functions=results)
    write(output.with_suffix('.json'), result)
    render(result, output.with_suffix('.md'))
    print(json.dumps(result['coverage']), flush=True)
    return result


def edges(checkout, output, cache, offline=False):
    """Locate published numeric/error transitions; no hidden route is assumed.

    The supplied brackets come from previous original Bessel domain readings
    and the textbook EXP overflow neighborhood. All endpoints are re-read.
    The bisection proves local adjacency, not global monotonicity or source.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    plan=[dict(function='EXP',low=709.,high=710.,arguments=[],index=0)]
    for function in ('COS','SIN'):
        plan.append(dict(function=function,low=float(2**27-2),high=float(2**27+4),arguments=[],index=0))
    for function in ('BESSELJ','BESSELY'):
        for order in (0,1,2,3,10):
            plan.append(dict(function=function,low=float(2**27-2),high=float(2**27+4),
                             arguments=[client.number(float(order))],index=0))
    results=[]
    predicate=lambda a:False if a['type']=='number'else True if a==dict(type='error',code='#NUM!')else None
    def row(p,b):
        args=list(p['arguments']);args.insert(p['index'],client.number(f'{b:016x}'));return args
    # Several orders share a function and one Excel batch. All independent
    # brackets advance together, rather than paying one queue round per order.
    for function in dict.fromkeys(p['function']for p in plan):
        states=[dict(plan=p,lo=int(bits(p['low']),16),hi=int(bits(p['high']),16),steps=0)
                for p in plan if p['function']==function]
        endpoints=service.ask(function,[row(s['plan'],b)for s in states for b in (s['lo'],s['hi'])])['rows']
        for s,a,b in zip(states,endpoints[::2],endpoints[1::2]):
            if predicate(a)is not False or predicate(b)is not True:
                s['unresolved']='endpoints do not bracket a numeric-to-NUM transition'
            s['low_answer'],s['high_answer']=a,b
        while True:
            active=[s for s in states if 'unresolved'not in s and s['hi']-s['lo']>1]
            if not active:break
            middles=[(s['lo']+s['hi'])//2 for s in active]
            replies=service.ask(function,[row(s['plan'],b)for s,b in zip(active,middles)])['rows']
            for s,b,answer in zip(active,middles,replies):
                outcome=predicate(answer);s['steps']+=1
                if outcome is None:s['unresolved']='unclassifiable midpoint; no side guessed'
                elif outcome:s['hi'],s['high_answer']=b,answer
                else:s['lo'],s['low_answer']=b,answer
            write(output,dict(complete=False,edges=results,states=states,receipts=service.receipts))
        completed=[s for s in states if 'unresolved'not in s]
        neighbors=service.ask(function,[row(s['plan'],b)for s in completed
                                        for b in range(s['lo']-2,s['hi']+3)])['rows']if completed else[]
        for i,s in enumerate(completed):
            answers=neighbors[i*6:(i+1)*6]
            if list(map(predicate,answers))!=[False]*3+[True]*3:
                s['unresolved']='neighbors contradict a local numeric-to-NUM transition'
            s['neighbors']=answers
        for s in states:
            p=s['plan']
            if 'unresolved'in s:record=dict(**p,status='unresolved bracket',reason=s['unresolved'])
            else:record=dict(**p,status='adjacent observed edge',
                edge=dict(low=f"{s['lo']:016x}",high=f"{s['hi']:016x}",low_answer=s['low_answer'],
                          high_answer=s['high_answer'],questions=s['steps']+2),neighbors=s['neighbors'])
            results.append(record)
            print(function,p['arguments'],record['status'],flush=True)
        write(output,dict(complete=False,edges=results,receipts=service.receipts))
    result=dict(format='dnacalc-fingerprint-edges/1',complete=True,edges=results,receipts=service.receipts,
                limitations='Exact local numeric/#NUM transitions only, under one-transition bracket assumption with two exterior neighbors rechecked. No global monotonicity, hidden algorithm, or threshold rounding expression inferred.')
    write(output,result)
    return result


def bessel_bands(checkout,output,cache,offline=False,count=256):
    """Adjacent-input windows at the published seed-family boundary x=8.

    Exact residual jumps diagnose an approximation change. A jump is an
    observation; attributing it to an internal branch remains an inference.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    center=int(bits(8.),16);bands=[]
    for function in ('BESSELJ','BESSELY'):
        rows=[[client.number(f'{b:016x}'),client.number(float(order))]for order in (0,1)
              for b in range(center-count,center+count)]
        answers=service.ask(function,rows)['rows']
        for j,order in enumerate((0,1)):
            q=rows[j*2*count:(j+1)*2*count];a=answers[j*2*count:(j+1)*2*count]
            if any(x['type']!='number'for x in a):raise ValueError('Bessel boundary window contains a non-numeric answer')
            residuals=[]
            for args,answer in zip(q,a):
                lo,hi=elementary_bounds(function,args)
                residuals.append((Q(value(answer['bits'])))-(lo+hi)/2)
            unit=Q(math.ulp(value(a[count]['bits'])))
            jump=(residuals[count]-residuals[count-1])/unit
            bands.append(dict(function=function,order=order,
                low=f'{center-1:016x}',high=f'{center:016x}',
                low_answer=a[count-1],high_answer=a[count],adjacent_residual_jump_ulp=float(jump),
                left_spectrum=spectrum(a[:count]),right_spectrum=spectrum(a[count:]),
                inference='A large residual discontinuity supports a change of approximation at8; output readings alone do not identify the internal source branch.'))
            print(function,order,'adjacent residual jump',float(jump),'ULP',flush=True)
    result=dict(format='dnacalc-bessel-boundary-windows/1',complete=True,bands=bands,receipts=service.receipts,
        limitations='Published candidate boundary at8, orders0/1, positive normal scalar cells. Mathematical MPFR references are diagnostics. No table coefficient or hidden branch selector is fit.')
    write(output,result)
    return result


def root_study(checkout, output, cache, offline=False, count=256):
    """Unit-factor GEOMEAN probes remove all product-rounding uncertainty.

    Candidate order is fixed by code shape: product plus root, public POWER,
    public LN/divide/EXP, then public LN/reciprocal/multiply/EXP. Live public
    helpers are conditional leaves and are not substituted into product code.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    counts=tuple(range(2,33))
    centers=(0.13,1.3,7.0,1000.)
    inputs=[f'{b:016x}' for center in centers
            for b in range(int(bits(center),16)-count//2,int(bits(center),16)+count//2)]
    logarithms=service.ask('LN',[[client.number(x)]for x in inputs])['rows']
    roots=[]
    for n in counts:
        questions=[[client.number(x)]+[client.number(1.)]*(n-1)for x in inputs]
        answers=service.ask('GEOMEAN',questions)['rows']
        pows=service.ask('POWER',[[client.number(x),client.number(1./n)]for x in inputs])['rows']
        candidates=[dict(name='correct real root, one final rounding',kind='mathematical'),
                    dict(name='correct real root, extended64 then binary64',kind='mathematical'),
                    dict(name='public POWER(p, binary64(1/n))',kind='conditional helper')]
        if n==2:candidates.append(dict(name='public SQRT(p)',kind='conditional helper'))
        for association in ('divide','reciprocal_multiply'):
            for model in ('binary64','stored_pc64'):
                candidates.append(dict(name=f'public EXP({model}(public LN(p) {association} n))',
                                       kind='conditional helper',association=association,model=model))
        predictions=[]
        real=[publication(elementary_bounds('GEOMEAN',row))for row in questions]
        predictions.extend([[p[i]if p is not None else None for p in real]for i in (0,1)])
        predictions.append([a.get('bits')for a in pows])
        if n==2:predictions.append([a.get('bits')for a in service.ask('SQRT',[[client.number(x)]for x in inputs])['rows']])
        for candidate in candidates[len(predictions):]:
            def rounded(q):
                if candidate['model']=='stored_pc64':q=EXTENDED64.round(q)
                return BINARY64.round(q)
            reciprocal=rounded(Q(1,n))
            qrows=[]
            for a in logarithms:
                if a['type']!='number':raise ValueError('non-numeric logarithm in root study')
                q=Q(value(a['bits']))
                q=rounded(q/n if candidate['association']=='divide'else q*reciprocal)
                qrows.append([client.number(bits(float(q)))])
            predictions.append([a.get('bits')for a in service.ask('EXP',qrows)['rows']])
        for candidate,prediction in zip(candidates,predictions):
            tally=Counter();first=None
            for x,a,p in zip(inputs,answers,prediction):
                if a['type']!='number' or p is None:key='unresolved'
                else:key='matched'if a['bits']==p else'mismatched'
                tally[key]+=1
                if key=='mismatched'and first is None:first=dict(input=x,excel=a['bits'],candidate=p)
            candidate.update(counts=tally,first_contradiction=first,consistent=tally['matched']==len(inputs))
        roots.append(dict(count=n,observations=len(inputs),candidates=candidates))
        write(output,dict(complete=False,roots=roots,receipts=service.receipts))
        print('GEOMEAN unit-factor n',n,[(c['name'],dict(c['counts']))for c in candidates],flush=True)
    result=dict(format='dnacalc-geomean-unit-factor-study/1',complete=True,roots=roots,receipts=service.receipts,
        input_centers=list(centers),consecutive_inputs_per_center=count,
        limits='All but one argument equal one, so every product stage is exact. Public helper compositions are specified conditional programs; their failure does not refute an unpublished extended helper. Four fixed windows are discovery, not universal equivalence. No fitted constant, store mask, tolerance, row selector or accuracy ranking.')
    write(output,result)
    return result


def product_study(checkout,output,cache,offline=False,count=128):
    """Exact power-of-two redistributions test a product-before-root tree.

    Every input program below forms precisely the same p in strict binary64,
    stored PC64 and resident PC64, in any product order. A different answer
    therefore refutes a pure product plus count-dependent root on this scope,
    without assuming which root helper Excel uses.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    centers=(.13,1.3,7.,1000.)
    xs=[value(f'{b:016x}')for x in centers for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    patterns=[('p,1',lambda p:[p,1.]),('2p,1/2',lambda p:[2*p,.5]),
              ('1/2,2p',lambda p:[.5,2*p]),('1,p',lambda p:[1.,p]),
              ('p,1,1',lambda p:[p,1.,1.]),('2p,1/2,1',lambda p:[2*p,.5,1.]),
              ('p,2,1/2',lambda p:[p,2.,.5]),('4p,1/4,1',lambda p:[4*p,.25,1.]),
              ('p,4,1/4',lambda p:[p,4.,.25]),('1,p,1',lambda p:[1.,p,1.]),
              ('1,1,p',lambda p:[1.,1.,p]),('p,1,1,1',lambda p:[p,1.,1.,1.]),
              ('p,2,1/2,1',lambda p:[p,2.,.5,1.])]
    for n in (9,10,17,32):
        patterns.extend([(f'p,ones n={n}',lambda p,n=n:[p]+[1.]*(n-1)),
                         (f'2p,half,ones n={n}',lambda p,n=n:[2*p,.5]+[1.]*(n-2)),
                         (f'half,ones,2p n={n}',lambda p,n=n:[.5]+[1.]*(n-2)+[2*p])])
    questions=[]
    for name,pattern in patterns:
        for p in xs:
            arguments=pattern(p)
            for reverse in (False,True):
                product=Q(1)
                for x in arguments[::-1]if reverse else arguments:
                    exact=product*Q(x)
                    if BINARY64.round(exact)!=exact or EXTENDED64.round(exact)!=exact:
                        raise ValueError('product-invariance control has an inexact multiplication')
                    product=exact
                if product!=Q(p):raise ValueError('product-invariance control changes its mathematical input')
            questions.append([client.number(x)for x in arguments])
    # Retain the original thirteen-pattern request identity, then extend it.
    answers=[]
    for first in range(0,len(questions),13*len(xs)):
        answers.extend(service.ask('GEOMEAN',questions[first:first+13*len(xs)])['rows'])
    baseline={};records=[]
    for j,(name,pattern)in enumerate(patterns):
        n=len(pattern(xs[0]));group=answers[j*len(xs):(j+1)*len(xs)]
        if n not in baseline:baseline[n]=group
        tally=Counter();first=None
        for p,a,b in zip(xs,group,baseline[n]):
            key='equal'if a==b else'different';tally[key]+=1
            if key=='different'and first is None:first=dict(input=bits(p),baseline=b,redistributed=a)
        records.append(dict(pattern=name,count=n,counts=tally,first_difference=first))
        print('GEOMEAN product invariant',name,dict(tally),flush=True)
    result=dict(format='dnacalc-geomean-product-invariance/1',complete=True,records=records,receipts=service.receipts,
        consecutive_inputs_per_center=count,centers=centers,
        limitations='Separate scalar Value2 cells, normal finite numbers only. A difference refutes the specified pure product-plus-root composition; it does not identify a logarithm, branch, or alternate root. Identical publications do not prove the composition generally. No expected answer, threshold or coefficient fitted.')
    write(output,result)
    return result


def product_store_study(checkout,output,cache,offline=False):
    """Measure GEOMEAN's product independently of its unknown root helper.

    For each declared multiplication tie, ask the original factors AND two
    unit-factor calls at the complete strict/storedPC64 products. These root
    controls are new independent Excel questions, never substituted answers.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    seed_bits='401f9ff8b1ae9622'
    seed=Q(value(seed_bits))
    programs=[]
    # The known square tie was independently found in R25 variance research.
    # Fixed exact powers of two move its root without changing that tie.
    for k in range(-30,31):
        scale=Q(2)**k
        xs=[seed,seed,scale]
        strict=stored=Q(1)
        for x in xs:
            strict=BINARY64.round(strict*x)
            stored=BINARY64.round(EXTENDED64.round(stored*x))
        if strict==stored:raise ValueError('declared product tie disappeared')
        programs.append(dict(args=[bits(float(x))for x in xs],strict=bits(float(strict)),stored_pc64=bits(float(stored)),scale_exponent=k))
    questions=[]
    for p in programs:
        questions.append([client.number(x)for x in p['args']])
        for name in ('strict','stored_pc64'):
            questions.append([client.number(p[name]),client.number(1.),client.number(1.)])
    answers=service.ask('GEOMEAN',questions)['rows']
    tally=Counter();records=[]
    for i,p in enumerate(programs):
        target,strict,stored=answers[3*i:3*i+3]
        if strict==stored:outcome='root publications indistinguishable'
        elif target==strict:outcome='strict product only'
        elif target==stored:outcome='storedPC64 product only'
        else:outcome='neither product root control'
        tally[outcome]+=1
        records.append(dict(**p,outcome=outcome,original=target,strict_root_control=strict,stored_root_control=stored))
    result=dict(format='dnacalc-geomean-product-store-isolation/1',complete=True,counts=tally,records=records,receipts=service.receipts,
        limitations='One independently known multiplication tie under61 predeclared exact scale shifts. Root controls are measured through the same GEOMEAN argument count; no root implementation is assumed. Agreement is local evidence about this product tree, not a complete GEOMEAN kernel. No program selected by match counts.')
    write(output,result);print(json.dumps(dict(tally)),flush=True)
    return result


def product_stage_study(checkout,output,cache,offline=False,count=256):
    """Non-power third factors separate per-step stores from resident PC64.

    Three complete product programs are fixed before asking Excel. Every
    candidate's binary64 product gets its own live unit-factor root control.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    x=Q(value('401f9ff8b1ae9622'));programs=[];questions=[]
    models=('strict','stored_each','resident_final')
    for center in (1.3,3.7):
        first=int(bits(center),16)-count//2
        for b in range(first,first+count):
            args=[x,x,Q(value(f'{b:016x}'))]
            products={}
            for model in models:
                product=Q(1)
                for a in args:
                    product=EXTENDED64.round(product*a)if model=='resident_final'else BINARY64.round(
                        EXTENDED64.round(product*a)if model=='stored_each'else product*a)
                products[model]=bits(float(BINARY64.round(product)))
            encoded=[bits(float(a))for a in args]
            programs.append(dict(args=encoded,products=products,center=center))
            questions.append([client.number(a)for a in encoded])
            for model in models:questions.append([client.number(products[model]),client.number(1.),client.number(1.)])
    replies=service.ask('GEOMEAN',questions)['rows']
    counts=Counter();records=[]
    for i,program in enumerate(programs):
        target,*roots=replies[4*i:4*i+4]
        compatible=[m for m,r in zip(models,roots)if target==r]
        label='neither'if not compatible else'only '+compatible[0]if len(compatible)==1 else'ambiguous '+','.join(compatible)
        counts[label]+=1
        records.append(dict(**program,original=target,root_controls=dict(zip(models,roots)),compatible=compatible))
    result=dict(format='dnacalc-geomean-intermediate-product-stores/1',complete=True,models=list(models),counts=counts,
        centers=[1.3,3.7],consecutive_inputs_per_center=count,records=records,receipts=service.receipts,
        limitations='Known square tie followed by a non-power factor, positive normal scalar inputs, three factors. Distinguishes declared product programs with a binary64 product at the root boundary. A resident product consumed by an unpublished extended root is not identified by these controls. No source threshold, coefficient or per-row program selected.')
    write(output,result);print(json.dumps(dict(counts)),flush=True)
    return result


def poisson_study(checkout,output,cache,offline=False,count=512,centers=None):
    """Fixed zero/one-count textbook products; independent EXP readings.

    Measured public EXP is an explicitly conditional binary64 leaf. Surviving
    products identify only this stated consumer on the sampled windows.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    centers=(.13,1.3,3.7,10.,50.)if centers is None else tuple(centers)
    if count<1 or not centers or any(not math.isfinite(x)or x<=0 for x in centers):raise ValueError('nonempty positive Poisson mapping required')
    xs=[f'{b:016x}'for x in centers for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    exponential=service.ask('EXP',[[client.number(bits(-value(x)))]for x in xs])['rows']
    jobs=[(0,False),(0,True),(1,False),(1,True),(2,False),(2,True)]
    questions=[[client.number(float(k)),client.number(x),client.boolean(cdf)]for k,cdf in jobs for x in xs]
    answers=service.ask('POISSON.DIST',questions)['rows']
    records=[]
    def rounding(q,model):return BINARY64.round(q if model=='binary64'else EXTENDED64.round(q))
    for j,(k,cdf)in enumerate(jobs):
        profiles=[]
        for shape in ('finite_sum_then_product','recurrence_from_exp'):
            for model in ('binary64','stored_pc64'):
                tally=Counter();first=None;bands=[Counter()for _ in centers]
                for input_index,(x,e,a) in enumerate(zip(xs,exponential,answers[j*len(xs):(j+1)*len(xs)])):
                    if e['type']!='number'or a['type']!='number':tally['unresolved']+=1;continue
                    mean,leaf=Q(value(x)),Q(value(e['bits']))
                    op=lambda q:rounding(q,model)
                    if shape=='finite_sum_then_product':
                        term=Q(1);s=term if cdf else None
                        for i in range(1,k+1):
                            term=op(op(term*mean)/i)
                            if cdf:s=op(s+term)
                        result=op(leaf*(s if cdf else term))
                    else:
                        term=leaf;s=term if cdf else None
                        for i in range(1,k+1):
                            term=op(op(term*mean)/i)
                            if cdf:s=op(s+term)
                        result=s if cdf else term
                    predicted=bits(float(result));key='matched'if predicted==a['bits']else'mismatched';tally[key]+=1
                    bands[input_index//count][key]+=1
                    if key=='mismatched'and first is None:first=dict(input=x,expo=e['bits'],excel=a['bits'],candidate=predicted)
                profiles.append(dict(shape=shape,model=model,counts=tally,first_contradiction=first,
                                     consistent=tally['matched']==len(xs),
                                     bands=[dict(center=bits(x),counts=dict(c))for x,c in zip(centers,bands)]))
        records.append(dict(count=k,cumulative=cdf,profiles=profiles))
        print('POISSON',k,cdf,[(p['shape'],p['model'],dict(p['counts']))for p in profiles],flush=True)
    logarithms=service.ask('LN',[[client.number(x)]for x in xs]+[[client.number(2.)]])['rows']
    log_two=Q(value(logarithms[-1]['bits']))
    log_programs=[]
    for k in (1,2):
        for model in ('binary64','stored_pc64'):
            op=lambda q:rounding(q,model)
            exponent_rows=[]
            for x,log in zip(xs,logarithms[:-1]):
                t=op(op(k*Q(value(log['bits'])))-Q(value(x)))
                if k==2:t=op(t-log_two)
                exponent_rows.append([client.number(bits(float(t)))])
            predicted=service.ask('EXP',exponent_rows)['rows']
            target=answers[jobs.index((k,False))*len(xs):(jobs.index((k,False))+1)*len(xs)]
            tally=Counter();first=None
            for x,p,a in zip(xs,predicted,target):
                key='matched'if p==a else'mismatched';tally[key]+=1
                if key=='mismatched'and first is None:first=dict(input=x,excel=a,candidate=p)
            log_programs.append(dict(count=k,model=model,counts=tally,first_contradiction=first,
                                     expression='EXP(k*public LN(mean)-mean-public LN(k!)); every arithmetic result stored'))
            print('POISSON log program',k,model,dict(tally),flush=True)
    result=dict(format='dnacalc-poisson-textbook-consumers/1',complete=True,program_order='Finite sum then product; recurrence from exponential. Binary64 then storedPC64. Fixed structural order, not accuracy rank.',
        centers=centers,consecutive_inputs_per_center=count,records=records,log_programs=log_programs,receipts=service.receipts,
        count_one_rows=[dict(input=x,public_exp=e,pmf=answers[jobs.index((1,False))*len(xs)+i],cdf=answers[jobs.index((1,True))*len(xs)+i])for i,(x,e)in enumerate(zip(xs,exponential))],
        limitations='Positive finite small means and counts0..2 only. EXP readings are public scalar binary64 leaves, not recovered private Ext80 seeds. A mismatch rejects only the specified consumer. No threshold, polynomial, recurrence selector or kernel is fit to counts.')
    write(output,result)
    return result


def integer_gamma_horner_graph(count,association=0):
    """The finite gamma sum in backward Horner order, before the exponential."""
    g=Builder(2);s=g.one
    for k in range(count,0,-1):
        divisor=g.constant(bits(float(k)))
        if association==0:t=g.op('div',g.op('mul',0,s),divisor)
        elif association==1:t=g.op('mul',g.op('div',0,divisor),s)
        else:t=g.op('mul',0,g.op('div',s,divisor))
        s=g.op('add',g.one,t)
    final=g.op('mul',1,s)
    return g.graph(final,final)


def integer_gamma_horner_study(study,cache,library,output):
    """Test the Horner shape on all existing finite-sum mapping observations."""
    import numpy as np
    document=json.loads(Path(study).read_text());service=CachedService(None,cache,True)
    exponential={};observations=[]
    for receipt in document['receipts']:
        path=Path(cache)/receipt['file']
        if sha(path)!=receipt['sha256']:raise ValueError('original service receipt changed')
        with gzip.open(path,'rt')as stream:doc=json.load(stream)
        request=doc['request']
        if request['function']not in('EXP','POISSON.DIST'):continue
        reply=service.ask(request['function'],request['rows'])
        for args,answer in zip(request['rows'],reply['rows']):
            if answer['type']!='number':raise ValueError('finite numeric mapping required')
            if request['function']=='EXP':exponential[args[0]['bits']]=answer['bits']
            elif args[2].get('type')=='boolean'and args[2].get('value')is True:
                k=int(value(args[0]['bits']));mean=value(args[1]['bits'])
                if 1<=k<=18 and k+1<=mean<=30.:observations.append((k,args[1]['bits'],answer['bits']))
    if not observations:raise ValueError('empty integer finite-sum observation scope')
    native=Native(library);profiles=[]
    for association in range(3):
        for model,cw in [('strict',0x123f),('stored_pc64',0x133f),('resident_pc64',0x133f)]:
            first=None;checked=0;contradictions=Counter();examples=[]
            for count in sorted({r[0]for r in observations}):
                rows=[r for r in observations if r[0]==count]
                inputs=np.array([[int(x,16),int(exponential[bits(-value(x))],16)]for _,x,_ in rows],dtype=np.uint64)
                graph=integer_gamma_horner_graph(count,association)
                prediction=native.grid(graph,inputs,cw,0 if model=='resident_pc64'else(1<<graph['sites'])-1)[:,0]
                for (_,x,wanted),got in zip(rows,prediction):
                    checked+=1
                    if int(got)!=int(wanted,16):
                        contradiction=dict(count=count,input=x,expected=wanted,prediction=f'{int(got):016x}')
                        if first is None:first=contradiction
                        contradictions['at_lower_edge'if value(x)==count+1 else 'interior']+=1
                        if len(examples)<8:examples.append(contradiction)
            profiles.append(dict(association=association,model=model,rows=checked,first_contradiction=first,
                contradiction_counts=contradictions,examples=examples))
    result=dict(format='dnacalc-integer-gamma-horner-discovery/1',complete=True,profiles=profiles,
        survivors=[p for p in profiles if p['first_contradiction']is None],receipts=service.receipts,
        source_bindings={str(p):sha(p)for p in(Path(__file__),Path(study),Path(library))},
        limits='Backward finite polynomial from the mapping, conditioned on separately observed EXP(-mean). The declared integer finite-sum region is count+1<=mean<=30; its edges still require adjacent and new-count controls. All prior mapping answers are discovery, not fresh validation.')
    write(output,result);print(json.dumps(dict(profiles=profiles)),flush=True)


def integer_gamma_sum_validation(checkout,cache,library,output,polynomial=False,seed=35099,points=64,max_count=19,min_count=1):
    """Freeze Horner predictions for fresh counts, branch edges and gamma aliases."""
    import random
    import numpy as np
    if not 1<=min_count<=max_count<=29 or points<1:raise ValueError('nonempty bounded integer-gamma mapping required')
    rng=random.Random(seed);native=Native(library)
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=180),cache)
    cases=[]
    for count in range(min_count,max_count+1):
        a=float(count+1)
        means=[a+(30.-a)*rng.random()for _ in range(points)]
        means += [value(f'{int(bits(edge),16)+offset:016x}')for edge in(a,30.)for offset in(-2,-1,0,1,2)]
        if polynomial and count==1:
            means += [value(f'{int(bits(3.25),16)+offset:016x}')for offset in range(-64,65)]
        for mean in means:
            scope='interior'if a<mean<=30. and a<=20 else 'edge_or_large_shape'
            cases.append(dict(count=count,shape=a,mean=bits(mean),scope=scope))
    frozen={str(p):sha(p)for p in(Path(__file__),Path(library))}
    write(str(output)+'.questions.json',dict(format='dnacalc-integer-gamma-sum-validation-questions/1',
        cases=cases,seed=seed,source_bindings=frozen,
        selection=f'Fixed counts{min_count}through{max_count},{points} uniform mapping points per count from seed{seed}; repeated edge controls and adjacent doubles at the mathematical lower edge and observed upper edge30. Coefficient separation also uses129 adjacent inputs near3.25 at count1 when polynomial=True. No answer-guided selection.'))
    exps=service.ask_rows('EXP',[[client.number(bits(-value(c['mean'])))]for c in cases])
    if any(e['type']!='number'for e in exps):raise ValueError('numeric EXP required')
    predictions=[dict(c,public_exp=e,candidates={})for c,e in zip(cases,exps)]
    coefficient_probes=[]
    for count in range(min_count,max_count+1):
        selected=[i for i,c in enumerate(cases)if c['count']==count]
        inputs=np.array([[int(cases[i]['mean'],16),int(exps[i]['bits'],16)]for i in selected],dtype=np.uint64)
        for association in (('factorial_coefficients',)if polynomial else range(3)):
            graph=integer_gamma_polynomial_graph(count,association)if polynomial else integer_gamma_horner_graph(count,association)
            for model,cw in [('strict',0x123f),('stored_pc64',0x133f),('resident_pc64',0x133f)]:
                values=native.grid(graph,inputs,cw,0 if model=='resident_pc64'else(1<<graph['sites'])-1)[:,0]
                for i,q in zip(selected,values):predictions[i]['candidates'][f'{association}/{model}']=f'{int(q):016x}'
            if polynomial and count<max_count and any(cases[i]['scope']=='interior'for i in selected):
                from plateau_probe import probe
                scoped=[j for j,i in enumerate(selected)if cases[i]['scope']=='interior']
                shape=dict(graph=graph,constants={f'F{k}':graph['coefficient_nodes'][f'F{k}']for k in((0,1)if count==1 else(count,))},
                           cw=0x133f,mask=(1<<graph['sites'])-1)
                coefficient_probes.extend(dict(c,count=count)for c in probe(shape,inputs[scoped],native,keep=6))
                frozen[str(ROOT/'racer/plateau_probe.py')]=sha(ROOT/'racer/plateau_probe.py')
    write(str(output)+'.predictions.json',dict(format='dnacalc-integer-gamma-sum-validation-predictions/1',
        records=predictions,coefficient_probes=coefficient_probes,receipts=service.receipts.copy(),source_bindings=frozen))
    poisson=service.ask_rows('POISSON.DIST',[[client.number(float(c['count'])),client.number(c['mean']),client.boolean(True)]for c in cases])
    chi=service.ask_rows('CHIDIST',[[client.number(2.*value(c['mean'])),client.number(2.*c['shape'])]for c in cases])
    gamma=service.ask_rows('GAMMA.DIST',[[client.number(c['mean']),client.number(c['shape']),client.number(1.),client.boolean(True)]for c in cases])
    counts=Counter();contradictions={}
    for row,p,ch,ga in zip(predictions,poisson,chi,gamma):
        row['answers']=dict(POISSON=p,CHIDIST=ch,GAMMA=ga)
        if any(x['type']!='number'for x in(p,ch,ga)):raise ValueError('finite gamma aliases required')
        if row['scope']!='interior':continue
        counts['interior_rows']+=1;counts['poisson_equals_chidist']+=p==ch
        for program,prediction in row['candidates'].items():
            counts[program]+=p==dict(type='number',bits=prediction)
            if p!=dict(type='number',bits=prediction)and program not in contradictions:
                contradictions[program]=dict(count=row['count'],mean=row['mean'],expected=p,prediction=prediction)
        q=row['candidates'][('factorial_coefficients'if polynomial else'0')+'/stored_pc64'];pv=value(q)
        counts['gamma_complement_of_poisson']+=ga==dict(type='number',bits=bits(.5+(.5-value(p['bits']))))
        counts['gamma_horner_complement']+=ga==dict(type='number',bits=bits(.5+(.5-pv)))
    if frozen!={p:sha(p)for p in frozen}:raise ValueError('gamma validation source changed')
    observed={(row['count'],row['mean']):row['answers']['POISSON']for row in predictions}
    for c in coefficient_probes:
        for neighbor in c['neighbors']:
            samples=neighbor['examples'];neighbor['baseline_agrees']=sum(observed[c['count'],r['inputs'][0]]==dict(type='number',bits=r['baseline'])for r in samples)
            neighbor['alternative_agrees']=sum(observed[c['count'],r['inputs'][0]]==dict(type='number',bits=r['alternative'])for r in samples)
    result=dict(format='dnacalc-integer-gamma-sum-validation/1',complete=True,records=predictions,
        counts=counts,coefficient_probes=coefficient_probes,first_contradictions=contradictions,receipts=service.receipts,source_bindings=frozen,
        limits='Fresh validation of the finite polynomial in the declared integer shape region. Edge and shapes above20 answers are retained separately, never counted as interior passes. Separately observed EXP is a conditional binary64 leaf. Multiple equal associations remain undetermined until new questions distinguish them; this does not identify the general gamma path.')
    write(output,result);print(json.dumps(counts),flush=True)



def gamma_inverse_plateaus(checkout,cache,output,seed=35103):
    """Test the shared-forward lattice-ceiling premise on fresh quantiles."""
    import random
    rng=random.Random(seed);client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache)
    rows=[]
    for shape in (2.,3.,5.,10.,20.):
        for _ in range(16):rows.append(dict(p=bits(.55+.4*rng.random()),shape=bits(shape),scale=bits(1.)))
    rows += [dict(p=bits(.75),shape=bits(5.),scale=bits(1.))]
    source=sha(__file__)
    write(str(output)+'.questions.json',dict(rows=rows,seed=seed,source_sha256=source,
        selection='Five fixed integer shapes,16 fresh probability inputs each, plus the known lost Gamma inverse control. No quantile output selected by agreement.'))
    inverse=service.ask_rows('GAMMA.INV',[[client.number(r['p']),client.number(r['shape']),client.number(r['scale'])]for r in rows])
    questions=[]
    for r,a in zip(rows,inverse):
        if a['type']!='number':raise ValueError('finite inverse required')
        center=int(a['bits'],16)
        for offset in range(-16,17):questions.append([client.number(f'{center+offset:016x}'),client.number(r['shape']),client.number(r['scale']),client.boolean(True)])
    forward=service.ask_rows('GAMMA.DIST',questions)
    counts=Counter();records=[]
    for i,(r,a)in enumerate(zip(rows,inverse)):
        answers=forward[i*33:(i+1)*33]
        if any(x['type']!='number'for x in answers):raise ValueError('finite forward required')
        values=[value(x['bits'])for x in answers];target=value(r['p'])
        counts['rows']+=1;counts['returned_forward_equals_p']+=values[16]==target
        counts['lattice_ceiling_at_return']+=values[15]<target<=values[16]
        crossings=[j-16 for j,v in enumerate(values)if v>=target and(j==0 or values[j-1]<target)]
        records.append(dict(inputs=r,inverse=a,forward=answers,first_ge_offset=crossings[0]if crossings else None,
                            nonmonotone=any(x>y for x,y in zip(values,values[1:]))))
    if source!=sha(__file__):raise ValueError('inverse plateau source changed')
    result=dict(format='dnacalc-gamma-inverse-forward-plateaus/1',complete=True,counts=counts,
        records=records,receipts=service.receipts,source_sha256=source,
        limits='Tests the existing inverse premise that the result is the least binary64 x whose own public CDF reaches p. A contradiction refutes that premise, not the newly identified forward polynomial. No inverse rule is selected by a match count and no private CDF is inferred solely from these public outputs.')
    write(output,result);print(json.dumps(counts),flush=True)

def integer_gamma_polynomial_graph(count,form):
    """Fixed factorial-polynomial trees suggested by new-count contradictions."""
    g=Builder(2)
    if form in('factorial_coefficients','scaled_integer_coefficients'):
        scale=math.factorial(count)if form.startswith('scaled')else 1
        coefficients=[g.constant(bits(float(Q(scale,math.factorial(i)))))for i in range(count+1)]
        coefficient_nodes={f'F{i}':node for i,node in enumerate(coefficients)}
        s=coefficients[-1]
        for c in coefficients[-2::-1]:s=g.op('add',g.op('mul',s,0),c)
        if scale!=1:s=g.op('div',s,g.constant(bits(float(scale))))
    elif form in('reciprocal_after_product','reciprocal_before_product','reciprocal_on_x'):
        s=g.one
        for k in range(count,0,-1):
            c=g.constant(bits(1./k))
            if form=='reciprocal_after_product':t=g.op('mul',g.op('mul',s,0),c)
            elif form=='reciprocal_before_product':t=g.op('mul',g.op('mul',s,c),0)
            else:t=g.op('mul',s,g.op('mul',0,c))
            s=g.op('add',g.one,t)
    elif form=='coefficient_power_sum':
        s=g.one;p=g.one
        for k in range(1,count+1):
            p=g.op('mul',p,0);c=g.constant(bits(float(Q(1,math.factorial(k)))))
            s=g.op('add',s,g.op('mul',p,c))
    else:raise ValueError('unknown integer polynomial form')
    final=g.op('mul',1,s);graph=g.graph(final,final)
    graph['polynomial_store_mask']=(1<<g.nodes[s]['site'])|(1<<g.nodes[final]['site'])
    if form in('factorial_coefficients','scaled_integer_coefficients'):graph['coefficient_nodes']=coefficient_nodes
    return graph


def integer_gamma_polynomial_study(study,library,output):
    """Refute whole fixed polynomial trees; new-count data are now discovery."""
    import numpy as np
    doc=json.loads(Path(study).read_text());native=Native(library)
    rows=[r for r in doc['records']if r['scope']=='interior']
    if not rows:raise ValueError('empty polynomial scope')
    profiles=[]
    forms=('factorial_coefficients','scaled_integer_coefficients','reciprocal_after_product',
           'reciprocal_before_product','reciprocal_on_x','coefficient_power_sum')
    for form in forms:
        for model in('strict','stored_pc64','resident_pc64','resident_polynomial'):
            checked=0;first=None;mismatches=Counter()
            for count in sorted({r['count']for r in rows}):
                selected=[r for r in rows if r['count']==count]
                graph=integer_gamma_polynomial_graph(count,form)
                inputs=np.array([[int(r['mean'],16),int(r['public_exp']['bits'],16)]for r in selected],dtype=np.uint64)
                mask=0 if model=='resident_pc64'else graph['polynomial_store_mask']if model=='resident_polynomial'else(1<<graph['sites'])-1
                results=native.grid(graph,inputs,0x123f if model=='strict'else 0x133f,mask)[:,0]
                for row,got in zip(selected,results):
                    checked+=1;wanted=row['answers']['POISSON']['bits']
                    if int(got)!=int(wanted,16):
                        mismatches[str(count)]+=1
                        if first is None:first=dict(count=count,mean=row['mean'],expected=wanted,prediction=f'{int(got):016x}')
            profiles.append(dict(form=form,model=model,rows=checked,first_contradiction=first,mismatches_by_count=mismatches))
    result=dict(format='dnacalc-integer-gamma-polynomial-study/1',complete=True,profiles=profiles,
        survivors=[p for p in profiles if p['first_contradiction']is None],
        source_bindings={str(p):sha(p)for p in(Path(__file__),Path(study),Path(library))},
        limits='Six fixed mathematical factorial-polynomial trees, each with four coherent store policies, on all new-count observations. No coefficient adjusted and no count-specific selector. All answers are discovery after the divided-Horner failure; any survivor needs genuinely fresh validation.')
    write(output,result);print(json.dumps(dict(profiles=profiles)),flush=True)


def binomial_study(checkout,output,cache,offline=False,count=512):
    """Small fixed textbook trees, before introducing a distribution kernel.

    Complement, multiplication and addition use one coherent declared storage
    policy. Public POWER is a separate measured leaf, never an assumed private
    implementation. No tree is selected by its match count.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    centers=(.13,.3,.7)
    xs=[f'{b:016x}'for x in centers for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    jobs=((2,0,False),(2,1,False),(2,2,False),(2,1,True),(5,0,False),(5,5,False))
    questions=[[client.number(float(k)),client.number(float(n)),client.number(x),client.boolean(cdf)]
               for n,k,cdf in jobs for x in xs]
    answers=service.ask('BINOM.DIST',questions)['rows']
    # The effective complement is the published binary64 subtraction here;
    # the resident-complement candidate below keeps its own extended result.
    power_questions=[[client.number(bits(float(BINARY64.round(1-Q(value(x))))if complement else value(x))),client.number(float(n))]
                     for n,complement in ((2,True),(2,False),(5,True),(5,False))for x in xs]
    powers=service.ask('POWER',power_questions)['rows']
    records=[]
    for j,(n,k,cdf)in enumerate(jobs):
        shapes=('product_loop','square_then_multiply')if n==5 else ('direct_polynomial','sum_two_masses')if cdf else ('direct_polynomial',)
        profiles=[]
        for shape in shapes:
            for model in ('binary64','stored_pc64','resident_pc64'):
                op=lambda q:BINARY64.round(q)if model=='binary64'else BINARY64.round(EXTENDED64.round(q))if model=='stored_pc64'else EXTENDED64.round(q)
                tally=Counter();first=None
                for x,a in zip(xs,answers[j*len(xs):(j+1)*len(xs)]):
                    p=Q(value(x));q=op(1-p)
                    if n==5:
                        base=q if k==0 else p
                        if shape=='product_loop':
                            result=Q(1)
                            for _ in range(5):result=op(result*base)
                        else:
                            square=op(base*base);result=op(op(square*square)*base)
                    elif cdf:
                        if shape=='direct_polynomial':result=op(1-op(p*p))
                        else:result=op(op(q*q)+op(op(2*p)*q))
                    elif k==0:result=op(q*q)
                    elif k==1:result=op(op(2*p)*q)
                    else:result=op(p*p)
                    prediction=bits(float(BINARY64.round(result)))
                    key='matched'if a==dict(type='number',bits=prediction)else'mismatched';tally[key]+=1
                    if key=='mismatched'and first is None:first=dict(input=x,excel=a,candidate=prediction)
                profiles.append(dict(shape=shape,model=model,counts=tally,first_contradiction=first,
                                     consistent=tally['matched']==len(xs)))
        if not cdf and k in (0,n):
            leaf_index={(2,0):0,(2,2):1,(5,0):2,(5,5):3}[n,k]
            predicted=powers[leaf_index*len(xs):(leaf_index+1)*len(xs)]
            target=answers[j*len(xs):(j+1)*len(xs)]
            tally=Counter('matched'if p==a else'mismatched'for p,a in zip(predicted,target))
            first=next((dict(input=x,excel=a,candidate=p)for x,p,a in zip(xs,predicted,target)if p!=a),None)
            profiles.append(dict(shape='public_POWER_of_published_base',model='measured_binary64_leaf',counts=tally,
                                 first_contradiction=first,consistent=tally['matched']==len(xs)))
        records.append(dict(trials=n,successes=k,cumulative=cdf,profiles=profiles))
        print('BINOM',n,k,cdf,[(p['shape'],p['model'],dict(p['counts']))for p in profiles],flush=True)
    # The next textbook family is a logarithm followed by an exponential.
    # This is declared for every job, irrespective of earlier model counts.
    log_questions=[[client.number(bits(float(BINARY64.round(1-Q(value(x))))if complement else value(x)))]
                   for complement in (True,False)for x in xs]+[[client.number(2.)]]
    logarithms=service.ask('LN',log_questions)['rows']
    if any(a['type']!='number'for a in logarithms):raise ValueError('non-numeric binomial logarithm control')
    logq=[Q(value(a['bits']))for a in logarithms[:len(xs)]]
    logp=[Q(value(a['bits']))for a in logarithms[len(xs):2*len(xs)]]
    log2=Q(value(logarithms[-1]['bits']))
    complement_logs={'public LN(binary64(1-p))':logq}
    ideal={name:[]for name in ('mathematical LOG1P(-p) single','mathematical LOG1P(-p) double','mathematical LOG1P(-p) resident64')}
    for x in xs:
        lo,hi=elementary_bounds('LOG1P',[client.number(-value(x))])
        for name,rounder in zip(ideal,(BINARY64.round,lambda q:BINARY64.round(EXTENDED64.round(q)),EXTENDED64.round)):
            lower,upper=rounder(lo),rounder(hi)
            if lower!=upper:raise ValueError('private LOG1P candidate not certified by the enclosure')
            ideal[name].append(lower)
    complement_logs.update(ideal)
    exponents={};programs=[]
    for j,(n,k,cdf)in enumerate(jobs):
        for leaf_name,leaf_values in complement_logs.items():
            if leaf_name!='public LN(binary64(1-p))'and(cdf or k==n):continue
            for model in ('binary64','stored_pc64','resident_pc64'):
                op=lambda q:BINARY64.round(q)if model=='binary64'else BINARY64.round(EXTENDED64.round(q))if model=='stored_pc64'else EXTENDED64.round(q)
                arguments=[]
                for lp,lq in zip(logp,leaf_values):
                    t=op(2*lp)if cdf else op(n*lq)if k==0 else op(n*lp)if k==n else op(op(log2+lp)+lq)
                    raw=bits(float(BINARY64.round(t)));exponents.setdefault(raw,len(exponents));arguments.append(raw)
                programs.append((j,model,arguments,leaf_name))
    unique=list(exponents);exponential=[]
    for first in range(0,len(unique),8192):
        exponential.extend(service.ask('EXP',[[client.number(x)]for x in unique[first:first+8192]])['rows'])
    leaves=dict(zip(unique,exponential));log_profiles=[]
    for j,model,arguments,leaf_name in programs:
        n,k,cdf=jobs[j];tally=Counter();first=None
        op=lambda q:BINARY64.round(q)if model=='binary64'else BINARY64.round(EXTENDED64.round(q))if model=='stored_pc64'else EXTENDED64.round(q)
        for x,argument,a in zip(xs,arguments,answers[j*len(xs):(j+1)*len(xs)]):
            predicted=leaves[argument]
            if cdf and predicted['type']=='number':predicted=dict(type='number',bits=bits(float(BINARY64.round(op(1-Q(value(predicted['bits'])))))))
            key='matched'if predicted==a else'mismatched';tally[key]+=1
            if key=='mismatched'and first is None:first=dict(input=x,excel=a,candidate=predicted,exponent=argument)
        log_profiles.append(dict(trials=n,successes=k,cumulative=cdf,model=model,complement_logarithm=leaf_name,counts=tally,
            first_contradiction=first,consistent=tally['matched']==len(xs),
            shape='1-EXP(2*public LN(p))'if cdf else'EXP(n*complement_logarithm)'if k==0 else'EXP(n*public LN(p))'if k==n else'EXP((public LN(2)+public LN(p))+complement_logarithm)'))
        print('BINOM log',n,k,cdf,leaf_name,model,dict(tally),flush=True)
    result=dict(format='dnacalc-binomial-textbook-consumers/1',complete=True,
        program_order='Direct polynomial, plain product loop, square then multiply; binary64, storedPC64, residentPC64; public POWER control. Fixed structural order.',
        centers=centers,consecutive_inputs_per_center=count,records=records,log_profiles=log_profiles,receipts=service.receipts,
        limitations='Only n=2 and n=5 tail masses and n=2 one-count CDF. Public POWER and LN controls use a binary64 complement. Polynomial resident arithmetic retains the complement until final publication. Log trees use public binary64 LN(p) and EXP leaves, with a binary64 exponent at the service boundary. Mathematical LOG1P(-p) leaves are labelled conditional candidates, never private Excel observations. A surviving tree is conditional on this scope; no hidden general kernel or branch is identified. No coefficients or thresholds fit.')
    write(output,result)
    return result


def log_consumer_study(checkout,output,cache,offline=False,count=4096):
    """Textbook change-of-base consumers with independently measured LN leaves.

    Default decimal base first, then binary, ternary and the binary64 e anchor.
    This structural order is fixed before measurement, irrespective of counts.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    centers=(1.3,1.625);bases=(10.,2.,3.,math.e)
    xs=[f'{b:016x}'for x in centers for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    logarithms=service.ask('LN',[[client.number(x)]for x in xs])['rows']
    denominators=service.ask('LN',[[client.number(b)]for b in bases])['rows']
    if any(a['type']!='number'for a in logarithms+denominators):raise ValueError('non-numeric logarithm control')
    records=[]
    for base,denominator in zip(bases,denominators):
        target=service.ask('LOG',[[client.number(x),client.number(base)]for x in xs])['rows']
        profiles=[];division_predictions={}
        for shape in ('divide','reciprocal_then_multiply'):
            for model in ('binary64','stored_pc64','resident_pc64'):
                op=lambda q:BINARY64.round(q)if model=='binary64'else BINARY64.round(EXTENDED64.round(q))if model=='stored_pc64'else EXTENDED64.round(q)
                tally=Counter();first=None;predictions=[];d=Q(value(denominator['bits']))
                for x,leaf,a in zip(xs,logarithms,target):
                    numerator=Q(value(leaf['bits']))
                    q=op(numerator/d)if shape=='divide'else op(numerator*op(1/d))
                    predicted=bits(float(BINARY64.round(q)));predictions.append(predicted)
                    key='matched'if a==dict(type='number',bits=predicted)else'mismatched';tally[key]+=1
                    if key=='mismatched'and first is None:first=dict(input=x,excel=a,candidate=predicted)
                profiles.append(dict(shape=shape,model=model,counts=tally,first_contradiction=first,
                                     consistent=tally['matched']==len(xs)))
                if shape=='divide':division_predictions[model]=predictions
        discriminants=Counter()
        for a,s,d in zip(target,division_predictions['binary64'],division_predictions['stored_pc64']):
            if s==d:discriminants['indistinguishable']+=1
            else:discriminants['single_only'if a==dict(type='number',bits=s)else'double_only'if a==dict(type='number',bits=d)else'neither']+=1
        alias=None
        if base==10.:
            values=service.ask('LOG10',[[client.number(x)]for x in xs])['rows']
            alias=dict(shape='public LOG10(x)',counts=Counter('matched'if a==b else'mismatched'for a,b in zip(values,target)))
        records.append(dict(base=bits(base),denominator=denominator,profiles=profiles,
                            division_rounding=discriminants,log10_alias=alias,
                            ideal_leaf=rounding_evidence('LOG',[[client.number(x),client.number(base)]for x in xs],target)))
        print('LOG consumers',base,[(p['shape'],p['model'],dict(p['counts']))for p in profiles],dict(discriminants),flush=True)
        write(output,dict(complete=False,records=records,receipts=service.receipts,requested_bases=list(bases)))
    result=dict(format='dnacalc-log-change-of-base-consumers/1',complete=True,records=records,receipts=service.receipts,
        centers=centers,consecutive_inputs_per_center=count,
        limitations='Public LN outputs are measured binary64 leaves. Agreement identifies only the stated consumer on these windows and bases; private extended leaves and general store counts remain open. Equivalent arithmetic labels are not independent rounding evidence. No coefficient, threshold or selector is fit, and no candidate ranked by matches.')
    write(output,result)
    return result


def primitive_study(checkout,output,cache,offline=False,count=256):
    """Public leaves and small private consumers, in fixed source-code order.

    EXPON density and zero-count Poisson expose an exponential without any
    intervening arithmetic. All-success binomial exposes a power consumer.
    Kahan's expm1 identity is measured independently of a distribution body.
    These are scoped observations, never a count-based choice of a kernel.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    def window(centers):
        return [f'{b:016x}'for x in centers for b in
                range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    def ask(name,rows):
        answers=[]
        for start in range(0,len(rows),8192):
            answers.extend(service.ask(name,rows[start:start+8192])['rows'])
        return answers
    def profile(label,predictions,answers,inputs):
        tally=Counter();contradictions=[]
        for x,p,a in zip(inputs,predictions,answers):
            key='matched'if p==a else'mismatched';tally[key]+=1
            if p!=a:contradictions.append(dict(input=x,excel=a,candidate=p))
        return dict(program=label,counts=tally,consistent=tally['matched']==len(answers),
                    first_contradiction=next(iter(contradictions),None))
    def number(q):return dict(type='number',bits=bits(float(q)))
    def op(q,model):
        return BINARY64.round(q)if model=='binary64'else BINARY64.round(EXTENDED64.round(q))
    result=dict(format='dnacalc-private-primitive-consumers/1',complete=False,
                program_order='Published scalar leaves; zero-count and density identities; all-success EXP(n*LN(p)); Kahan correction, product first then division first. Binary64 then storedPC64. Fixed code-shape order, never match-count rank.',
                consecutive_inputs_per_center=count,leaves=[],exponential_consumers=[],powers=[],expm1=[])
    for name,centers in [('EXP',(-700.,-50.,-1.,-.7,.7,1.,50.,700.)),
                         ('LN',(.001,.13,.3,.7,.99,1.3,3.7,10.)),
                         ('SQRT',(.13,.7,1.3,3.7,10.))]:
        xs=window(centers);rows=[[client.number(x)]for x in xs]
        answers=ask(name,rows)
        result['leaves'].append(dict(function=name,centers=centers,
                                    reference=rounding_evidence(name,rows,answers)))
        print('primitive leaf',name,result['leaves'][-1]['reference'],flush=True)
        write(output,dict(result,receipts=service.receipts))
    xs=window((2**-30,.13,.7,1.,3.7,50.,700.))
    negative=[[client.number(bits(-value(x)))]for x in xs]
    exponential=ask('EXP',negative)
    for name,rows in [('EXPON.DIST',[[client.number(x),client.number(1.),client.boolean(False)]for x in xs]),
                      ('POISSON.DIST',[[client.number(0.),client.number(x),client.boolean(False)]for x in xs])]:
        answers=ask(name,rows)
        result['exponential_consumers'].append(profile(name+' = public EXP(-x)',exponential,answers,xs))
        print('private exponential',name,result['exponential_consumers'][-1],flush=True)
    ps=window((.001,.13,.3,.7,.99))
    logarithms=ask('LN',[[client.number(p)]for p in ps])
    if any(a['type']!='number'for a in logarithms):raise ValueError('non-numeric power logarithm')
    unique={};programs=[]
    for n in (1,2,3,5,10,32,100):
        answers=ask('BINOM.DIST',[[client.number(float(n)),client.number(float(n)),client.number(p),client.boolean(False)]for p in ps])
        for model in ('binary64','stored_pc64'):
            arguments=[bits(float(op(n*Q(value(a['bits'])),model)))for a in logarithms]
            for t in arguments:unique.setdefault(t,None)
            programs.append((n,model,arguments,answers))
    keys=list(unique)
    leaves=ask('EXP',[[client.number(t)]for t in keys]);lookup=dict(zip(keys,leaves))
    for n,model,arguments,answers in programs:
        predictions=[lookup[t]for t in arguments]
        row=profile('BINOM(n,n,p,FALSE) = EXP(n*LN(p)); '+model,predictions,answers,ps)
        row.update(trials=n,model=model);result['powers'].append(row)
        print('private power',n,model,dict(row['counts']),flush=True)
    # Only the small negative arm needs LN(u); large |t| uses the plain
    # difference, and u==1 must return t before logarithm/division.
    answers=ask('EXPON.DIST',[[client.number(x),client.number(1.),client.boolean(True)]for x in xs])
    log_inputs=sorted({a['bits']for a in exponential if a['type']=='number'and value(a['bits'])!=0.})
    logs=dict(zip(log_inputs,ask('LN',[[client.number(x)]for x in log_inputs])))
    for shape in ('product_then_divide','divide_then_product'):
        for model in ('binary64','stored_pc64'):
            predictions=[]
            for x,u in zip(xs,exponential):
                if u['type']!='number':raise ValueError('non-numeric expm1 exponential')
                t=-Q(value(x));u=Q(value(u['bits']))
                if u==1:z=t
                elif abs(t)>=1:z=op(u-1,model)
                else:
                    l=Q(value(logs[bits(float(u))]['bits']));d=op(u-1,model)
                    z=op(op(d*t,model)/l,model)if shape=='product_then_divide'else op(d*op(t/l,model),model)
                predictions.append(number(-z))
            row=profile('EXPON CDF = -Kahan_expm1(-x); '+shape+'; '+model,predictions,answers,xs)
            result['expm1'].append(row)
            print('private expm1',shape,model,dict(row['counts']),flush=True)
    result.update(complete=True,receipts=service.receipts,
        limitations='Scalar Value2 windows only. Public leaves have independent Excel readings and MPFR directed mathematical diagnostics; mathematical values never become Excel answers. Complete consumer equality is conditional on the measured binary64 leaf boundary. These observations do not recover resident extended private values, global branches, or a general algorithm. No coefficient, tolerance, threshold or selector fitted.')
    write(output,result)
    return result


def bisect_observed_edges(service, function, brackets, make_row, side):
    """Find adjacent false/true observations, retaining every original reply.

    Endpoints and every interior observation are checked; an indistinguishable
    answer stops that bracket. Uniqueness requires a separately justified
    monotone predicate. No unobserved midpoint is assigned a side.
    """
    states=[dict(b,low=int(b['low'],16),high=int(b['high'],16),trace=[])for b in brackets]
    for s in states:
        if not 0<=s['low']<s['high']<0x7ff0000000000000:
            raise ValueError('positive finite increasing edge bracket required')
    def query(points):
        rows=[make_row(s,f'{x:016x}')for s,x in points]
        answers=service.ask(function,rows)['rows']
        for (s,x),answer in zip(points,answers):
            result=side(s,x,answer)
            if result is not True and result is not False:
                raise ValueError('edge predicate cannot distinguish a measured input')
            s['trace'].append(dict(input=f'{x:016x}',answer=answer,side=result))
        return answers
    endpoints=[(s,x)for s in states for x in (s['low'],s['high'])]
    query(endpoints)
    for s in states:
        if s['trace'][0]['side'] or not s['trace'][1]['side']:
            raise ValueError('edge endpoints must be false then true')
    while any(s['high']-s['low']>1 for s in states):
        active=[s for s in states if s['high']-s['low']>1]
        points=[(s,(s['low']+s['high'])//2)for s in active]
        query(points)
        for s,x in points:
            s['high'if s['trace'][-1]['side']else'low']=x
    return [dict(s,low=f'{s["low"]:016x}',high=f'{s["high"]:016x}')for s in states]


def binomial_admission_study(checkout,output,cache,offline=False):
    """Isolate the count admission boundary before any mass primitive.

    The bracket comes from the R27 numeric/error observations, not a fitted
    cutoff. Endpoint probabilities and cumulative all-success calls make the
    answer independent of LN, EXP and the interior probability algorithm.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    probabilities=(0.,1.,1.-2.**-24)
    brackets=[dict(low=bits(2.**30),high=bits(2.**31),probability=bits(p),
                   cumulative=cumulative)
              for p in probabilities for cumulative in (False,True)]
    def row(s,b):
        return [client.number(b),client.number(b),client.number(s['probability']),
                client.boolean(s['cumulative'])]
    def side(s,x,a):
        if a['type']=='number':return False
        if a==dict(type='error',code='#NUM!'):return True
        return None
    edges=bisect_observed_edges(service,'BINOM.DIST',brackets,row,side)
    window=[]
    for edge in sorted({int(e['high'],16)for e in edges}):
        window.extend(f'{b:016x}'for b in range(edge-16,edge+16))
    controls=[]
    # Adjacent raw count operands separate truncation from nearest-integer
    # conversion. k=0/1/all-success and p=0/1 expose validation before leaves.
    for n in window+[bits(x)for x in (0.,1.,2.**32-1,2.**32,2.**40)]:
        for k in ('zero','one','all'):
            for p in probabilities:
                for cumulative in (False,True):
                    b=n if k=='all'else bits(0. if k=='zero'else 1.)
                    controls.append(dict(trials=n,successes=b,probability=bits(p),cumulative=cumulative))
    rows=[[client.number(r['successes']),client.number(r['trials']),
           client.number(r['probability']),client.boolean(r['cumulative'])]for r in controls]
    aliases={name:service.ask(name,rows)['rows']for name in ('BINOM.DIST','BINOMDIST')}
    result=dict(format='dnacalc-binomial-admission-edge/1',complete=True,edges=edges,
                controls=[dict(r,answers={name:a[i]for name,a in aliases.items()})for i,r in enumerate(controls)],
                receipts=service.receipts,
                limits='Measured count admission, separate from numerical primitive identification. Binary64 inputs are read back exactly. No candidate threshold chosen by pass counts; the adjacent transition is directly observed. Other discrete functions need their own admission observations.')
    write(output,result)
    print('BINOM admission',[(e['probability'],e['cumulative'],e['low'],e['high'])for e in edges],flush=True)
    return result


def primitive_stage_graph(kind, constant='rom', stored_reduction=False):
    """One declared ISA primitive; no outer consumer or store-mask search."""
    g=Builder(1)
    if kind in ('LN','LOG2'):
        c=(g.one if kind=='LOG2'else
           g.op('fldln2',store=False)if constant=='rom'else g.constant(bits(math.log(2.))))
        out=g.op('fyl2x',c,0,store=False)
    else:
        if kind=='EXP':
            c=g.op('fldl2e',store=False)if constant=='rom'else g.constant(bits(math.log2(math.e)))
            t=g.op('mul',0,c,store=stored_reduction)
        elif kind=='EXP2':t=0
        else:raise ValueError('unknown primitive stage')
        whole=g.op('frndint',t,store=False)
        fraction=g.op('sub',t,whole,store=False)
        small=g.op('f2xm1',fraction,store=False)
        mantissa=g.op('add',g.one,small,store=False)
        out=g.op('fscale',mantissa,whole,store=False)
    return g.graph(out,out)


def positive_extended_cell(native, graph, answer, center,decreasing=False,max_radius_bits=18,output_span=None):
    """Conditional inverse of one fixed positive extended-input primitive.

    Find both ends of an output's preimage on the PC64 input lattice. The
    local monotonicity of this fixed instruction chain is a declared premise;
    Excel is never queried with an inferred extended operand. No logarithm
    producer is assumed. Empty preimages are kept as contradictions. An
    optional exact output span constrains the resident PC64 result instead
    of its final binary64 publication, using the same inverse search.
    """
    from x87 import encode,decode
    significand=EXTENDED64.floor_coordinate(EXTENDED64.round(center))
    radius=1<<18
    lo,hi=significand-radius,significand+radius
    if not 0<lo<hi<EXTENDED64.last:
        raise ValueError('inverse bracket leaves positive finite extended values')
    def pack(sig):return encode(EXTENDED64.value(sig))
    def raw_result(sig):
        return native.evaluate(graph,[],0x133f,extended=pack(sig))[:10]
    def evaluate(sig):
        raw=raw_result(sig)
        return (native.publish(raw)if output_span is None else
                EXTENDED64.floor_coordinate(decode(raw)))
    direction=-1 if decreasing else 1
    if output_span is None:
        target_low=target_high=direction*int(answer,16)
    else:
        coordinates=EXTENDED64.coordinates(output_span)
        if coordinates is None:return None
        target_low,target_high=(coordinates if direction==1 else(-coordinates[1],-coordinates[0]))
    if not 18<=max_radius_bits<=62:raise ValueError('bounded inverse radius required')
    while not(direction*evaluate(lo)<target_low and target_high<direction*evaluate(hi)):
        if radius>=1<<max_radius_bits:raise ValueError('inverse bracket does not enclose requested output')
        radius*=2;lo,hi=significand-radius,significand+radius
        if not 0<lo<hi<EXTENDED64.last:raise ValueError('expanded inverse bracket leaves positive finite values')
    def first(target,strict):
        a,b=lo,hi
        while b-a>1:
            m=(a+b)//2;got=direction*evaluate(m)
            if got>target if strict else got>=target:b=m
            else:a=m
        return b
    lower,upper=first(target_low,False),first(target_high,True)-1
    if lower>upper:return None
    if not(target_low<=direction*evaluate(lower)<=target_high and target_low<=direction*evaluate(upper)<=target_high):
        raise ValueError('inverse endpoint does not publish requested output')
    result=dict(low=pack(lower).hex(),high=pack(upper).hex(),bracket_radius=radius)
    if output_span is None:
        result.update(predecessor_output=f'{evaluate(lower-1):016x}',successor_output=f'{evaluate(upper+1):016x}')
    else:
        result.update(predecessor_output_extended=raw_result(lower-1).hex(),
                      successor_output_extended=raw_result(upper+1).hex(),
                      output_span=output_span.json())
    return result


def resident_stage_study(checkout,output,cache,library,offline=False,count=256):
    """Isolate reduction precision and invert root observations one stage at a time.

    GEOMEAN(x,1,...,1) has an exact product. Arity2/4/8 makes division
    by the count exact in PC64 and removes reciprocal-rounding uncertainty.
    Shared-log interval intersections test an unknown resident logarithm,
    rather than enumerating complete logarithm/root programs.
    """
    from x87 import decode
    library=Path(library)
    if not library.exists():build_library(library)
    native=Native(library)
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    graphs={
        'LN_ROM80':primitive_stage_graph('LN'),
        'LN_constant53':primitive_stage_graph('LN','binary64'),
        'EXP_ROM80_resident_reduction':primitive_stage_graph('EXP'),
        'EXP_constant53_resident_reduction':primitive_stage_graph('EXP','binary64'),
        'EXP_ROM80_stored_reduction':primitive_stage_graph('EXP',stored_reduction=True),
        'EXP2':primitive_stage_graph('EXP2'),'LOG2':primitive_stage_graph('LOG2')}
    def window(centers):
        return [f'{b:016x}'for center in centers for b in
                range(int(bits(center),16)-count//2,int(bits(center),16)+count//2)]
    exp_inputs=window((.7,1.,50.,700.))
    log_inputs=window((1.3,3.7,7.,1000.))
    observations={};primitive_checks=[]
    for name,xs,programs in [('EXP',exp_inputs,list(graphs)[2:5]),('LN',log_inputs,list(graphs)[:2])]:
        answers=service.ask(name,[[client.number(x)]for x in xs])['rows']
        observations[name]=dict(inputs=xs,answers=answers)
        for program in programs:
            tally=Counter();first=None
            for x,a in zip(xs,answers):
                raw=native.evaluate(graphs[program],[int(x,16)],0x133f,mask=1 if program.endswith('stored_reduction')else 0)[:10]
                predicted=f'{native.publish(raw):016x}'
                key='matched'if a==dict(type='number',bits=predicted)else'mismatched'
                tally[key]+=1
                if key=='mismatched'and first is None:first=dict(input=x,excel=a,predicted=predicted,resident_output=raw.hex())
            primitive_checks.append(dict(program=program,counts=tally,first_contradiction=first))
    power_answers=service.ask('POWER',[[client.number(math.e),client.number(x)]for x in exp_inputs])['rows']
    power_counts=Counter('matched'if a==b else'mismatched'for a,b in zip(power_answers,observations['EXP']['answers']))
    power_routes={route:Counter()for route in ('integer','fractional')}
    power_contradictions=[]
    for x,a,b in zip(exp_inputs,power_answers,observations['EXP']['answers']):
        route='integer'if value(x).is_integer()else'fractional'
        power_routes[route]['matched'if a==b else'mismatched']+=1
        if a!=b:power_contradictions.append(dict(exponent=x,route=route,power=a,exp=b))
    # Natural and base2 argument inversions are distinct conditional primitive
    # studies. Neither is allowed to invent a per-row logarithm or correction.
    xs=window((1.3,7.,1000.))
    roots={n:service.ask('GEOMEAN',[[client.number(x)]+[client.number(1.)]*(n-1)for x in xs])['rows']for n in (2,4,8)}
    intervals=[];tallies={kind:Counter()for kind in ('EXP','EXP2')}
    for index,x in enumerate(xs):
        log=decode(native.evaluate(graphs['LN_ROM80'],[int(x,16)],0x133f)[:10])
        log2=decode(native.evaluate(graphs['LOG2'],[int(x,16)],0x133f)[:10])
        for kind,center in [('EXP',log),('EXP2',log2)]:
            graph=graphs['EXP_ROM80_resident_reduction'if kind=='EXP'else'EXP2']
            cells=[];lower=[];upper=[]
            for n in (2,4,8):
                a=roots[n][index]
                if a['type']!='number':raise ValueError('non-numeric root in positive stage window')
                cell=positive_extended_cell(native,graph,a['bits'],center/n)
                cells.append(dict(count=n,excel=a,cell=cell))
                if cell:
                    lower.append(decode(bytes.fromhex(cell['low']))*n)
                    upper.append(decode(bytes.fromhex(cell['high']))*n)
            common_low=max(lower)if len(lower)==3 else None
            common_high=min(upper)if len(upper)==3 else None
            feasible=common_low is not None and common_low<=common_high
            label='shared_log_possible'if feasible else'shared_log_impossible'
            tallies[kind][label]+=1
            candidate=center
            if feasible:
                tallies[kind]['registered_resident_log_in_intersection'if common_low<=candidate<=common_high else'registered_resident_log_excluded']+=1
                published=BINARY64.round(candidate)
                tallies[kind]['binary64_log_in_intersection'if common_low<=published<=common_high else'binary64_log_excluded']+=1
            intervals.append(dict(input=x,primitive=kind,cells=cells,shared_log_possible=feasible,
                                  common_log_interval=[str(common_low),str(common_high)]if feasible else None))
    result=dict(format='dnacalc-resident-primitive-stage-inverses/1',complete=True,
        primitive_checks=primitive_checks,power_e_equals_public_exp=power_counts,
        power_routes=power_routes,power_contradictions=power_contradictions,
        root_stage_counts=tallies,root_intervals=intervals,
        graphs=graphs,consecutive_inputs_per_center=count,receipts=service.receipts,
        native_bindings={str(p):sha(p)for p in (library,ROOT/'racer/x87_native.c',ROOT/'racer/x87_programs.py',ROOT/'checker/x87.py')},
        reduction_constants=dict(ln2_rom80='ac79cfd1f71772b1fe3f',log2e_rom80='bcf0175c293baab8ff3f',
                                 ln2_binary64=bits(math.log(2.)),log2e_binary64=bits(math.log2(math.e))),
        limitations='Positive normal scalar windows only. Public EXP/LN agreement is a conditional ISA-chain observation, not a global software implementation or microcode identification. Root inverse cells condition on one fixed monotone EXP or EXP2 primitive with resident PC64 input and no later correction. Empty shared-log intersections refute that factorization even with an arbitrary common logarithm. Nonempty intersections only bound a hidden stage; they never identify a logarithm polynomial/table or authorize fitted values. Arity2/4/8 eliminates reciprocal precision, not private route differences. No whole-program enumeration, count ranking or answer-dependent algorithm selection.')
    write(output,result)
    print('primitive stages',primitive_checks,dict(power_counts),{k:dict(v)for k,v in tallies.items()},flush=True)
    return result


def reduction_constant_study(checkout,output,cache,library,offline=False,count=256):
    """Intersect stage preimages to constrain one shared reduction constant.

    LN uses a fused native FYL2X primitive with an unknown PC64 coefficient.
    EXP uses one PC64 multiplication followed by the fixed EXP2 primitive.
    Constants are recovered as necessary intervals, then compared to the
    published x87 ROM values. No fitted coefficient is installed.
    """
    from exact import Span
    from x87 import decode,encode
    native=Native(library)
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    ln2=native.leaf('fldln2',bytes(10),bytes(10),0x133f)
    l2e=native.leaf('fldl2e',bytes(10),bytes(10),0x133f)
    ln_graph=Builder(2);out=ln_graph.op('fyl2x',1,0,store=False);ln_graph=ln_graph.graph(out,out)
    exp_graph=primitive_stage_graph('EXP2')
    class Coefficient:
        def __init__(self,x):self.x=x
        def evaluate(self,graph,inputs,cw,extended):return native.evaluate(graph,[self.x],cw,extended=extended)
        def publish(self,raw):return native.publish(raw)
    def neighbor(raw,delta):
        n=int.from_bytes(raw[:8],'little')+delta
        if not 1<<63<=n<1<<64:raise ValueError('constant inverse crosses extended binade')
        return n.to_bytes(8,'little')+raw[8:]
    records=[];families=[]
    for function,centers,coefficient in [('LN',(1.3,3.7,7.,1000.),ln2),('EXP',(.7,1.,50.,700.),l2e)]:
        xs=[f'{b:016x}'for center in centers for b in range(int(bits(center),16)-count//2,int(bits(center),16)+count//2)]
        answers=service.ask(function,[[client.number(x)]for x in xs])['rows']
        common=None;first=True
        for x,a in zip(xs,answers):
            if a['type']!='number':raise ValueError('constant probe must be positive numeric')
            q=Q(value(x));center=decode(coefficient)
            if function=='LN':
                cell=positive_extended_cell(Coefficient(int(x,16)),ln_graph,a['bits'],center)
            else:
                cell=positive_extended_cell(native,exp_graph,a['bits'],EXTENDED64.round(q*center))
            if cell is None:raise ValueError('empty stage preimage for reduction-constant control')
            lo,hi=(bytes.fromhex(cell[k])for k in ('low','high'))
            if function=='LN':
                domain=Span(decode(lo),decode(hi))
            else:
                # Invert the single RN64 multiply at both extreme input
                # cells. Exact midpoint ownership comes from tie parity.
                a0=(decode(neighbor(lo,-1))+decode(lo))/2/q
                b0=(decode(hi)+decode(neighbor(hi,1)))/2/q
                domain=Span(a0,b0,not(int.from_bytes(lo[:8],'little')&1),not(int.from_bytes(hi[:8],'little')&1))
            common=domain if first else common.meet(domain)if common is not None else None
            first=False
            records.append(dict(function=function,input=x,excel=a,cell=cell,coefficient_interval=domain.json()))
        if common is None:
            candidates=[];total=0
        else:
            lower=EXTENDED64.ceiling_coordinate(common.low)
            upper=EXTENDED64.floor_coordinate(common.high)
            if not common.left and EXTENDED64.value(lower)==common.low:lower+=1
            if not common.right and EXTENDED64.value(upper)==common.high:upper-=1
            total=max(0,upper-lower+1)
            candidates=[encode(EXTENDED64.value(k)).hex()for k in range(lower,min(upper+1,lower+16))]
        families.append(dict(function=function,observations=len(xs),
            common_coefficient_interval=common.json()if common is not None else None,
            pc64_candidate_count=total,first_pc64_candidates=candidates,
            published_rom=coefficient.hex(),rom_in_intersection=common is not None and common.contains(decode(coefficient))))
    result=dict(format='dnacalc-primitive-reduction-constant-preimages/1',complete=True,
        families=families,records=records,receipts=service.receipts,
        native_bindings={str(p):sha(p)for p in (Path(library),ROOT/'racer/x87_native.c',ROOT/'racer/x87_programs.py',ROOT/'checker/x87.py')},
        limits='Necessary coefficient intervals under the explicitly fixed positive native FYL2X and EXP2 primitives and local monotonicity. EXP assumes one resident RN64 multiply, not a binary64 spill. Each coefficient must be common to every input. Published ROM values are controls, not fitted implementations. Public scalar probes do not identify GEOMEAN/private helper code, a software polynomial/table, or global route boundaries. Reuses original readings; no additional calls are counted.')
    write(output,result);print('reduction constants',families,flush=True)
    return result


def refine_log_constant(checkout,prior,output,cache,library,offline=False,spread=False):
    """Ask only preselected coefficient discriminants in fixed input windows.

    The prior stage interval supplies two endpoint coefficients, not a fitted
    numerical implementation. Inputs are chosen where those endpoints publish
    different answers before any new Excel reading. The common constant must
    satisfy every old and new observation.
    """
    from exact import Span
    from x87 import decode,encode
    prior=Path(prior);old=json.loads(prior.read_text())
    family=next(r for r in old['families']if r['function']=='LN')
    candidates=family['first_pc64_candidates']
    if not 1<len(candidates)==family['pc64_candidate_count']<=16:
        raise ValueError('refinement needs the complete bounded prior PC64 interval')
    low,high=bytes.fromhex(candidates[0]),bytes.fromhex(candidates[-1])
    native=Native(library);client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    centers=tuple(1.+i/256 for i in range(1,256))if spread else(1.3,3.7,7.,1000.)
    count=512 if spread else 65536
    selected=[]
    for center in centers:
        start=int(bits(center),16)-count//2
        for b in range(start,start+count):
            x=encode(Q(value(f'{b:016x}')))
            a,z=(native.publish(native.leaf('fyl2x',c,x,0x133f))for c in (low,high))
            if a!=z:selected.append(dict(input=f'{b:016x}',lower_prediction=f'{a:016x}',upper_prediction=f'{z:016x}'))
    # This frozen selection is written before a new observation is requested.
    selection=dict(format='dnacalc-log-constant-refinement-selection/1',prior_sha256=sha(prior),
        centers=centers,consecutive_inputs_per_center=count,coefficient_endpoints=[low.hex(),high.hex()],
        selected=selected,selection='All endpoint disagreements in the fixed consecutive windows; chosen before new answers.')
    selection_path=Path(str(output)+'.selection.json');write(selection_path,selection)
    answers=[]
    for start in range(0,len(selected),8192):
        answers.extend(service.ask('LN',[[client.number(r['input'])]for r in selected[start:start+8192]])['rows'])
    g=Builder(2);out=g.op('fyl2x',1,0,store=False);graph=g.graph(out,out)
    class Coefficient:
        def __init__(self,x):self.x=x
        def evaluate(self,graph,inputs,cw,extended):return native.evaluate(graph,[self.x],cw,extended=extended)
        def publish(self,raw):return native.publish(raw)
    original=family['common_coefficient_interval']
    common=Span(Q(*map(int,original['lower'])),Q(*map(int,original['upper'])),original['lower_closed'],original['upper_closed'])
    records=[]
    for r,a in zip(selected,answers):
        if a['type']!='number':raise ValueError('nonnumeric LN refinement answer')
        cell=positive_extended_cell(Coefficient(int(r['input'],16)),graph,a['bits'],decode(low))
        if cell is None:raise ValueError('empty coefficient cell')
        domain=Span(decode(bytes.fromhex(cell['low'])),decode(bytes.fromhex(cell['high'])))
        common=common.meet(domain)if common is not None else None
        records.append(dict(r,excel=a,cell=cell))
    remaining=[c for c in candidates if common is not None and common.contains(decode(bytes.fromhex(c)))]
    result=dict(format='dnacalc-log-constant-refinement/1',complete=True,prior_sha256=sha(prior),
        selection_sha256=sha(selection_path),original_candidates=candidates,remaining_pc64_candidates=remaining,
        common_coefficient_interval=common.json()if common is not None else None,
        records=records,receipts=service.receipts,
        limits='Conditional coefficient identification for the declared native FYL2X leaf only. No table/polynomial or private Excel logarithm is identified; no coefficient is fitted into product code. Selected questions are new observed coverage, not kernel gains.')
    write(output,result);print('LN coefficient refinement',len(selected),'questions; remaining',remaining,flush=True)
    return result


def root_anchor_study(checkout,output,cache,offline=False):
    """Exact product/root anchors separate the private root from public pow.

    On product2^(n*j), the mathematical nth root is the exactly representable
    2^j. No logarithm, exponential approximation or recurrence is needed for
    that control. Neighbor exponents retain their directed mathematical bounds.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    keys=[]
    for n in (*range(2,13),16,32):
        for j in (-100,-10,-1,0,1,10,100):
            if -1022<=n*j<=1023:keys.append((n,n*j,'exact_integer_log2_root'))
        for k in range(-n+1,n):
            if k:keys.append((n,k,'fractional_log2_root'))
    questions=[[client.number(math.ldexp(1.,k))]+[client.number(1.)]*(n-1)for n,k,_ in keys]
    answers=service.ask('GEOMEAN',questions)['rows']
    powers=service.ask('POWER',[[row[0],client.number(1./n)]for (n,_,_),row in zip(keys,questions)])['rows']
    records=[];counts={}
    for (n,k,role),row,a,p in zip(keys,questions,answers,powers):
        ideal=publication(elementary_bounds('GEOMEAN',row))
        if ideal is None:raise ValueError('root anchor reference not enclosed')
        target=dict(type='number',bits=ideal[0])
        tally=counts.setdefault(str(n),Counter())
        tally[role+'_matched'if a==target else role+'_mismatched']+=1
        tally['public_POWER_equal'if a==p else'public_POWER_different']+=1
        records.append(dict(count=n,product_power=k,role=role,input=row[0]['bits'],
                            excel=a,public_power=p,mathematical_single=ideal[0],mathematical_double=ideal[1]))
    result=dict(format='dnacalc-geomean-exact-power-anchors/1',complete=True,records=records,counts=counts,
        receipts=service.receipts,
        limits='Exact unit-factor products. Integer log2 roots are exact mathematical controls; noninteger roots are directed MPFR controls, never Excel answers. Public POWER uses a stored binary64 reciprocal and its own wrapper. Differences isolate routes but do not identify a private polynomial/table, a Newton loop or its stopping rule. No fitted coefficient, tolerance or selector.')
    write(output,result);print('root anchors',{n:dict(c)for n,c in counts.items()},flush=True)
    return result


def poisson_log_difference_study(checkout,original,output,cache,library,offline=False):
    """Invert one subtraction before the resident PMF1 exponential argument."""
    from exact import Span
    from x87 import decode,load64
    prior=json.loads(Path(original).read_text())
    if not prior.get('complete'):raise ValueError('incomplete exponential preimages')
    rows=[r for r in prior['records']if r['stage']=='zero_to_one']
    inputs=[r['input']for r in rows]
    if len(set(inputs))!=len(inputs):raise ValueError('repeated original PMF1 input')
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    logs=dict(zip(inputs,service.ask('LN',[[client.number(x)]for x in inputs])['rows']))
    native=Native(library);coefficient=native.leaf('fldln2',bytes(10),bytes(10),0x133f)
    records=[];counts=Counter()
    for r in rows:
        x=r['input'];mean=Q(value(x));cell=r['cell']
        lo,hi=(decode(bytes.fromhex(cell[k]))for k in ('low','high'))
        argument=Span(-hi,-lo)
        # t=RN64(L-mean): take every rounding cell of the permitted t lattice,
        # then undo exactly one subtraction. Midpoint ownership is retained.
        domain=EXTENDED64.pullback(argument).linear(Q(1),mean)
        observed=logs[x]
        if observed['type']!='number':raise ValueError('positive-mean public log must be numeric')
        resident=native.leaf('fyl2x',coefficient,load64(int(x,16)),0x133f)
        public_possible=domain.contains(Q(value(observed['bits'])))
        resident_possible=domain.contains(decode(resident))
        pc64=EXTENDED64.coordinates(domain);binary64=BINARY64.coordinates(domain)
        counts['public_LN_possible'if public_possible else'public_LN_excluded']+=1
        counts['resident_LN_possible'if resident_possible else'resident_LN_excluded']+=1
        counts['binary64_log_possible'if binary64 is not None else'binary64_log_impossible']+=1
        records.append(dict(index=r['index'],input=x,excel=r['excel'],cell=cell,
            private_log_interval=domain.json(),public_LN=observed,resident_LN=resident.hex(),
            public_LN_possible=public_possible,resident_LN_possible=resident_possible,
            pc64_log_count=0 if pc64 is None else pc64[1]-pc64[0]+1,
            binary64_log_count=0 if binary64 is None else binary64[1]-binary64[0]+1))
    result=dict(format='dnacalc-poisson-log-difference-preimages/1',complete=True,
        original_sha256=sha(original),original_receipts=prior['original_receipts'],records=records,counts=counts,receipts=service.receipts,
        native_bindings={str(p):sha(p)for p in (Path(library),ROOT/'racer/x87_native.c',ROOT/'racer/x87_programs.py')},
        limits='Necessary private-log intervals for one declared t=RN64(L-mean) stage and the earlier conditional resident public EXP primitive. Public LN is independently observed; resident FYL2X is an explicit comparison. No fitting, private logarithm identification or full PMF algorithm adoption follows from nonempty intervals.')
    write(output,result);print(json.dumps(counts));return result


def poisson_log_primitive_contrasts(original,output,library):
    """Test fixed native LN leaves against previously derived necessary cells.

    This isolates the coefficient precision, operation precision and return
    publication. It does not select the least-bad leaf or infer a private LN.
    """
    from exact import Span
    from x87 import decode,load64
    path=Path(original)
    prior=json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes())
    if prior.get('complete') is not True:raise ValueError('incomplete original log constraints')
    native=Native(library);rom=native.leaf('fldln2',bytes(10),bytes(10),0x133f)
    published=load64(native.publish(rom))
    leaves=[(name,coefficient,word,stored)
            for name,coefficient in [('rom_ln2',rom),('binary64_ln2',published)]
            for word in (0x123f,0x133f) for stored in (False,True)]
    counts={f'{name}/pc{53 if word==0x123f else 64}/{"stored" if stored else "resident"}':Counter()
            for name,_,word,stored in leaves}
    records=[]
    for row in prior['records']:
        d=row['private_log_interval']
        cell=Span(Q(*map(int,d['lower'])),Q(*map(int,d['upper'])),d['lower_closed'],d['upper_closed'])
        predictions={}
        for name,coefficient,word,stored in leaves:
            key=f'{name}/pc{53 if word==0x123f else 64}/{"stored" if stored else "resident"}'
            raw=native.leaf('fyl2x',coefficient,load64(int(row['input'],16)),word)
            if stored:raw=load64(native.publish(raw))
            possible=cell.contains(decode(raw));counts[key]['possible' if possible else 'excluded']+=1
            predictions[key]=dict(native=raw.hex(),possible=possible)
        records.append(dict(input=row['input'],excel=row['excel'],private_log_interval=d,predictions=predictions))
    result=dict(format='dnacalc-poisson-native-log-contrasts/1',complete=True,
        source_bindings={str(path.resolve()):sha(path),str(Path(library).resolve()):sha(library),
                         str(ROOT/'racer/x87_native.c'):sha(ROOT/'racer/x87_native.c'),
                         str(ROOT/'racer/x87_programs.py'):sha(ROOT/'racer/x87_programs.py')},
        original_receipts=prior['original_receipts'],counts=counts,records=records,
        limits='All eight leaves are fixed before evaluation. A single exclusion refutes that leaf '
        'under the original resident public EXP and one-subtraction premises. Nonempty cells identify '
        'neither a logarithm nor the private exponential. Original PMF and public LN readings are '
        'unchanged; no new Excel question, fitted constant or product change.')
    write(output,result);print(json.dumps(counts));return result


def poisson_one_primitive_contrasts(original,output,library):
    """Separate signed EXP reduction from the PMF1 log/subtract producer."""
    path=Path(original)
    prior=json.loads(gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes())
    if prior.get('complete') is not True:raise ValueError('incomplete original PMF1 observations')
    graphs={}
    # Source forms fixed in this order, before consulting any target matches.
    for log in ('public_LN','resident_rom_LN','resident_binary64_ln2_LN'):
        for stored_argument in (True,False):
            for reciprocal in (False,True):
                for reduction in ('rom_log2e','binary64_log2e'):
                    g=Builder(2)
                    if log=='public_LN':leaf=1
                    else:
                        coefficient=(g.op('fldln2',store=False) if log=='resident_rom_LN'
                                     else g.constant(bits(math.log(2.))))
                        leaf=g.op('fyl2x',coefficient,0,store=False)
                    argument=g.op('sub',leaf,0,store=False)
                    if stored_argument:g.nodes[argument]['site']=0;g.sites=1
                    if reciprocal:argument=g.op('neg',argument,store=False)
                    coefficient=(g.op('fldl2e',store=False) if reduction=='rom_log2e'
                                 else g.constant(bits(math.log2(math.e))))
                    t=g.op('mul',argument,coefficient,store=False)
                    whole=g.op('frndint',t,store=False)
                    small=g.op('f2xm1',g.op('sub',t,whole,store=False),store=False)
                    result=g.op('fscale',g.op('add',g.one,small,store=False),whole,store=False)
                    if reciprocal:result=g.op('div',g.one,result,store=False)
                    key=f'{log}/{"stored" if stored_argument else "resident"}/{"reciprocal" if reciprocal else "signed"}/{reduction}'
                    graphs[key]=g.graph(result,result)
    native=Native(library);counts={k:Counter() for k in graphs};records=[]
    for row in prior['records']:
        if row['public_LN']['type']!='number' or row['excel']['type']!='number':raise ValueError('nonnumeric original primitive control')
        inputs=[int(row['input'],16),int(row['public_LN']['bits'],16)];predictions={}
        for name,graph in graphs.items():
            prediction=f'{native.publish(native.evaluate(graph,inputs,0x133f,mask=1)[:10]):016x}'
            counts[name]['matched' if prediction==row['excel']['bits'] else 'mismatched']+=1
            predictions[name]=prediction
        records.append(dict(input=row['input'],public_LN=row['public_LN'],excel=row['excel'],predictions=predictions))
    result=dict(format='dnacalc-poisson-one-primitive-contrasts/1',complete=True,
        source_bindings={str(path.resolve()):sha(path),str(Path(library).resolve()):sha(library),
                         str(ROOT/'racer/x87_native.c'):sha(ROOT/'racer/x87_native.c'),
                         str(ROOT/'racer/x87_programs.py'):sha(ROOT/'racer/x87_programs.py')},
        original_receipts=prior['original_receipts'],graphs=graphs,records=records,counts=counts,
        limits='Count1 only: EXP(L(mean)-mean), without factorial or normalizer. Three fixed '
        'log leaves, two argument publications, signed or reciprocal exponential, and two '
        'fixed reduction constants. A single contradiction refutes its complete declared '
        'composition. No least-bad programme, fitted coefficient, cutoff or private helper '
        'is selected. Original PMF1/LN readings and environments remain unchanged.')
    write(output,result);print(json.dumps(dict(programs=len(graphs),complete_survivors=[k for k,v in counts.items() if not v['mismatched']])));return result


def poisson_exp_argument_study(checkout,original,output,cache,library,offline=False):
    """Project a single PMF exponential input onto the binary64 lattice.

    The original mass is the target. No log, factorial or mass program is
    selected. New public EXP calls bracket targets whose primitive preimage
    contains no binary64 argument.
    """
    from x87 import decode
    prior=json.loads(Path(original).read_text())
    if not prior.get('complete'):raise ValueError('incomplete Poisson stage observations')
    native=Native(library);g=Builder(1)
    t=g.op('mul',0,g.op('fldl2e',store=False),store=False);whole=g.op('frndint',t,store=False)
    small=g.op('f2xm1',g.op('sub',t,whole,store=False),store=False)
    positive=g.op('fscale',g.op('add',g.one,small,store=False),whole,store=False)
    out=g.op('div',g.one,positive,store=False);graph=g.graph(out,out)
    rows=[];questions=set();counts=Counter()
    for i,r in enumerate(prior['records']):
        answer=r['published_output'];row=dict(index=i,input=r['input'],stage=r['stage'],excel=answer)
        if answer['type']!='number' or value(answer['bits'])<=0.:
            row['status']='no positive numeric target';counts[row['status']]+=1;rows.append(row);continue
        cell=positive_extended_cell(native,graph,answer['bits'],Q(-math.log(value(answer['bits']))),decreasing=True)
        if cell is None:raise ValueError('empty extended exponential preimage needs a separate certificate')
        lo,hi=(decode(bytes.fromhex(cell[k]))for k in ('low','high'))
        a=BINARY64.ceiling_coordinate(lo);b=BINARY64.floor_coordinate(hi)
        count=max(0,b-a+1)
        row.update(cell=cell,binary64_argument_count=count,
            first_binary64_magnitudes=[bits(float(BINARY64.value(k)))for k in range(a,min(b+1,a+8))])
        row['status']='binary64 argument possible'if count else'no binary64 argument'
        if not count:
            lower=BINARY64.floor_coordinate(lo);upper=BINARY64.ceiling_coordinate(hi)
            if upper-lower!=1:raise ValueError('missing argument neighbors are not adjacent')
            row['magnitude_neighbors']=[bits(float(BINARY64.value(k)))for k in (lower,upper)]
            questions.update(row['magnitude_neighbors'])
        counts[row['status']]+=1;rows.append(row)
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    magnitudes=sorted(questions)
    actual=dict(zip(magnitudes,service.ask('EXP',[[client.number(-value(x))]for x in magnitudes])['rows']))if magnitudes else{}
    for row in rows:
        if'magnitude_neighbors'not in row:continue
        lo,hi=[actual[x]for x in row['magnitude_neighbors']]
        bracket=lo['type']==hi['type']=='number'and int(hi['bits'],16)<int(row['excel']['bits'],16)<int(lo['bits'],16)
        row.update(public_EXP_neighbors=[lo,hi],public_bracket_verified=bracket)
        counts['public bracket verified'if bracket else'public bracket contradicted']+=1
    result=dict(format='dnacalc-poisson-exponential-input-preimages/1',complete=True,
        original_sha256=sha(original),original_receipts=prior['receipts'],graph=graph,records=rows,counts=counts,receipts=service.receipts,
        native_bindings={str(p):sha(p)for p in (Path(library),ROOT/'racer/x87_native.c',ROOT/'racer/x87_programs.py')},
        limits='One conditional public EXP primitive, negative as reciprocal of positive exp. Original PMF outputs are inverted without assuming their exponent producer. A verified adjacent public bracket excludes a stored binary64 exponential argument within this primitive; it does not identify a resident argument, private exponential or full PMF algorithm. No coefficients or count-selected correction adopted.')
    write(output,result);print(json.dumps(counts));return result


def root_store_contrasts(original,output,library):
    """Remove one publication at a time from the paired primitive composition."""
    prior=json.loads(Path(original).read_text())
    if not prior.get('complete'):raise ValueError('incomplete paired primitives')
    native=Native(library);graphs={};counts={k:Counter()for k in ('stored','resident_log','resident_reciprocal','resident_argument')}
    masks={'stored':7,'resident_log':6,'resident_reciprocal':5,'resident_argument':3}
    for n in sorted({r['count']for r in prior['records']}):
        for negative in (False,True):
            g=Builder(1);log=g.op('fyl2x',g.op('fldln2',store=False),0)
            reciprocal=g.op('div',g.one,g.constant(bits(n)))
            t=g.op('mul',log,reciprocal)
            if negative:t=g.op('neg',t,store=False)
            t=g.op('mul',t,g.op('fldl2e',store=False),store=False)
            whole=g.op('frndint',t,store=False);small=g.op('f2xm1',g.op('sub',t,whole,store=False),store=False)
            out=g.op('fscale',g.op('add',g.one,small,store=False),whole,store=False)
            if negative:out=g.op('div',g.one,out,store=False)
            graphs[f'{n}/{negative}']=g.graph(out,out)
    records=[]
    for r in prior['records']:
        graph=graphs[f'{r["count"]}/{value(r["input"])<1.}']
        predictions={name:f'{native.publish(native.evaluate(graph,[int(r["input"],16)],0x133f,mask=mask)[:10]):016x}'for name,mask in masks.items()}
        for name,prediction in predictions.items():counts[name]['matched'if r['excel']==dict(type='number',bits=prediction)else'mismatched']+=1
        records.append(dict(count=r['count'],input=r['input'],excel=r['excel'],predictions=predictions))
    result=dict(format='dnacalc-root-single-store-contrasts/1',complete=True,original_sha256=sha(original),
        masks=masks,graphs=graphs,records=records,counts=counts,
        native_bindings={str(p):sha(p)for p in (Path(library),ROOT/'racer/x87_native.c',ROOT/'racer/x87_programs.py')},
        limits='One publication removed at a time, with all other paired stages fixed. Native FYL2X and the actual public EXP core (negative as the resident reciprocal of positive exp) are declared conditional primitives; no Cartesian store-mask or whole-program search. Original paired-study scope and provenance are unchanged. These contrasts distinguish sampled stages within this model, not all possible private algorithms.')
    write(output,result);print(json.dumps(counts));return result


def root_public_stages(checkout,original,output,cache,offline=False):
    """Pair each measured reduced root with independent public leaf calls.

    The quotient/publication shape is fixed by the exact-power stage study,
    before these answers. This neither searches whole programs nor substitutes
    a calculated value for an Excel observation.
    """
    prior=json.loads(Path(original).read_text())
    if not prior.get('complete'):raise ValueError('incomplete normalization study')
    roots={}
    for row in prior['records']:
        item=next(c for c in row['conventions']if c['convention']=='adjusted_frexp_truncation')
        key=row['count'],item['reduced'];answer=item['measured_normalized_root']
        if key in roots and roots[key]!=answer:raise ValueError('inconsistent original root control')
        roots[key]=answer
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    inputs=sorted({x for _,x in roots})
    logs=dict(zip(inputs,service.ask('LN',[[client.number(x)]for x in inputs])['rows']))
    singletons=dict(zip(inputs,service.ask('GEOMEAN',[[client.number(x)]for x in inputs])['rows']))
    arguments={};reciprocals={}
    for n,x in roots:
        a=logs[x]
        if a['type']!='number':raise ValueError('reduced log must be numeric')
        # Both assignments publish from PC64. The independently observed leaf
        # is binary64; it is not an inferred resident logarithm.
        reciprocal=BINARY64.round(EXTENDED64.round(Q(1,n)))
        argument=BINARY64.round(EXTENDED64.round(Q(value(a['bits']))*reciprocal))
        reciprocals[n]=bits(float(reciprocal));arguments[n,x]=bits(float(argument))
    exp_inputs=sorted(set(arguments.values()))
    exps=dict(zip(exp_inputs,service.ask('EXP',[[client.number(x)]for x in exp_inputs])['rows']))
    records=[];counts=Counter()
    for(n,x),answer in sorted(roots.items()):
        argument=arguments[n,x];prediction=exps[argument]
        counts['matched'if answer==prediction else'mismatched']+=1
        records.append(dict(count=n,input=x,excel=answer,public_LN=logs[x],reciprocal=reciprocals[n],
            stored_argument=argument,public_EXP=prediction,matched=answer==prediction))
    singleton_counts=Counter('identity'if a==dict(type='number',bits=x)else'different'for x,a in singletons.items())
    result=dict(format='dnacalc-root-paired-public-stages/1',complete=True,original_sha256=sha(original),
        original_receipts=prior['receipts'],records=records,counts=counts,
        singleton_controls=[dict(input=x,excel=a)for x,a in singletons.items()],singleton_counts=singleton_counts,
        public_log_calls=len(inputs),public_exp_calls=len(exp_inputs),receipts=service.receipts,
        limits='Adjusted normal reduced products only, with original unit-factor GEOMEAN controls. Every public LN and EXP is independently measured at its actual binary64 input. Assignments use PC64 then binary64 publication. Agreement establishes the declared sampled composition, not an arbitrary resident helper, table/polynomial internals or overflow fallback. One-element controls separately test identity.')
    write(output,result);print(json.dumps(dict(counts=counts,singletons=singleton_counts,public_log_calls=len(inputs),public_exp_calls=len(exp_inputs))))
    return result


def root_log_stage(checkout,original,output,cache,library,offline=False):
    """Invert one quotient stage, leaving each private logarithm arbitrary.

    Positive fractional power anchors have q=0 in the measured normalization.
    Group by the same exact product. Contrast one quotient publication and
    division versus a stored reciprocal under a declared exponential primitive.
    """
    from exact import Span
    from x87 import decode
    prior=json.loads(Path(original).read_text())
    if not prior.get('complete'):raise ValueError('incomplete original anchors')
    rows=[r for r in prior['records']if r['role']=='fractional_log2_root'and r['product_power']>0]
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=180),cache,offline)
    powers=sorted({r['product_power']for r in rows})
    public=dict(zip(powers,service.ask('LN',[[client.number(math.ldexp(1.,k))]for k in powers])['rows']))
    native=Native(library);ln2=decode(native.leaf('fldln2',bytes(10),bytes(10),0x133f))
    records=[];models=[]
    for quotient in ('divide','stored_reciprocal'):
        for publish in (False,True):
            groups={};used=[]
            for row in rows:
                n=row['count'];g=Builder(1)
                operand=n if quotient=='divide'else 1./n
                t=g.op('div'if quotient=='divide'else'mul',0,g.constant(bits(operand)),store=publish)
                t=g.op('mul',t,g.op('fldl2e',store=False),store=False)
                whole=g.op('frndint',t,store=False);fraction=g.op('sub',t,whole,store=False)
                small=g.op('f2xm1',fraction,store=False)
                out=g.op('fscale',g.op('add',g.one,small,store=False),whole,store=False)
                graph=g.graph(out,out)
                class Bound:
                    def evaluate(self,graph,inputs,cw,extended):
                        return native.evaluate(graph,inputs,cw,mask=int(publish),extended=extended)
                    def publish(self,raw):return native.publish(raw)
                cell=positive_extended_cell(Bound(),graph,row['excel']['bits'],ln2*row['product_power'])
                domain=Span(decode(bytes.fromhex(cell['low'])),decode(bytes.fromhex(cell['high'])))if cell else None
                k=row['product_power'];old=groups.get(k)
                groups[k]=domain if k not in groups else old.meet(domain)if old is not None and domain is not None else None
                used.append(k)
                records.append(dict(quotient=quotient,published=publish,count=n,product_power=k,
                    input=row['input'],excel=row['excel'],cell=cell,graph=graph))
            summaries=[]
            for k,domain in sorted(groups.items()):
                answer=public[k]
                summaries.append(dict(product_power=k,observations=used.count(k),
                    common_private_log=domain.json()if domain is not None else None,public_LN=answer,
                    public_LN_possible=domain is not None and answer['type']=='number'and domain.contains(Q(value(answer['bits'])))))
            models.append(dict(quotient=quotient,published=publish,groups=summaries))
    result=dict(format='dnacalc-private-log-quotient-stage/1',complete=True,
        original_sha256=sha(original),original_receipts=prior['receipts'],records=records,models=models,
        receipts=service.receipts,native_bindings={str(p):sha(p)for p in (Path(library),ROOT/'racer/x87_native.c',ROOT/'racer/x87_programs.py')},
        limits='Complete positive fractional original anchor subset, q=0. Each same-input log is arbitrary and shared only across counts. The native public EXP model and monotone preimages are conditional premises. An impossible interval rejects this quotient/primitive factorization, not all logarithms or exponentials; possible intervals do not identify a private implementation.')
    write(output,result)
    print(json.dumps([dict(quotient=m['quotient'],published=m['published'],groups=len(m['groups']),
        impossible=sum(g['common_private_log']is None for g in m['groups']),
        excludes_public=sum(not g['public_LN_possible']for g in m['groups']))for m in models]))
    return result


def root_primitive_spill_study(checkout,output,cache,library,count=256):
    """Perturb one primitive store at a time; leave the log producer unknown.

    This inverts the same original roots, with a store either at the reduced
    fraction or at F2XM1's return. There is no Cartesian store-mask search and
    no candidate logarithm is fitted to the intervals.
    """
    from x87 import decode
    native=Native(library);client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,True)
    xs=[f'{b:016x}'for x in (1.3,7.,1000.)for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    roots={n:service.ask('GEOMEAN',[[client.number(x)]+[client.number(1.)]*(n-1)for x in xs])['rows']for n in (2,4,8)}
    log_graphs={k:primitive_stage_graph('LN'if k=='EXP'else'LOG2')for k in ('EXP','EXP2')}
    class Stored:
        def evaluate(self,graph,inputs,cw,extended):return native.evaluate(graph,inputs,cw,mask=1,extended=extended)
        def publish(self,raw):return native.publish(raw)
    stored=Stored();records=[];models=[]
    for kind in ('EXP','EXP2'):
        for site,operation in [('fraction','sub'),('f2xm1_return','f2xm1')]:
            graph=primitive_stage_graph(kind)
            nodes=[n for n in graph['nodes']if n['op']==operation]
            if len(nodes)!=1:raise ValueError('primitive store site is not unique')
            nodes[0]['site']=0;graph['sites']=1
            tally=Counter()
            for index,x in enumerate(xs):
                center=decode(native.evaluate(log_graphs[kind],[int(x,16)],0x133f)[:10])
                lows=[];highs=[];cells=[]
                for n in (2,4,8):
                    a=roots[n][index]
                    cell=positive_extended_cell(stored,graph,a['bits'],center/n)
                    cells.append(dict(count=n,excel=a,cell=cell))
                    if cell:
                        lows.append(decode(bytes.fromhex(cell['low']))*n)
                        highs.append(decode(bytes.fromhex(cell['high']))*n)
                possible=len(lows)==3 and max(lows)<=min(highs)
                tally['shared_log_possible'if possible else'shared_log_impossible']+=1
                if possible:
                    tally['registered_log_possible'if max(lows)<=center<=min(highs)else'registered_log_excluded']+=1
                records.append(dict(input=x,primitive=kind,store=site,cells=cells,shared_log_possible=possible,
                                    common_log_interval=[str(max(lows)),str(min(highs))]if possible else None))
            models.append(dict(primitive=kind,store=site,graph=graph,counts=tally))
    result=dict(format='dnacalc-private-root-single-primitive-stores/1',complete=True,
        models=models,records=records,receipts=service.receipts,
        limits='Only one binary64 spill is changed inside each fixed EXP/EXP2 primitive. All three arities share an otherwise arbitrary resident logarithm; n2/4/8 scaling is exact. A feasible intersection is a necessary condition, never an identified logarithm. Local monotonicity is conditional. Refutations apply to these single-stage perturbations only, not to an entire root or polynomial/table family. Original readings reused; no new Excel calls or product change.')
    write(output,result);print('single primitive stores',[(m['primitive'],m['store'],dict(m['counts']))for m in models],flush=True)
    return result


def bessel_admission_study(checkout,output,cache,offline=False):
    """Separate logical Value2 cells, blanks and error order from omitted args.

    Existing workbook captures settle BESSELJ's literal logical and omitted
    arguments. This asks only the distinct cell transport and error controls.
    """
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    rows=[];roles=[]
    for logical in (False,True):
        for position in (0,1):
            row=[client.number(1.1),client.number(3.)];row[position]=client.boolean(logical)
            rows.append(row);roles.append('logical_cell_position_'+str(position))
    for position in (0,1):
        row=[client.number(1.1),client.number(3.)];row[position]=client.blank()
        rows.append(row);roles.append('blank_cell_position_'+str(position))
    for row,role in [([client.error('#DIV/0!'),client.boolean(True)],'first_error_before_second_logical'),
                     ([client.boolean(True),client.error('#DIV/0!')],'first_logical_before_second_error'),
                     ([client.boolean(True),client.number(-1.)],'logical_before_numeric_domain')]:
        rows.append(row);roles.append(role)
    write(str(output)+'.questions.json',dict(format='dnacalc-bessel-admission-questions/1',
        complete=True,rows=rows,roles=roles,selection='Fixed logical-cell, blank-cell and argument-order controls; no existing literal/omitted workbook question repeated.'))
    answers=service.ask('BESSELJ',rows)['rows']
    result=dict(format='dnacalc-bessel-admission/1',complete=True,
        records=[dict(role=role,args=row,answer=answer) for role,row,answer in zip(roles,rows,answers)],
        receipts=service.receipts,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Original scalar Value2 controls only. Missing and blank are distinct; omitted '
        'arguments remain bound to their existing workbook captures. No numeric Bessel '
        'kernel, other Bessel-family entry or reference-coercion algorithm is inferred.')
    write(output,result);print(json.dumps([(r['role'],r['answer']) for r in result['records']]));return result


def hypergeom_return_study(checkout,output,cache,offline=False):
    """Cancel identical mass factors to constrain one return/store boundary.

    At zero/all population successes, the mathematical mass is exactly one.
    Dyadic sampling probabilities also remove complement-rounding uncertainty.
    A public BINOM reading supplies a conditional stored-mass control; it is
    never substituted as an answer for the private HYPGEOM call.
    """
    from exact import Span,pull_rounding
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    keys=[]
    for total in (2,4,8,16,32,64,128,256,1024):
        for draws in sorted({1,2,3,total//2,total-1}):
            if 0<draws<total:keys.append((draws,total))
    binoms=service.ask('BINOM.DIST',[[client.number(d),client.number(n),client.number(d/n),client.boolean(False)]for d,n in keys])['rows']
    queries=[]
    for draws,total in keys:
        for population_successes in (0,total):
            successes=0 if population_successes==0 else draws
            for cumulative in (False,True):
                queries.append([client.number(successes),client.number(draws),client.number(population_successes),client.number(total),client.boolean(cumulative)])
    answers=service.ask('HYPGEOM.DIST',queries)['rows']
    legacy=service.ask('HYPGEOMDIST',[r[:4]for r in queries[::2]])['rows']
    records=[];counts=Counter();alias=Counter();one=client.number(1.)
    for index,((draws,total),b)in enumerate(zip(keys,binoms)):
        group=answers[4*index:4*index+4]
        counts['zero_all_same'if group[0]==group[2]else'zero_all_differ']+=1
        for side in (0,1):
            a,cdf=group[2*side:2*side+2]
            alias['equal'if legacy[2*index+side]==a else'different']+=1
            counts['pmf_mathematical_one'if a==one else'pmf_differs_from_one']+=1
            counts['cdf_mathematical_one'if cdf==one else'cdf_differs_from_one']+=1
            domains={}
            if a['type']=='number'and b['type']=='number'and value(b['bits'])>0:
                stored=Q(value(b['bits']));got=Q(value(a['bits']))
                stored_cell=BINARY64.pullback(Span(stored,stored))
                quotient_cell=pull_rounding(Span(got,got),'stored_pc64')
                for orientation in ('resident_over_stored','stored_over_resident'):
                    domain=(quotient_cell.linear(stored)if orientation=='resident_over_stored'else quotient_cell.reciprocal(stored))
                    domain=stored_cell.meet(domain)
                    domains[orientation]=domain.json()if domain is not None else None
                    counts[orientation+'_possible'if domain is not None else orientation+'_impossible']+=1
            records.append(dict(draws=draws,total=total,population_successes=0 if side==0 else total,
                public_binom=b,pmf=a,cdf=cdf,legacy=legacy[2*index+side],conditional_resident_mass_intervals=domains))
    result=dict(format='dnacalc-hypergeom-identical-factor-return-stage/1',complete=True,
        records=records,counts=counts,alias_counts=alias,receipts=service.receipts,
        limits='Dyadic n/N, zero/all population successes, exact identity factors. The return-stage intervals assume the private identical mass has the measured public BINOM binary64 publication and one RN64 then RN53 division. A contradiction rejects that return-stage premise, not a whole mass implementation. Mathematical one controls do not replace Excel outputs. No fitted correction or change is adopted.')
    write(output,result);print('hypergeom return stage',dict(counts),dict(alias),flush=True)
    return result


def multiplication_seed_interval(seed,output,multiplier):
    """Unknown resident seed publishing to seed, consumed by one stored multiply."""
    from exact import Span,pull_rounding
    if multiplier<=0:raise ValueError('positive multiplier required')
    source=BINARY64.pullback(Span(seed,seed))
    target=pull_rounding(Span(output,output),'stored_pc64').linear(1/multiplier)
    return source.meet(target)


def poisson_seed_study(checkout,output,cache,offline=False,count=256):
    """Infer one seed-consumer boundary without specifying the seed producer.

    P0, P1 and P2 separately constrain the unknown resident seed at the
    next multiply. Mode pairs at integer means have multiplier exactly one,
    removing multiplication precision from that necessary recurrence test.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    xs=[f'{b:016x}'for x in (.13,.7,1.,3.7)for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    rows=[[client.number(float(k)),client.number(x),client.boolean(False)]for k in (0,1,2)for x in xs]
    answers=service.ask('POISSON.DIST',rows)['rows']
    aliases=service.ask('POISSON',rows)['rows']
    exps=service.ask('EXP',[[client.number(bits(-value(x)))]for x in xs])['rows']
    records=[];counts=Counter();n=len(xs)
    for i,x in enumerate(xs):
        p0,p1,p2=[answers[k*n+i]for k in (0,1,2)]
        counts['zero_equals_EXP'if p0==exps[i]else'zero_differs_from_EXP']+=1
        for stage,seed,target,multiplier in [('zero_to_one',p0,p1,Q(value(x))),('one_to_two',p1,p2,Q(value(x))/2)]:
            if any(a['type']!='number'for a in (seed,target)):
                counts[stage+'_nonnumeric']+=1;domain=None
            else:
                domain=multiplication_seed_interval(Q(value(seed['bits'])),Q(value(target['bits'])),multiplier)
                counts[stage+'_possible'if domain is not None else stage+'_impossible']+=1
            records.append(dict(input=x,stage=stage,published_seed=seed,published_output=target,
                                resident_seed_interval=domain.json()if domain is not None else None))
    means=(*range(1,65),128,256,500,1000)
    mode_rows=[[client.number(float(k)),client.number(float(mean)),client.boolean(False)]for mean in means for k in (mean-1,mean)]
    modes=service.ask('POISSON.DIST',mode_rows)['rows'];mode_records=[]
    for i,mean in enumerate(means):
        a,b=modes[2*i:2*i+2]
        counts['exact_mode_step_equal'if a==b else'exact_mode_step_different']+=1
        mode_records.append(dict(mean=mean,previous_count=mean-1,previous=a,current=b))
    result=dict(format='dnacalc-poisson-resident-seed-multiply-stage/1',complete=True,
        records=records,mode_records=mode_records,counts=counts,
        alias_counts=Counter('equal'if a==b else'different'for a,b in zip(answers,aliases)),receipts=service.receipts,
        limits='One positive stored-PC64 multiplication at a time, with any resident seed whose publication is the independently measured prior mass. Empty intervals refute that factorization, not every PMF algorithm. Integer-mean mode pairs have exact multiplier one. Public EXP controls establish only the zero-count route on these windows. No recurrence, correction, threshold or leaf is adopted from counts; all outputs remain original Excel observations.')
    write(output,result);print('poisson seed stage',dict(counts),flush=True)
    return result


def root_equivalence_study(checkout,output,cache,offline=False):
    """Equal rational log2 arguments expose a shared primitive independently.

    Product2^k and arityn have exact log2 ratio k/n. Equivalent fractions
    must reach one common exponential under an exact-log2/direct-division
    factorization, regardless of the exponential's numerical implementation.
    """
    from exact import Span,pull_rounding
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    fractions=(Q(1,7),Q(1,5),Q(1,4),Q(1,3),Q(1,2),Q(2,3),Q(1),Q(3,2),Q(5,2))
    keys=[]
    for q in fractions:
        for sign in (1,-1):
            for m in range(1,32//q.denominator+1):
                n=q.denominator*m;k=sign*q.numerator*m
                if n>=2:keys.append((q*sign,n,k))
    questions=[[client.number(math.ldexp(1.,k))]+[client.number(1.)]*(n-1)for q,n,k in keys]
    answers=service.ask('GEOMEAN',questions)['rows']
    powers=service.ask('POWER',[[row[0],client.number(1./n)]for (q,n,k),row in zip(keys,questions)])['rows']
    records=[];groups={};by_key={}
    for (q,n,k),a,p in zip(keys,answers,powers):
        if a['type']!='number':raise ValueError('nonnumeric exact-ratio root')
        r=dict(ratio=str(q),count=n,product_power=k,excel=a,public_power=p)
        records.append(r);groups.setdefault(str(q),[]).append(r);by_key[q,n]=a
    group_records=[]
    for ratio,rows in groups.items():
        outputs={r['excel']['bits']for r in rows};public={r['public_power']['bits']for r in rows}
        group_records.append(dict(ratio=ratio,observations=len(rows),private_outputs=sorted(outputs),
                                  public_outputs=sorted(public),consistent_shared_exact_ratio=len(outputs)==1))
    reciprocals=[];counts=Counter()
    for (q,n),a in by_key.items():
        if q<=0:continue
        b=by_key[-q,n];positive=Q(value(a['bits']));negative=Q(value(b['bits']))
        source=BINARY64.pullback(Span(positive,positive))
        inverse=pull_rounding(Span(negative,negative),'stored_pc64').reciprocal(Q(1))
        intersection=source.meet(inverse)
        counts['resident_reciprocal_possible'if intersection is not None else'resident_reciprocal_impossible']+=1
        reciprocals.append(dict(ratio=str(q),count=n,positive=a,negative=b,
                                 resident_positive_interval=intersection.json()if intersection is not None else None))
    result=dict(format='dnacalc-private-root-exact-ratio-equivalence/1',complete=True,
        records=records,groups=group_records,reciprocals=reciprocals,counts=counts,receipts=service.receipts,
        limits='Exact powers of two and unit-factor products. Equal-ratio agreement is a necessary condition, not a logarithm or exponential implementation identification. Disagreement refutes a common exponential fed by exact log2 and direct count division, even if that exponential is unknown. Reciprocal intervals allow any resident positive output with the observed publication; they test one final reciprocal only. No table, coefficient, threshold or programme is selected by counts.')
    write(output,result);print('root equivalence',[(g['ratio'],len(g['private_outputs']),len(g['public_outputs']))for g in group_records],dict(counts),flush=True)
    return result


def root_scale_study(checkout,output,cache,offline=False,count=128):
    """Move only the exact binary exponent by multiples of the root count."""
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    xs=[f'{b:016x}'for x in (1.3,7.,1000.)for b in range(int(bits(x),16)-count//2,int(bits(x),16)+count//2)]
    scales=(0,-20,-1,1,20);records=[];counts={}
    for n in (2,3,4,5,7,9,10,11,16,32):
        rows=[[client.number(math.ldexp(value(x),n*s))]+[client.number(1.)]*(n-1)for s in scales for x in xs]
        answers=service.ask('GEOMEAN',rows)['rows'];base=answers[:len(xs)];tally=Counter()
        for j,s in enumerate(scales[1:],1):
            for i,x in enumerate(xs):
                a,b=base[i],answers[j*len(xs)+i]
                if a['type']!='number'or b['type']!='number':raise ValueError('nonnumeric scaled root')
                expected=bits(math.ldexp(value(a['bits']),s))
                equal=b['bits']==expected;tally['equal'if equal else'different']+=1
                records.append(dict(count=n,output_scale=s,input=x,scaled_input=bits(math.ldexp(value(x),n*s)),
                    base=a,scaled=b,exact_rescaling_of_observed_base=expected))
        counts[str(n)]=tally
    result=dict(format='dnacalc-root-exact-exponent-transfer/1',complete=True,
        consecutive_inputs_per_center=count,records=records,counts=counts,receipts=service.receipts,
        limits='Exact power-of-two changes to the unit-factor product, with scale a multiple of arity. The compared base and scaled outputs are independent original observations. Agreement supports only the measured transfer; it does not identify an exponent-reduction algorithm or permit an unmeasured dispatch. No numerical implementation changes.')
    write(output,result);print('root exponent transfer',{k:dict(v)for k,v in counts.items()},flush=True)
    return result


def root_output_edges(checkout,output,cache,offline=False):
    """Bisect published output transitions without assuming a private primitive.

    The rational midpoint-power supplies the initial neighborhood only. The
    observed crossing and consecutive window are retained even if the local
    function is nonmonotone or has multiple crossings.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    brackets=[]
    for n in (2,3,4,5,7,9,10,11,16,32):
        for y in (1.,1.125,1.25,1.375,1.5,1.625,1.75,1.875,2.):
            target=int(bits(y),16);previous=value(f'{target-1:016x}')
            midpoint=(Q(y)+Q(previous))/2
            ideal_input=midpoint**n;center=int(bits(float(ideal_input)),16)
            brackets.append(dict(count=n,target=bits(y),midpoint=str(midpoint),
                                 ideal_input=str(ideal_input),low=f'{center-1024:016x}',high=f'{center+1024:016x}'))
    def row(s,x):return [client.number(x)]+[client.number(1.)]*(s['count']-1)
    def side(s,x,a):return int(a['bits'],16)>=int(s['target'],16)if a['type']=='number'else None
    edges=bisect_observed_edges(service,'GEOMEAN',brackets,row,side)
    requests=[(e,f'{b:016x}')for e in edges for b in range(int(e['high'],16)-64,int(e['high'],16)+65)]
    windows=[]
    for start in range(0,len(requests),4096):
        part=requests[start:start+4096]
        windows.extend(service.ask('GEOMEAN',[row(e,x)for e,x in part])['rows'])
    records=[];counts=Counter()
    for i,e in enumerate(edges):
        xs=[x for _,x in requests[i*129:(i+1)*129]];answers=windows[i*129:(i+1)*129]
        if any(a['type']!='number'for a in answers):raise ValueError('nonnumeric root edge window')
        transitions=[];decreases=[]
        for j in range(1,len(xs)):
            a,b=answers[j-1:j+1]
            if side(e,None,a)!=side(e,None,b):transitions.append(dict(low=xs[j-1],high=xs[j],low_answer=a,high_answer=b))
            if int(a['bits'],16)>int(b['bits'],16):decreases.append(dict(low=xs[j-1],high=xs[j],low_answer=a,high_answer=b))
        lo,hi=(Q(value(e[k]))for k in ('low','high'))
        exact=Q(e['ideal_input'])
        unit=Q(math.ulp(value(e['high'])))
        offset=[(lo-exact)/unit,(hi-exact)/unit]
        counts['unique_local_crossing'if len(transitions)==1 else'multiple_local_crossings']+=1
        counts['monotone_local_window'if not decreases else'nonmonotone_local_window']+=1
        records.append(dict(e,window=dict(inputs=xs,answers=answers,transitions=transitions,decreases=decreases),
                            reference_input_ulp=str(unit),input_offset_from_exact_midpoint_power_in_input_ulps=list(map(str,offset))))
    result=dict(format='dnacalc-root-observed-output-edges/1',complete=True,
        records=records,counts=counts,receipts=service.receipts,
        limits='One directly observed adjacent false/true crossing per bracket. The full local window records every reversal and threshold crossing; no uniqueness or global monotonicity is assumed. Exact mathematical midpoint powers locate questions but never supply Excel answers. Edge offsets are observations, not fitted coefficients or an identified private algorithm.')
    write(output,result);print('root output edges',dict(counts),flush=True)
    return result


def root_normalization_study(checkout,original,output,cache,offline=False):
    """Isolate binary-exponent reduction using independently measured core roots.

    Test the three fixed exponent conventions without specifying the private
    root body. Each normalized product gets its own original GEOMEAN reading.
    The final power-of-two reconstruction is exact on this normal cohort.
    """
    original=Path(original);doc=json.loads(original.read_text())
    if not doc['complete']or doc['format']!='dnacalc-root-exact-exponent-transfer/1':raise ValueError('complete exponent-transfer study required')
    client,_=load_service(checkout)
    prior=CachedService(client.Client(timeout=180),original.parent/'root-scales-service',True)
    rows=[]
    for receipt in doc['receipts']:
        raw_path=prior.directory/receipt['file']
        if sha(raw_path)!=receipt['sha256']:raise ValueError('original exponent-transfer reading changed')
        with gzip.open(raw_path,'rt')as stream:raw=json.load(stream)
        accept(raw['reply'],'GEOMEAN',raw['request']['rows'])
        rows.extend(zip(raw['request']['rows'],raw['reply']['rows']))
    conventions=('frexp_truncation','adjusted_frexp_truncation','adjusted_frexp_floor')
    controls={};descriptions=[]
    for i,(args,a)in enumerate(rows):
        n=len(args);p=args[0]['bits'];b=int(p,16);e=((b>>52)&0x7ff)-1023
        if not 0<(b>>52)&0x7ff<0x7ff or b>>63:raise ValueError('normal positive product required')
        items=[]
        for convention in conventions:
            exponent=e+int(convention=='frexp_truncation')
            q=exponent//n if convention.endswith('_floor')else(abs(exponent)//n)*(-1 if exponent<0 else 1)
            reduced=bits(math.ldexp(value(p),-q*n));key=n,reduced
            controls[key]=None
            items.append(dict(convention=convention,exponent=exponent,quotient=q,reduced=reduced))
        descriptions.append(dict(index=i,count=n,input=p,excel=a,conventions=items))
    service=CachedService(client.Client(timeout=180),cache,offline)
    keys=sorted(controls)
    for start in range(0,len(keys),4096):
        part=keys[start:start+4096]
        answers=service.ask('GEOMEAN',[[client.number(x)]+[client.number(1.)]*(n-1)for n,x in part])['rows']
        for key,a in zip(part,answers):controls[key]=a
    counts={name:Counter()for name in conventions}
    for r in descriptions:
        for item in r['conventions']:
            a=controls[r['count'],item['reduced']]
            if a['type']!='number'or r['excel']['type']!='number':raise ValueError('nonnumeric normalized root')
            predicted=dict(type='number',bits=bits(math.ldexp(value(a['bits']),item['quotient'])))
            item.update(measured_normalized_root=a,prediction=predicted)
            counts[item['convention']]['matched'if predicted==r['excel']else'mismatched']+=1
    result=dict(format='dnacalc-private-root-exponent-normalization/1',complete=True,
        original_study_sha256=sha(original),original_receipts=doc['receipts'],
        convention_order=conventions,counts=counts,records=descriptions,receipts=service.receipts,
        limits='Exponent reduction only, with the private normalized root independently measured for each argument. Fixed source-shape order: frexp integer truncation, frexp adjusted to[1,2) then integer truncation, adjusted exponent then floor. No inner root formula or parameter is fitted. Agreement identifies the declared stage only on the normal unit-factor cohort; product accumulation, fallback and private root remain separate.')
    write(output,result);print('root normalization',{k:dict(v)for k,v in counts.items()},'controls',len(controls),flush=True)
    return result


def integer_power_study(checkout,output,cache,offline=False,count=2048):
    """Choose multiplication/publication discriminants before observing masses.

    Small integer multipliers fit exactly in PC64, so their windows cannot
    distinguish strict multiplication from the stored step. Large odd trial
    counts at near-one probabilities expose that step without underflow.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=180),cache,offline)
    def ask(name,rows):
        result=[]
        for start in range(0,len(rows),8192):
            result.extend(service.ask(name,rows[start:start+8192])['rows'])
        return result
    seeds=[]
    for k in (8,12,16,20,24,28,32,36,40):
        center=int(bits(1.-2.**-k),16)
        seeds.extend((k,f'{p:016x}')for p in range(center-count//2,center+count//2))
    logs=ask('LN',[[client.number(p)]for _,p in seeds])
    selected=[];controls=[]
    for (k,p),log in zip(seeds,logs):
        if log['type']!='number':raise ValueError('non-numeric power logarithm')
        l=Q(value(log['bits']))
        for offset in (-3,-1,1,3):
            n=2**(k+4)+offset;t=n*l
            a,b=BINARY64.round(t),BINARY64.round(EXTENDED64.round(t))
            row=dict(probability=p,trials=n,logarithm=log['bits'],
                     strict_exponent=bits(float(a)),stored_exponent=bits(float(b)))
            if a!=b:selected.append(row)
            elif len(controls)<256 and int(p,16)%127==0:controls.append(row)
    all_rows=selected+controls
    exponents=sorted({r[key]for r in all_rows for key in ('strict_exponent','stored_exponent')})
    leaves=dict(zip(exponents,ask('EXP',[[client.number(t)]for t in exponents])))
    answers=ask('BINOM.DIST',[[client.number(float(r['trials'])),client.number(float(r['trials'])),
                             client.number(r['probability']),client.boolean(False)]for r in all_rows])
    counts=Counter();records=[]
    for row,a in zip(all_rows,answers):
        strict,stored=(leaves[row[key]]for key in ('strict_exponent','stored_exponent'))
        label=('excel_non_numeric'if a['type']!='number'else
               'all_equal'if a==strict==stored else
               'candidate_equal_excel_differs'if strict==stored else
               'stored_only'if a==stored else'strict_only'if a==strict else'neither')
        counts[label]+=1
        records.append(dict(row,excel=a,strict=strict,stored=stored,classification=label))
    result=dict(format='dnacalc-integer-power-stores/1',complete=True,
        selection='All exponent disagreements under RN53 versus RN53(RN64), chosen from independently measured LN before requesting BINOM. Counts/offsets and probability windows fixed first; controls chosen by input index only.',
        consecutive_inputs_per_center=count,logarithm_calls=len(logs),
        multiply_discriminants=len(selected),controls=len(controls),counts=counts,records=records,receipts=service.receipts,
        excel_types=Counter(a['type']for a in answers),
        candidate_exact=dict(strict=sum(r['excel']==r['strict']for r in records),
                             stored=sum(r['excel']==r['stored']for r in records)),
        limits='All-success masses only, near-one finite normal inputs, large exact odd integers. EXP and LN are measured public binary64 leaves. This identifies only the stated consumer, not an unpublished resident logarithm or every power site. No fitted threshold, coefficient or count-driven selector.')
    write(output,result);print('integer power stores',dict(counts),flush=True)
    return result


def hypergeom_diagonal_study(checkout,output,cache,offline=False):
    """Falling-factorial products with no unknown binomial mass helper."""
    from exact import operation,round_model,Span,pull_rounding,inverse
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=300),cache,offline)
    groups=[]
    for exponent in (4,8,16,24,32,48,52):
        for offset in (0,1):
            total=(1<<exponent)+offset
            for successes in (total//2,total//2+1,total-1):
                groups.append(dict(total=total,successes=successes,draws=list(range(1,min(16,successes)+1))))
    profiles=[dict(form=form,model=model) for form in ('multiply_then_divide','ratio_then_multiply','separate_products')
              for model in ('strict','stored_pc64','resident_pc64')]
    questions=[];rows=[]
    for group in groups:
        for n in group['draws']:
            args=[bits(n),bits(n),bits(group['successes']),bits(group['total'])]
            questions.append([client.number(v) for v in args]+[client.boolean(False)])
            predictions=[]
            for profile in profiles:
                model=profile['model'];v=Q(1);num=Q(1);den=Q(1)
                for i in range(n):
                    a,b=Q(group['successes']-i),Q(group['total']-i)
                    if profile['form']=='multiply_then_divide':
                        v=operation(operation(v,a,'mul',model),b,'div',model)
                    elif profile['form']=='ratio_then_multiply':
                        v=operation(v,operation(a,b,'div',model),'mul',model)
                    else:
                        num=operation(num,a,'mul',model);den=operation(den,b,'mul',model)
                if profile['form']=='separate_products':v=operation(num,den,'div',model)
                predictions.append(bits(float(round_model(v,'strict'))))
            rows.append(dict(args=args,predictions=predictions))
    write(str(output)+'.questions.json',dict(format='dnacalc-hypergeom-diagonal-questions/1',complete=True,
        groups=groups,profiles=profiles,rows=rows,
        selection='Seven integer exponent edges, total=2^e or2^e+1, three fixed success counts and every diagonal sample1 through16. Nine literal product loops, source-code order, precede every answer.'))
    answers=service.ask('HYPGEOM.DIST',questions)['rows']
    aliases=service.ask('HYPGEOMDIST',[q[:4] for q,row in zip(questions,rows) if row['args'][0]==bits(1.)])['rows']
    counts=[Counter() for _ in profiles];paired=[Counter() for _ in profiles]
    records=[];last={};alias_index=0
    for row,answer in zip(rows,answers):
        n,_,k,total=(int(value(v)) for v in row['args'])
        pair=[];prior=last.get((k,total))
        for i,profile in enumerate(profiles):
            counts[i]['calls']+=1
            counts[i]['matched']+=answer==dict(type='number',bits=row['predictions'][i])
            counts[i]['failed']+=answer!=dict(type='number',bits=row['predictions'][i])
            domain=None
            if prior is not None and prior['type']=='number' and answer['type']=='number':
                # All three loops have a separately stated simple recurrence test.
                # This does not assume the full product program's seed producer.
                seed=BINARY64.pullback(Span(Q(value(prior['bits'])),Q(value(prior['bits']))))
                target=Span(Q(value(answer['bits'])),Q(value(answer['bits'])))
                model=profile['model'];a,b=Q(k-n+1),Q(total-n+1)
                if profile['form']=='ratio_then_multiply':
                    factor=operation(a,b,'div',model)
                    needed=inverse(target,factor,'mul',True,'strict' if model=='strict' else 'stored_pc64')
                elif profile['form']=='multiply_then_divide':
                    needed=inverse(target,b,'div',True,'strict' if model=='strict' else 'stored_pc64')
                    needed=inverse(needed,a,'mul',True,model) if needed else None
                else:
                    pair.append(None);continue
                domain=seed.meet(needed) if needed else None
                paired[i]['calls']+=1;paired[i]['excluded']+=domain is None
            pair.append(domain.json() if domain else None)
        record=dict(row,answer=answer,previous_diagonal=prior,recurrence_domains=pair)
        if n==1:
            record['legacy']=aliases[alias_index];alias_index+=1
        records.append(record);last[k,total]=answer
    result=dict(format='dnacalc-hypergeom-diagonal-products/1',complete=True,profiles=profiles,
        counts=counts,paired_counts=paired,records=records,receipts=service.receipts,
        questions_sha256=sha(str(output)+'.questions.json'),
        limits='Diagonal mass is mathematically the falling product of(K-i)/(N-i), without a binomial helper. '
               'All finite integer inputs are below2^53 and at most16 factors are evaluated. '
               'The mathematical identity supplies predictions only, never Excel answers. '
               'A complete loop needs every original answer, including errors. Recurrence cells '
               'allow an arbitrary real preceding value with the prior binary64 publication; '
               'an empty cell refutes only that recurrence. Resident full-loop and stored public '
               'recurrence models remain distinct. No program is selected by match count.')
    write(output,result);print(json.dumps(dict(counts=counts,paired_counts=paired)),flush=True)
    return result


def hypergeom_product_constraints(input_path,output):
    """Fixed small factorial-ratio trees, contradicted by recorded PMFs.

    This is a necessary-constraint search, with one counterexample per rejected
    complete program. Mathematical references never replace Excel readings.
    """
    import itertools
    rows=[json.loads(line)for line in Path(input_path).read_text().splitlines()]
    profiles=[dict(id=i,model=model,symmetric=symmetric,direction=direction,
                   combination=combination,outer=outer)
        for i,(model,symmetric,direction,combination,outer)in enumerate(itertools.product(
            ('binary64','stored_pc64','resident_pc64'),(True,False),
            ('ascending','descending'),('multiply_then_divide','ratio_then_multiply'),range(3)))]
    alive=set(range(len(profiles)));contradictions={};scope=Counter();mathematical=Counter();first_math=None;executions=0
    def rounded(x,model):
        return BINARY64.round(x)if model=='binary64'else(
            BINARY64.round(EXTENDED64.round(x))if model=='stored_pc64'else EXTENDED64.round(x))
    def publication(q):
        x=float(BINARY64.round(q))
        return '0x'+bits(0. if 0<abs(x)<sys.float_info.min else x)
    def choose(n,k,p):
        if k<0 or k>n:return Q(0)
        if p['symmetric']:k=min(k,n-k)
        v=Q(1);pairs=[(n-k+i,i)for i in range(1,k+1)]
        if p['direction']=='descending':pairs.reverse()
        for numerator,denominator in pairs:
            if p['combination']=='multiply_then_divide':
                v=rounded(rounded(v*numerator,p['model'])/denominator,p['model'])
            else:v=rounded(v*rounded(Q(numerator,denominator),p['model']),p['model'])
        return v
    for row in rows:
        args=row.get('args',[])
        arity=4 if row['function']=='HYPGEOMDIST'else 5
        if row['function']not in('HYPGEOM.DIST','HYPGEOMDIST')or len(args)!=arity or row.get('result_index')is not None:scope['surface']+=1;continue
        if not all(isinstance(v,str)and len(v)==18 and v.startswith('0x')for v in args[:4]):scope['non_numeric']+=1;continue
        xs=[value(v)for v in args[:4]]
        if not all(math.isfinite(x)and x>=0 and x==int(x)for x in xs):scope['outside_nonnegative_integer_domain']+=1;continue
        s,n,k,total=map(int,xs)
        if n>total or k>total:scope['invalid_support_parameters']+=1;continue
        if max(s,total)>10000:scope['research_resource_bound']+=1;continue
        flag=args[4]if len(args)==5 else False
        if isinstance(flag,bool):cumulative=flag
        elif isinstance(flag,str)and flag.startswith('0x')and len(flag)==18:cumulative=value(flag)!=0.0
        else:scope['non_numeric_cumulative_flag']+=1;continue
        if cumulative:scope['not_pmf']+=1;continue
        if not str(row['expected']).startswith('0x'):scope['non_numeric_expected']+=1;continue
        scope['scoped']+=1
        support=0<=s<=min(n,k)and 0<=n-s<=total-k
        numerator=math.comb(k,s)*math.comb(total-k,n-s)if support else 0
        reference=publication(Q(numerator,math.comb(total,n)))
        match=reference==row['expected'];mathematical['all']+=1;mathematical['exact']+=match
        if not match and first_math is None:first_math=dict(original=row,mathematical=reference)
        for identity in list(alive):
            p=profiles[identity];m=p['model'];executions+=1
            try:
                a,b,c=choose(k,s,p),choose(total-k,n-s,p),choose(total,n,p)
                candidate=(rounded(rounded(a*b,m)/c,m)if p['outer']==0 else
                    rounded(a*rounded(b/c,m),m)if p['outer']==1 else rounded(rounded(a/c,m)*b,m))
                candidate=publication(candidate)
            except ArithmeticError:candidate='refusal:nonfinite_combination_path'
            if candidate!=row['expected']:
                alive.remove(identity);contradictions[identity]=dict(original=row,prediction=candidate)
    result=dict(format='dnacalc-hypergeom-product-constraints/1',complete=True,input_rows=len(rows),scope=scope,
        profiles=profiles,contradictions=contradictions,survivors=sorted(alive),executions=executions,
        mathematical_control=dict(counts=mathematical,first_contradiction=first_math),
        input_sha256=sha(input_path),candidate_order='Three named coherent precisions; symmetry reduction; ascending then descending; multiply/divide then ratio/multiply; three written quotient associations. Fixed code-shape order, never match counts.',
        limits='Integer PMFs only, total/sample bounded10000 for research resources, not an Excel cutoff. Every source row is counted in scope/exclusion; exclusions are not passes. A first contradiction refutes a complete program, so later rows do not earn calculated-only pass rates. No coefficient, threshold or selector fitted. Original mathematical probability is a diagnostic, never an Excel answer. No product change follows without a coherent identified step.')
    write(output,result);print(json.dumps(dict(scope=scope,programs=len(profiles),survivors=len(alive),executions=executions)))


def append_question(book,name,formula,role,round_number=26,**fields):
    row=len(book)+1
    book.append(dict(id=f'r{round_number}-{name}-{row:05}',row=row,column=4,formula=formula,role=role,**fields))
    return 'D'+str(row)


def workbook_questions(checkout,count=64):
    """Questions for entries outside the existing scalar service, no answers.

    The original workbook generator and exact two-limb operand expression are
    reused. One scalar-coordinate window is separate from the regression book.
    """
    try:
        from .corpus_fuzz import exact_expression
    except ImportError:
        from corpus_fuzz import exact_expression
    client,policy=load_service(checkout)
    inventory=json.loads((Path(checkout)/'docs/generated/problem-functions.json').read_text())
    allowed=policy.catalogue_names(checkout)
    books={'scalar':[],'regression':[]};requests={};constants={}
    current='scalar'
    def emit(formula,role,**fields):
        return append_question(books[current],current,formula,role,**fields)
    def argument(arg):
        key=json.dumps(arg,sort_keys=True)
        bucket=constants.setdefault(current,{})
        if key not in bucket:
            if arg['type']=='number':formula=exact_expression(value(arg['bits']))
            elif arg['type']=='boolean':formula=str(arg['value']).upper()
            elif arg['type']=='text':formula='"'+arg['value'].replace('"','""')+'"'
            else:raise ValueError('unsupported workbook argument')
            bucket[key]=emit(formula,'operand_echo',intended_bits='0x'+arg['bits']if arg['type']=='number'else None,
                             intended_value=arg)
        return bucket[key]
    for entry in inventory['functions']:
        name=entry['function']
        if name in allowed or name in ('LINEST','LOGEST','TREND','GROWTH','MINVERSE'):continue
        args,index,origin=seed(name,checkout,client)
        if args is None or index is None:
            requests[name]=dict(status='no numerical scalar coordinate',reason=origin)
            continue
        rows=[];center=value(args[index]['bits'])
        for window,multiplier in enumerate((1.,1.25)):
            first=int(bits(center*multiplier),16)-count//2
            for b in range(first,first+count):
                actual=[dict(a)for a in args];actual[index]=client.number(f'{b:016x}')
                refs=list(map(argument,actual))
                call=name+'('+','.join(refs)+')'
                ref=emit(call,'question',function=name,window=window,input_bits=f'{b:016x}',varying_argument=index)
                rows.append(int(ref[1:]))
        requests[name]=dict(status='awaiting workbook',book=current,rows=rows,count_per_window=count,seed_origin=origin)
    # Reference/array functions have no generic numeric rounding coordinate.
    for name,formula in [('ROWS','ROWS(D1:D3)'),('COLUMNS','COLUMNS(D1:D3)'),
                         ('JIS','JIS("abc123")'),('TEXTSPLIT','INDEX(TEXTSPLIT("a,b",","),1,1)')]:
        emit(formula,'question',function=name,qualification='Structural/text control; no numerical stage-count claim.')
    current='regression'
    def vector(xs):
        # Keep each vector contiguous, even when a constant is already echoed
        # elsewhere. A range must contain precisely its intended observations.
        refs=[emit(exact_expression(x),'operand_echo',intended_bits='0x'+bits(x))for x in xs]
        return refs[0]+':'+refs[-1]
    for function in ('LINEST','LOGEST'):
        start=int('401f9ff8b1ae9622',16)-16
        for b in range(start,start+32):
            x=value(f'{b:016x}')
            xr=vector([x,0.,0.,.5])
            for basis in range(4):
                yr=vector([float(i==basis) if function=='LINEST'else (2. if i==basis else 1.)for i in range(4)])
                for constant in (False,True):
                    call=f'{function}({yr},{xr},{str(constant).upper()},TRUE)'
                    for row,column in ((1,1),(2,1),(4,2)):
                        emit(f'INDEX({call},{row},{column})','question',function=function,
                             purpose='Exact full-rank norm-store window, basis response, coefficient/SE/df together.',
                             input_bits=f'{b:016x}',basis=basis,constant=constant,result_index=[row,column])
    # New interleaved rank midpoints extend R25's unanswered grid. None of the
    # original decimal-array ranks is relabelled as an exact-reference answer.
    low,high=int('402000007ffffffd',16),int('40200001ffffffff',16)
    for function in ('LINEST','LOGEST'):
        for j in range(1,32,2):
            b=low+(high-low)*j//32
            y=vector([1.,2.,4.,7.]);a=vector([1.,2.,5.,4.]);c=vector([2.,4.,10.,value(f'{b:016x}')])
            for constant in (False,True):
                call=f'{function}({y},CHOOSE({{1,2}},{a},{c}),{str(constant).upper()},TRUE)'
                for row,column in ((4,2),(1,1),(1,2)):
                    emit(f'INDEX({call},{row},{column})','question',function=function,
                         purpose='Interleaved exact-reference rank bracket midpoint; response and df observed together.',
                         input_bits=f'{b:016x}',constant=constant,result_index=[row,column])
    for function in ('TREND','GROWTH'):
        y=vector([1.,2.,4.,7.]);x=vector([1.,2.,5.,4.])
        first=int(bits(3.7),16)-32
        for b in range(first,first+64):
            operand=emit(exact_expression(value(f'{b:016x}')),'operand_echo',intended_bits='0x'+f'{b:016x}')
            emit(f'INDEX({function}({y},{x},{operand},TRUE),1,1)','question',function=function,
                 input_bits=f'{b:016x}',purpose='Numeric selected regression output at consecutive prediction inputs.')
    for b in range(int(bits(1.3),16)-32,int(bits(1.3),16)+32):
        a=vector([value(f'{b:016x}'),1.]);c=vector([1.,2.])
        emit(f'INDEX(MINVERSE(CHOOSE({{1,2}},{a},{c})),1,1)','question',function='MINVERSE',
             input_bits=f'{b:016x}',purpose='Selected inverse of a nonsingular two-by-two matrix, exact input echoes.')
    return books,requests


def root_transfer_questions(count=16):
    """Unanswered range/array checks around the observed scalar arity change.

    All routes share explicit source operands. A decimal-array INDEX echo
    measures that parser separately; no intended decimal value is an answer.
    Separate LN, reciprocal, mean and EXP cells expose publication boundaries.
    """
    try:
        from .corpus_fuzz import exact_expression
    except ImportError:
        from corpus_fuzz import exact_expression
    book=[]
    def emit(formula,role,**fields):return append_question(book,'root_transfer',formula,role,**fields)
    for n in (9,10,11,32):
        for center in (1.3,1000.):
            first=int(bits(center),16)-count//2
            for raw in range(first,first+count):
                input_bits=f'{raw:016x}';p=value(input_bits)
                group=dict(count=n,input_bits=input_bits)
                refs=[emit(exact_expression(p),'operand_echo',intended_bits='0x'+input_bits,**group)]
                refs += [emit('1','operand_echo',intended_bits='0x3ff0000000000000',**group)for _ in range(n-1)]
                array='{'+','.join([format(p,'.17g')]+['1']*(n-1))+'}'
                echo=emit('INDEX('+array+',1,1)','operand_echo',intended_array_operand_bits='0x'+input_bits,**group)
                for route,arguments in [('range',refs[0]+':'+refs[-1]),('scalar_cells',','.join(refs)),('decimal_array',array)]:
                    emit('GEOMEAN('+arguments+')','question',function='GEOMEAN',route=route,
                         source_cells=refs,array_operand_echo=echo if route=='decimal_array'else None,
                         purpose='Compare input routes at fixed arities; scalar agreement is not range/array evidence.',**group)
                logarithm=emit('LN('+refs[0]+')','question',function='LN',route='stored_helper_leaf',**group)
                reciprocal=emit('1/'+str(n),'arithmetic_control',operation='/',route='stored_reciprocal',**group)
                mean=emit(logarithm+'*'+reciprocal,'arithmetic_control',operation='*',route='stored_mean',**group)
                emit('EXP('+mean+')','question',function='EXP',route='staged_composition',
                     logarithm=logarithm,reciprocal=reciprocal,mean=mean,**group)
                emit('POWER('+refs[0]+','+reciprocal+')','question',function='POWER',route='stored_reciprocal_control',**group)
    return {'root_transfer':book},{}


def rank_dependency_questions():
    """Vary observation count, pivot order and predictor count at measured R25 edges.

    The exact residual geometry is unchanged by zero rows without an intercept,
    or by repeating observations. No tolerance is inferred from these questions.
    """
    try:
        from .corpus_fuzz import exact_expression
    except ImportError:
        from corpus_fuzz import exact_expression
    intake=ROOT/'evidence/round27/rank-intake.json'
    data=json.loads(intake.read_text());observed={r['id']:r for r in data['observations']}
    book=[]
    def emit(formula,role,**fields):
        return append_question(book,'rank_dependencies',formula,role,round_number=27,**fields)
    vectors={}
    def vector(xs,**fields):
        identity=tuple(bits(x)for x in xs)
        if identity not in vectors:
            refs=[emit(exact_expression(x),'operand_echo',intended_bits='0x'+bits(x),**fields)for x in xs]
            vectors[identity]=refs[0]+':'+refs[-1]
        return vectors[identity]
    seeds=[]
    for interval in data['necessary_intervals']:
        if interval['function']!='LINEST'or interval['reverse']or interval['normalizer']!='column':continue
        for side in ('lower_witness','upper_witness'):
            observation=observed[interval[side]];seeds.append(observation)
            for function in ('LINEST','LOGEST'):
                for design in ('original','reverse','pad1','pad4','repeat2','repeat4','third_orthogonal'):
                    for scales in ((1.,1.),(2.**20,1.),(1.,2.**20)):
                        for swap in (False,True):
                            # Those exact reference controls are already answered
                            # in R25. Reuse their recorded facts, never ask again.
                            if function=='LINEST'and design in ('original','reverse')and scales==(1.,1.)and not swap:
                                continue
                            a=[1.,2.,5.,4.];b=[2.,4.,10.,value(observation['changed_operand'])]
                            columns=[[x*scales[0]for x in a],[x*scales[1]for x in b]]
                            response=[1.,2.,4.,7.]if function=='LINEST'else[2.,4.,16.,128.]
                            if swap:columns.reverse()
                            if design=='third_orthogonal':columns.append([3.,-4.,1.,0.])
                            if design=='reverse':
                                columns=[col[::-1]for col in columns];response.reverse()
                            if design.startswith('pad'):
                                count=int(design[3:]);columns=[col+[0.]*count for col in columns]
                                response += [0. if function=='LINEST'else 1.]*count
                            if design.startswith('repeat'):
                                count=int(design[6:]);columns=[col*count for col in columns];response*=count
                            context=dict(seed=observation['id'],changed_operand=observation['changed_operand'],
                                design=design,scales=scales,swap=swap,constant=observation['constant'],function=function)
                            y=vector(response,**context);xs=[vector(col,**context)for col in columns]
                            x='CHOOSE({'+','.join(str(i+1)for i in range(len(xs)))+'},'+','.join(xs)+')'
                            call=f'{function}({y},{x},{str(observation["constant"]).upper()},TRUE)'
                            for row,column in ((4,2),(1,1),(1,len(columns))):
                                emit(f'INDEX({call},{row},{column})','question',result_index=[row,column],
                                    purpose='Separate row-count, pivot-order, column-scale and predictor-count rank decisions; coefficients remain recorded together.',**context)
    return {'rank_dependencies':book},dict(intake_sha256=sha(intake),seeds=seeds,
        limits='Original measured drop/keep bracket endpoints choose inputs only. New reference/formula results all await Excel. Padding with an intercept changes centered geometry; it is a separate question, not an assumed invariant.')


def rank_bisection_questions(intake):
    """Ask Excel to bisect its own rank edge using exact integer offsets.

    A bit offset from8 is below2^53, so Excel can halve the integer interval
    exactly. Each decision uses the preceding observed degrees of freedom.
    The generated workbook carries questions and no predicted answers.
    """
    intake=Path(intake)
    data=json.loads(intake.read_text());observed={r['id']:r for r in data['observations']}
    book=[];seeds=[]
    def emit(formula,role,**fields):
        return append_question(book,'rank_edges',formula,role,round_number=28,**fields)
    def vector(values,**fields):
        refs=[emit(str(v),'operand_echo',intended_bits='0x'+bits(v),**fields)for v in values]
        return refs[0]+':'+refs[-1]
    for interval in data['necessary_intervals']:
        if interval['normalizer']!='column'or interval['reverse']:continue
        lower=observed[interval['lower_witness']];upper=observed[interval['upper_witness']]
        low,high=(int(r['changed_operand'].removeprefix('0x'),16)for r in (lower,upper))
        anchor=0x4020000000000000
        if not anchor<low<high<anchor+2**40:raise ValueError('rank bracket outside exact offset construction')
        if high-low>2**36:raise ValueError('rank bracket is too wide for36 bisections')
        name,constant=interval['function'],interval['constant']
        group=dict(function=name,constant=constant)
        seeds.append(dict(group,low_witness=lower,high_witness=upper))
        a=vector([1,2,5,4],**group);y=vector([1,2,4,7],**group)
        lo=emit(str(low-anchor),'operand_echo',intended_bits='0x'+bits(low-anchor),**group)
        hi=emit(str(high-anchor),'operand_echo',intended_bits='0x'+bits(high-anchor),**group)
        def call(offset,phase):
            b=[emit(str(v),'operand_echo',intended_bits='0x'+bits(v),phase=phase,**group)for v in (2,4,10)]
            b.append(emit('8+'+offset+'/2^49','operand_echo',phase=phase,
                          purpose='Effective adaptive rank operand; its bits await Excel.',**group))
            expression=f'{name}({y},CHOOSE({{1,2}},{a},{b[0]}:{b[-1]}),{str(constant).upper()},TRUE)'
            df=emit('INDEX('+expression+',4,2)','question',phase=phase,result_index=[4,2],
                    purpose='Observed rank predicate for an adaptive integer-offset bisection.',**group)
            return df,expression
        call(lo,'initial_low');call(hi,'initial_high')
        drop,keep=lower['df'],upper['df']
        if drop!=keep+1:raise ValueError('rank bracket does not distinguish one predictor')
        for step in range(36):
            mid=emit(lo+'+INT(('+hi+'-'+lo+')/2)','arithmetic_control',step=step,operation='integer_midpoint',**group)
            df,_=call(mid,'bisect_'+str(step))
            lo=emit(f'IF({df}={drop},{mid},IF({df}={keep},{lo},NA()))','arithmetic_control',step=step,operation='lower_bracket',**group)
            hi=emit(f'IF({df}={keep},{mid},IF({df}={drop},{hi},NA()))','arithmetic_control',step=step,operation='upper_bracket',**group)
        width=emit(hi+'-'+lo,'arithmetic_control',operation='final_bracket_width',**group)
        for offset,phase in ((lo,'final_low'),(hi,'final_high')):
            df,expression=call(offset,phase)
            for column in (1,2):
                emit(f'INDEX({expression},1,{column})','question',phase=phase,result_index=[1,column],
                     final_width=width,purpose='Coefficient together with the observed adjacent rank edge.',**group)
    return {'rank_edges':book},dict(intake_sha256=sha(intake),seeds=seeds,
        limits='Four rank brackets from original R26 observations. Each Excel midpoint is chosen from its own preceding measured df. Local monotonicity is a premise to check against the returned trace.36 exact integer bisections suffice for each declared bracket; final endpoints and width are separate observations. No tolerance, QR producer or coefficient is fitted.')


def make_workbooks(checkout,output,executable,powershell,root_transfer=False,rank_dependencies=False,rank_edges=None):
    checkout,output,executable,powershell=map(lambda p:Path(p).resolve(),(checkout,output,executable,powershell))
    output.mkdir(parents=True,exist_ok=False)
    if sum(map(bool,(root_transfer,rank_dependencies,rank_edges)))>1:raise ValueError('choose one workbook study')
    round_name='r28'if rank_edges else'r27'if rank_dependencies else'r26'
    books,requests=(rank_bisection_questions(rank_edges)if rank_edges else rank_dependency_questions()if rank_dependencies else
                    root_transfer_questions()if root_transfer else workbook_questions(checkout))
    generator=checkout/'crates/calc-func/tests/fixtures/make-catalogue-fixtures.ps1'
    provenance=subprocess.check_output([str(executable),'provenance'],timeout=60)
    (output/'prediction-binding.json').write_bytes(provenance)
    inventory=[]
    for name,questions in books.items():
        directory=output/name;directory.mkdir()
        case='dcj9k_'+round_name+'_'+name
        doc=directory/'workbook-questions.json';write(doc,dict(sheet=case,allow_refusals=True,questions=questions))
        section=directory/('sheet-'+case+'.txt')
        section.write_text('# '+round_name.upper()+' exact operands and unanswered fingerprint questions.\n'+'\n'.join(q['formula']for q in questions)+'\n')
        predictions=directory/'application-predictions.jsonl'
        # Existing registered Workbook prediction route. These results remain
        # predictions and are never passed to Excel as expected answers.
        subprocess.run([str(executable),'sheet-predictions',str(doc),str(predictions)],check=True,timeout=600)
        rows=[json.loads(line)for line in predictions.read_text().splitlines()]
        if len(rows)!=len(questions):raise ValueError('missing packet predictions')
        for q,r in zip(questions,rows):
            if q['id']!=r['id']or q['formula']!=r['formula']:raise ValueError('packet prediction identity changed')
            if q.get('intended_bits') and r['calc_func_prediction']!=q['intended_bits']:
                raise ValueError('prediction route changes an exact operand')
        subprocess.run([str(powershell),'-NoProfile','-File',str(generator),'-OutDir',str(directory/'workbooks'),
                        '-QuestionSectionFile',str(section),'-QuestionCase',case],check=True,timeout=300)
        inventory.append(dict(name=name,case=case,formula_cells=len(questions),
                              selected_calls=sum(q['role']=='question'for q in questions),
                              workbook=str((directory/'workbooks'/case/'oracle.xlsx').relative_to(output))))
    manifest=dict(format='dnacalc-'+round_name+'-fingerprint-questions/1',new_excel_answers=0,workbooks=inventory,
                  scalar_requests=requests,generator_sha256=sha(generator),producer_sha256=sha(__file__),
                  executable_sha256=sha(executable),
                  files_sha256={str(p.relative_to(output)):sha(p)for p in output.rglob('*')if p.is_file()},
                  limits=('New workbook only. Reference-cell inputs use exact integer limbs and echoes. Decimal arrays have their own INDEX operand echo, without an assumed binary64 value. LN, reciprocal, mean and EXP occupy separate cells. Every workbook outcome awaits capture; scalar observations do not establish input-route transfer. Predictions are application calculations, not Excel readings.'if root_transfer else
                          'New workbooks only. Every numerical operand uses exact small integer limbs and is echoed. Predictions are application calculations, not Excel readings. QR rank thresholds and complete arithmetic remain unknown; no tolerance or coefficient fitted.'))
    write(output/'manifest.json',manifest)
    if rank_edges:
        (output/'README.txt').write_text('R28 unanswered adaptive rank bisection workbook.\nCapture only the new workbook listed in manifest.json through the existing capture route.\nExcel chooses exact integer-offset midpoints from its own preceding degrees-of-freedom answers.\nReturn every helper cell, input echo, df, coefficient, final bracket width and original worker/canonical capture.\nPredictions use the application route and are never Excel answers. No rank tolerance is inferred or fitted.\n')
    elif rank_dependencies:
        (output/'README.txt').write_text('R27 unanswered rank dependency workbook.\nCapture only the new workbook listed in manifest.json through the existing capture route.\nR25 measured bracket endpoints are reused as inputs; no new answer is supplied.\nSeparate padding, repeated observations, reversed rows, predictor swap, exact column scaling and a third orthogonal predictor.\nAll numerical operands use exact integer limbs and have unchanged bit echoes.\nReturn canonical captures, raw worker values/settings and original workbook bytes.\n')
    elif root_transfer:
        (output/'README.txt').write_text('R26 unanswered GEOMEAN input-route and stored-mean workbook.\nCapture only the workbook listed in manifest.json through the existing cargo xtask capture/compare route.\nCompare range, scalar-cell and decimal-array GEOMEAN calls at arities9,10,11,32.\nThe decimal-array operand has its own INDEX echo; its intended bits are not assumed.\nLN, reciprocal, mean and EXP occupy separate cells; POWER uses the stored reciprocal.\nAll workbook outcomes await capture. Existing scalar readings do not answer these route questions.\nReturn canonical captures, raw worker values/settings and original workbook bytes.\n')
    else:
        (output/'README.txt').write_text('R26 unanswered scalar fingerprint and regression workbooks.\nCapture only the two workbooks listed in manifest.json through the existing cargo xtask capture/compare route.\nReturn complete canonical captures, raw worker values/settings, operand echoes and original workbook bytes.\nThe scalar book covers service-unavailable entries; text/structural controls do not identify numerical stages.\nThe regression book separates full-rank norm stores from rank-bracket decisions. No new answer is supplied.\n')
    print(json.dumps(inventory),flush=True)
    return manifest


def yield_single_coupon(checkout, output, cache, offline=False, month_ends=False, date_holdout=False):
    """Separate the one-coupon YIELD day fraction with public coupon calls."""
    import datetime as dt
    import calendar
    from exact import operation
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    serial=lambda date:(date-dt.date(1899,12,30)).days
    states=[]
    maturities=([dt.date(2000,2,29),dt.date(2001,3,30),dt.date(2024,10,30),dt.date(2025,3,15),dt.date(2031,11,30)]
                if date_holdout else[dt.date(2023,2,28),dt.date(2024,2,29),dt.date(2025,5,31),dt.date(2020,8,31),dt.date(2021,1,31)]
                if month_ends else[dt.date(2021,1,1)])
    for maturity in maturities:
      for frequency in (1,2,4):
        month_index=maturity.year*12+maturity.month-1-12//frequency
        year,month=divmod(month_index,12);month+=1
        at_end=maturity.day==calendar.monthrange(maturity.year,maturity.month)[1]
        day=calendar.monthrange(year,month)[1]if at_end else min(maturity.day,calendar.monthrange(year,month)[1])
        previous=dt.date(year,month,day)
        days=(maturity-previous).days
        offsets={0,1,days//2,*range(days-10,days)}
        for month_index in range(previous.year*12+previous.month-1,maturity.year*12+maturity.month):
            year,month=divmod(month_index+1,12);next_month=dt.date(year,month+1,1)
            date=next_month-dt.timedelta(days=1)
            if previous<=date<maturity:offsets.add((date-previous).days)
        for offset in sorted(offsets):
            for basis in range(5):
                states.append(dict(settlement=serial(previous+dt.timedelta(days=offset)),maturity=serial(maturity),
                                   frequency=frequency,basis=basis))
    rates=(0.,math.ldexp(1.,-20),math.nextafter(.0625,math.inf))if date_holdout else(math.ldexp(1.,-20),.03125,math.nextafter(.05,0.),math.nextafter(.05,math.inf),.125)if month_ends else(0.,.04125,.05)
    prices=(.5,math.nextafter(123.456,math.inf))if date_holdout else(64.,math.nextafter(100.,0.),math.nextafter(100.,math.inf),128.)if month_ends else(95.,100.,105.)
    redemptions=(80.,128.)if date_holdout else(100.,)
    cases=[dict(state,rate=rate,price=price,redemption=red)for state in states for rate in rates for price in prices for red in redemptions]
    source=Path(checkout)/'crates/calc-func/src/ported/bond_core_family.rs'
    frozen={str(source.resolve()):sha(source),str(Path(__file__).resolve()):sha(__file__)}
    snapshots=Path(output).resolve().parent/'producers';snapshots.mkdir(exist_ok=True)
    for path,digest in frozen.items():
        dest=snapshots/(digest+'-'+Path(path).name)
        if not dest.exists():dest.write_bytes(Path(path).read_bytes())
    write(str(output)+'.questions.json',dict(format='dnacalc-yield-single-coupon-questions/1',complete=True,
        states=states,cases=cases,source_bindings=frozen,
        models=['direct published day span','derived E-A','coupon-date E otherwise public span'],
        selection='Fixed final coupon periods, all five bases and three frequencies, coupon/month ends '
                  'and the last ten settlement days. Rates, prices and redemptions are listed explicitly. '
                  'Public COUP factors are observed before the YIELD targets. No stopping tolerance or coefficient fitted.',
        month_end_holdout=month_ends,date_holdout=date_holdout,rates=rates,prices=prices,redemptions=redemptions,maturities=[str(d)for d in maturities],
        mathematical_source='https://support.microsoft.com/en-us/excel/functions/yield-function'))
    factor_questions=[[client.number(s[k])for k in('settlement','maturity','frequency','basis')]for s in states]
    factors={name:service.ask(name,factor_questions)['rows']for name in('COUPNUM','COUPDAYS','COUPDAYBS','COUPDAYSNC')}
    for readings in factors.values():
        if any(r['type']!='number'for r in readings):raise ValueError('coupon controls must be numeric')
    if any(value(r['bits'])!=1. for r in factors['COUPNUM']):raise ValueError('declared final-coupon scope contradicted')
    questions=[[client.number(c[k])for k in('settlement','maturity','rate','price','redemption','frequency','basis')]for c in cases]
    answers=service.ask_rows('YIELD',questions)
    def model(c,a,e,dsc,word):
        if not dsc:return None
        op=lambda x,y,k:operation(Q(x),Q(y),k,word)
        coupon=op(op(100,c['rate'],'mul'),c['frequency'],'div')
        accrued=op(op(coupon,a,'mul'),e,'div')
        ratio=op(op(c['redemption'],coupon,'add'),op(c['price'],accrued,'add'),'div')
        term=op(ratio,-1,'add')
        term=op(op(term,c['frequency'],'mul'),e,'mul')
        return bits(float(op(term,dsc,'div')))
    records=[];counts=Counter();price_questions=[];price_indices=[]
    for index,(case,answer)in enumerate(zip(cases,answers)):
        state_index=index//(len(rates)*len(prices)*len(redemptions))
        a,e,dsc=[value(factors[name][state_index]['bits'])for name in('COUPDAYBS','COUPDAYS','COUPDAYSNC')]
        predictions={name+'_'+word:model(case,a,e,span,word)
                     for name,span in(('public_span',dsc),('derived_span',e-a),('coupon_date',e if a==0 else dsc))
                     for word in('strict','stored_pc64')}
        for label,prediction in predictions.items():counts[label+'_matches']+=answer.get('type')=='number'and answer['bits']==prediction
        counts['rows']+=1;counts['nonnumeric_targets']+=answer['type']!='number'
        counts['all_finite_models_refuted']+=answer.get('bits')not in[p for p in predictions.values()if p is not None]
        record=dict(case,answer=answer,factors={name:rows[state_index]for name,rows in factors.items()},predictions=predictions)
        records.append(record)
        if answer['type']=='number'and value(answer['bits'])>=0:
            args=[client.number(case[k])for k in('settlement','maturity','rate')]
            args += [client.number(answer['bits'])]
            args += [client.number(case[k])for k in('redemption','frequency','basis')]
            price_questions.append(args);price_indices.append(index)
    prices=service.ask_rows('PRICE',price_questions)
    for index,price in zip(price_indices,prices):
        records[index]['price_at_excel_yield']=price
        counts['price_controls']+=1
        counts['price_at_yield_matches_target']+=price.get('type')=='number'and price['bits']==bits(records[index]['price'])
    result=dict(format='dnacalc-yield-single-coupon/1',complete=True,records=records,counts=counts,
        receipts=service.receipts,source_bindings=frozen,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Public coupon factors and single-period closed-form programs only. Stored arithmetic words '
               'remain separate; a nonfinite model is not scored as a typed-error match. PRICE refuses negative '
               'input yields, so those YIELD observations have no paired PRICE call. This does not identify '
               'the multi-period Newton solver, its objective, or an unobserved day-count helper.')
    write(output,result);print(json.dumps(counts),flush=True)
    return result


def yield_coupon_days(checkout, original, output, cache, offline=False):
    """Observe direct day counts separately from the published coupon counts."""
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    original=Path(original);study=json.loads(original.read_text());records=study['records']
    states=list({tuple(r[k]for k in('settlement','maturity','frequency','basis')):None for r in records})
    arguments=[[client.number(v)for v in state]for state in states]
    dates={name:service.ask(name,arguments)['rows']for name in('COUPPCD','COUPNCD')}
    require=lambda condition:None if condition else (_ for _ in()).throw(ValueError('coupon-date control failed'))
    require(all(r['type']=='number'for rows in dates.values()for r in rows))
    requests=[];slots=[];direct={};state_dates={}
    for index,state in enumerate(states):
        settlement,maturity,frequency,basis=state
        prev,next_=[int(value(dates[name][index]['bits']))for name in('COUPPCD','COUPNCD')]
        require(next_==maturity and prev<=settlement<maturity)
        state_dates[state]={name:rows[index]for name,rows in dates.items()}
        direct[state]={}
        for label,start,end in(('a',prev,settlement),('e',prev,next_),('dsr',settlement,next_)):
            if basis in(0,4):
                slots.append((state,label));requests.append([client.number(start),client.number(end),client.boolean(basis==4)])
            else:direct[state][label]=dict(type='number',bits=bits(end-start))
    spans=service.ask('DAYS360',requests)['rows'];require(all(r['type']=='number'for r in spans))
    for(state,label),answer in zip(slots,spans):direct[state][label]=answer
    for record in records:
        state=tuple(record[k]for k in('settlement','maturity','frequency','basis'))
        record['direct_day_counts']=direct[state]
        record['coupon_dates']=state_dates[state]
    result=dict(format='dnacalc-yield-direct-day-counts/1',complete=True,records=records,
        original_sha256=sha(original),receipts=service.receipts,parent_receipts=study['receipts'],
        limits='Direct DAYS360 for bases 0/4, actual serial-date differences otherwise, with '
               'independent public previous/next coupon dates. Distinct from the published COUP spans. '
               'Original YIELD and PRICE readings are unchanged.')
    write(output,result);print(json.dumps(dict(rows=len(records),states=len(states),days360_calls=len(requests))),flush=True)
    return result


def yield_coupon_programs(original, output, library, admitted=None, day_boundaries=False, nasd_days=False,numerator_days=False):
    """Test the documented normalized formula and its elementary associations."""
    original=Path(original);study=json.loads(original.read_text());native=Native(library)
    actual={tuple(r[k]for k in('settlement','maturity','frequency')):r['factors']['COUPDAYS']['bits']
            for r in study['records']if r['basis']==1}
    has_direct=all('direct_day_counts'in r for r in study['records'])
    templates=[]
    for normalization in ('percent','face'):
      for accrual in ('fraction_first','coupon_first'):
       for numerator in ('difference','ratio_minus_one'):
        for period in (('published','actual_for_2_3','direct')if has_direct else('published','actual_for_2_3')):
         for remaining in (('numerator_direct',)if numerator_days else('nasd_direct',)if nasd_days else('published','complement','direct','coupon_or_direct')if day_boundaries else('published','complement','direct')if has_direct else('published','complement')):
          for association in ('factor_first','left_products','fraction_first','combined_denominator'):
            g=Builder(8);rate,price,red,frequency,a,e,dsc,actual_e=range(8)
            if period in('actual_for_2_3','direct'):e=actual_e
            hundred=g.constant(bits(100.))
            if normalization=='percent':
                price=g.op('div',price,hundred);red=g.op('div',red,hundred)
                coupon=g.op('div',rate,frequency)
            else:coupon=g.op('div',g.op('mul',hundred,rate),frequency)
            if remaining=='complement':dsc=g.op('sub',e,a)
            if accrual=='fraction_first':accr=g.op('mul',g.op('div',a,e),coupon)
            else:accr=g.op('div',g.op('mul',coupon,a),e)
            dirty=g.op('add',price,accr);cash=g.op('add',red,coupon)
            if numerator=='difference':term=g.op('div',g.op('sub',cash,dirty),dirty)
            else:term=g.op('sub',g.op('div',cash,dirty),g.one)
            if association=='factor_first':result=g.op('mul',term,g.op('div',g.op('mul',frequency,e),dsc))
            elif association=='left_products':result=g.op('div',g.op('mul',g.op('mul',term,frequency),e),dsc)
            elif association=='fraction_first':result=g.op('div',g.op('mul',term,frequency),g.op('div',dsc,e))
            else:result=g.op('div',term,g.op('div',dsc,g.op('mul',frequency,e)))
            templates.append(dict(normalization=normalization,accrual=accrual,numerator=numerator,
                                  period=period,remaining=remaining,association=association,graph=g.graph(result,result)))
    rows=[];programs=[]
    for r in study['records']:
        factor=lambda name:int(r['factors'][name]['bits'],16)
        actual_e=actual[tuple(r[k]for k in('settlement','maturity','frequency'))]if r['basis']in(2,3)else r['factors']['COUPDAYS']['bits']
        rows.append(([int(bits(r[k]),16)for k in('rate','price','redemption','frequency')]+
                     [factor('COUPDAYBS'),factor('COUPDAYS'),factor('COUPDAYSNC'),int(actual_e,16)],
                     r['answer']))
    allowed=None
    if admitted is not None:
        admitted_programs=json.loads(Path(admitted).read_text())['complete_programs']
        allowed=set()
        for p in admitted_programs:
            template=dict(p['template'])
            if day_boundaries:template['remaining']='coupon_or_direct'
            if nasd_days:template['remaining']='nasd_direct'
            if numerator_days:template['remaining']='numerator_direct'
            allowed.add(json.dumps(template,sort_keys=True)+'/'+p['word'])
    if allowed is not None and not allowed:raise ValueError('no complete prior programs to validate')
    for template in templates:
      graph=template['graph']
      for word,cw,mask in [('strict',0x123f,(1<<graph['sites'])-1),('stored_pc64',0x133f,(1<<graph['sites'])-1),('resident_pc64',0x133f,0)]:
        if allowed is not None and json.dumps({k:v for k,v in template.items()if k!='graph'},sort_keys=True)+'/'+word not in allowed:continue
        count=0;first=None
        for index,(original_operands,answer)in enumerate(rows):
            operands=list(original_operands);record=study['records'][index]
            if template['period']=='direct':
                operands[4]=int(record['direct_day_counts']['a']['bits'],16)
                operands[7]=int(record['direct_day_counts']['e']['bits'],16)
            if template['remaining']in('direct','coupon_or_direct','nasd_direct','numerator_direct'):operands[6]=int(record['direct_day_counts']['dsr']['bits'],16)
            if(nasd_days or numerator_days)and record['basis']==0:
                import datetime as dt,calendar
                s,m=[dt.date(1899,12,30)+dt.timedelta(days=int(record[k]))for k in('settlement','maturity')]
                last_feb=lambda d:d.month==2 and d.day==calendar.monthrange(d.year,2)[1]
                if last_feb(s)and last_feb(m):operands[6]=int(bits(value(f'{operands[6]:016x}')+30-m.day),16)
                elif numerator_days and last_feb(s)and m.day==31:operands[6]=int(bits(value(f'{operands[6]:016x}')+1),16)
            if template['remaining']=='coupon_or_direct'and record['settlement']==value(record['coupon_dates']['COUPPCD']['bits']):
                operands[6]=operands[7]if template['period']in('actual_for_2_3','direct')else operands[5]
            predicted=('0000000000000000'if(day_boundaries or nasd_days or numerator_days)and operands[6]==0 else
                       f'{native.publish(native.evaluate(graph,operands,cw,mask)[:10]):016x}')
            matches=answer['type']=='number'and answer['bits']==predicted
            count+=matches
            if not matches and first is None:first=dict(index=index,expected=answer,predicted=predicted)
        programs.append(dict(template={k:v for k,v in template.items()if k!='graph'},word=word,cw=cw,mask=mask,
                             matches=count,first_failure=first))
    bindings={str(Path(__file__).resolve()):sha(__file__),str(ROOT/'racer/x87_programs.py'):sha(ROOT/'racer/x87_programs.py'),
              str(ROOT/'racer/x87_native.c'):sha(ROOT/'racer/x87_native.c')}
    snapshots=Path(output).resolve().parent/'producers';snapshots.mkdir(exist_ok=True)
    for path,digest in bindings.items():
        dest=snapshots/(digest+'-'+Path(path).name)
        if not dest.exists():dest.write_bytes(Path(path).read_bytes())
    result=dict(format='dnacalc-yield-coupon-programs/1',complete=True,original_sha256=sha(original),rows=len(rows),
                templates=templates,programs=programs,complete_programs=[p for p in programs if p['matches']==len(rows)],
                source_bindings=bindings,library_sha256=sha(library),admitted_programs_sha256=sha(admitted)if admitted else None,
                day_boundaries=day_boundaries,day_boundary_rule='Use E on the previous coupon date; a zero remaining direct day span publishes zero.'if day_boundaries else None,
                nasd_days=nasd_days,nasd_rule='For US basis only, when both original endpoints are February ends, adjust the final day to30. Zero remaining-day span publishes zero.'if nasd_days else None,
                numerator_days=numerator_days,numerator_rule='Use the existing US numerator-position day rule: adjust two February ends, but test end-day31 against the original start day before adjusting a February start. Zero remaining-day span publishes zero.'if numerator_days else None,
                limits='Fixed elementary associations of the public documented formula. Actual coupon length '
                'for bases 2/3 is independently taken from paired basis-1 COUPDAYS. No multi-period solver '
                'or different coupon date convention is inferred. Programs are not ranked by partial matches.')
    write(output,result);print(json.dumps(dict(rows=len(rows),programs=len(programs),complete=result['complete_programs'])),flush=True)
    return result


def duration_zero_coupon(checkout,output,cache,library,offline=False,count=64,holdout=False):
    """Separate the weighted redemption and final ratio on exact coupon dates."""
    import datetime as dt
    from exact import operation
    if not 8<=count<=1024:raise ValueError('bounded adjacent window required')
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library);states=[];cases=[]
    settlement=43831
    ns=(2,6,10,11,19,23,29)if holdout else(3,5,7,9,13,17,31)
    for frequency in(1,2,4):
      for n in ns:
        months=2020*12+n*(12//frequency);year,month=divmod(months,12)
        maturity=(dt.date(year,month+1,1)-dt.date(1899,12,30)).days
        state=dict(settlement=settlement,maturity=maturity,frequency=frequency,n=n)
        states.append(state)
        for center in(.0001,.005,.05,.125,.25,.5,.75,1.):
            anchor=int(bits(center),16)+(8192 if holdout else 0)
            for step in range(-count//2,count-count//2):
                yld=value(f'{anchor+step:016x}');q=1.+yld/frequency
                cases.append(dict(state,yield_bits=bits(yld),base_strict=bits(q),
                    base_stored_pc64=bits(float(operation(Q(1),Q(yld/frequency),'add','stored_pc64')))))
    source=Path(checkout)/'crates/calc-func/src/ported/bond_core_family.rs'
    bindings={str(source.resolve()):sha(source),str(Path(__file__).resolve()):sha(__file__),
              str(ROOT/'racer/x87_programs.py'):sha(ROOT/'racer/x87_programs.py'),
              str(ROOT/'racer/x87_native.c'):sha(ROOT/'racer/x87_native.c')}
    snapshots=Path(output).resolve().parent/'producers';snapshots.mkdir(exist_ok=True)
    for path,digest in bindings.items():
        dest=snapshots/(digest+'-'+Path(path).name)
        if not dest.exists():dest.write_bytes(Path(path).read_bytes())
    write(str(output)+'.questions.json',dict(format='dnacalc-duration-zero-coupon-questions/1',complete=True,
        states=states,cases=cases,holdout=holdout,count=count,source_bindings=bindings,
        selection='Fixed odd/even coupon counts, all three frequencies and eight ordinary/exponent-edge yield windows. '
        'Zero coupon removes the entire coupon sum. Public POWER, PRICE, SUM and coupon factors precede targets. '
        'Face1/100 and weighted-numerator/weighted-term/linear-time programs are fixed before target answers.'))
    factors=[[client.number(s[k])for k in('settlement','maturity','frequency')]+[client.number(0.)]for s in states]
    for name in('COUPNUM','COUPDAYBS'):
        answers=service.ask(name,factors)['rows']
        for state,answer in zip(states,answers):
            wanted=state['n']if name=='COUPNUM'else 0
            if answer!=dict(type='number',bits=bits(wanted)):raise ValueError('on-coupon schedule control changed')
    sums=service.ask_rows('SUM',[[client.number(1.),client.number(value(c['yield_bits'])/c['frequency'])]for c in cases])
    if any(r!=dict(type='number',bits=c['base_stored_pc64'])for c,r in zip(cases,sums)):
        raise ValueError('public base addition differs from its independent PC64 prediction')
    power_requests=[];power_slots=[]
    for c in cases:
        slots={}
        for prefix in('strict','stored_pc64'):
            base=c['base_'+prefix]
            if prefix=='stored_pc64'and base==c['base_strict']:slots[prefix]=slots['strict']
            else:
                slots[prefix]=len(power_requests)
                power_requests.append([client.number(base),client.number(c['n'])])
        power_slots.append(slots)
    powers=service.ask_rows('POWER',power_requests)
    if any(r['type']!='number'for r in powers):raise ValueError('nonnumeric discount-factor control')
    duration_args=[[client.number(c[k])for k in('settlement','maturity')]+[client.number(0.),client.number(c['yield_bits']),
                   client.number(c['frequency']),client.number(0.)]for c in cases]
    price_args=[args[:4]+[client.number(100.)]+args[4:]for args in duration_args]
    prices=service.ask_rows('PRICE',price_args)
    duration=service.ask_rows('DURATION',duration_args)
    modified=service.ask_rows('MDURATION',duration_args)
    templates=[]
    for face in(100.,1.):
      for association in('weighted_cash','weighted_present_value','time_identity'):
        g=Builder(3);discount,n,frequency=range(3);cash=g.constant(bits(face))
        term=g.op('div',cash,discount)
        if association=='weighted_cash':numerator=g.op('div',g.op('mul',n,cash),discount)
        elif association=='weighted_present_value':numerator=g.op('mul',n,term)
        else:numerator=None
        result=g.op('div',n,frequency)if numerator is None else g.op('div',g.op('div',numerator,term),frequency)
        templates.append(dict(face=face,association=association,graph=g.graph(result,result)))
    programs=[]
    for template in templates:
      graph=template['graph']
      for base in('strict','stored_pc64'):
        programs.append(dict(template={k:v for k,v in template.items()if k!='graph'},base=base,word='strict',
                             cw=0x123f,mask=(1<<graph['sites'])-1,matches=0,first_failure=None))
        for mask in range(1<<graph['sites']):
            programs.append(dict(template={k:v for k,v in template.items()if k!='graph'},base=base,word='pc64',
                                 cw=0x133f,mask=mask,matches=0,first_failure=None))
    counts=Counter();records=[]
    def integer_power(base,n,word):
        acc=Q(1);b=Q(base)
        while n:
            if n&1:acc=operation(acc,b,'mul',word)
            n>>=1
            if n:b=operation(b,b,'mul',word)
        return bits(float(acc))
    for index,(c,slots,d,md,price,base_sum)in enumerate(zip(cases,power_slots,duration,modified,prices,sums)):
        predictions={};matched=[]
        for ti,template in enumerate(templates):
            graph=template['graph']
            for base in('strict','stored_pc64'):
                operands=[int(powers[slots[base]]['bits'],16),int(bits(c['n']),16),int(bits(c['frequency']),16)]
                strict=native.evaluate(graph,operands,0x123f,(1<<graph['sites'])-1)
                masks=native.masks(graph,operands,0x133f)
                predictions[ti,base,'strict',(1<<graph['sites'])-1]=f'{native.publish(strict[:10]):016x}'
                for mask in range(1<<graph['sites']):predictions[ti,base,'pc64',mask]=f'{native.publish(masks[20*mask:20*mask+10]):016x}'
        cursor=0
        for ti,template in enumerate(templates):
          for base in('strict','stored_pc64'):
            for word,mask in [('strict',(1<<template['graph']['sites'])-1)]+[('pc64',m)for m in range(1<<template['graph']['sites'])]:
                prediction=predictions[ti,base,word,mask];program=programs[cursor]
                good=d.get('type')=='number'and d['bits']==prediction
                program['matches']+=good
                if good:matched.append(cursor)
                elif program['first_failure']is None:program['first_failure']=dict(index=index,answer=d,prediction=prediction)
                cursor+=1
        counts['rows']+=1;counts['all_declared_programs_refuted']+=not matched
        controls={}
        for word in('strict','stored_pc64'):
            discount=Q(value(powers[slots['strict']]['bits']))
            controls['price_'+word]=bits(float(operation(Q(100),discount,'div',word)))
            counts['price_'+word+'_matches']+=price.get('type')=='number'and price['bits']==controls['price_'+word]
            if d.get('type')=='number':
                for base in('strict','stored_pc64'):
                    key='mduration_'+word+'_base_'+base
                    controls[key]=bits(float(operation(Q(value(d['bits'])),Q(value(c['base_'+base])),'div',word)))
                    counts[key+'_matches']+=md.get('type')=='number'and md['bits']==controls[key]
        integer_powers={word:integer_power(value(c['base_strict']),c['n'],word)for word in('strict','stored_pc64')}
        for word,prediction in integer_powers.items():counts['integer_power_'+word+'_matches']+=powers[slots['strict']]['bits']==prediction
        records.append(dict(c,powers={base:powers[slots[base]]for base in slots},public_sum=base_sum,
                            duration=d,modified_duration=md,price=price,controls=controls,matching_programs=matched,
                            integer_power_predictions=integer_powers))
    result=dict(format='dnacalc-duration-zero-coupon/1',complete=True,counts=counts,templates=templates,programs=programs,
        complete_programs=[p for p in programs if p['matches']==len(cases)],records=records,receipts=service.receipts,
        source_bindings=bindings,library_sha256=sha(library),questions_sha256=sha(str(output)+'.questions.json'),
        limits='Conditional weighted-redemption models constrained by observed POWER and exact on-coupon schedules. '
        'Does not identify nonzero coupon sums, fractional discount exponents or a private POWER implementation. '
        'Every nonnumeric answer remains a failure; partial model matches are not a selection rule.')
    write(output,result);print(json.dumps(dict(counts=counts,complete_programs=result['complete_programs'])),flush=True)
    return result


def duration_two_coupon_transfer(checkout,output,cache,offline=False,trials=65536):
    """Nonzero coupons whose complete coupon contribution is exact dyadic arithmetic."""
    from exact import operation
    import itertools,datetime as dt
    if not 1024<=trials<=1048576:raise ValueError('bounded primitive-discriminant mining required')
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    seed=0x33d012;mask=(1<<64)-1;seen=set();discriminants=[]
    for _ in range(trials):
        seed^=(seed<<13)&mask;seed^=seed>>7;seed^=(seed<<17)&mask;seed&=mask
        numerator=5*((1<<25)//5+1+seed%((1<<25)//5-1))
        if numerator in seen:continue
        seen.add(numerator);q=Q(numerator,1<<25);square=q*q
        if BINARY64.round(square)!=square:raise ValueError('declared exact square does not fit')
        exact=Q(100)/square;strict=BINARY64.round(exact);stored=BINARY64.round(EXTENDED64.round(exact))
        if strict!=stored:discriminants.append(dict(base_bits=bits(float(q)),square_bits=bits(float(square)),
            redemption_strict=bits(float(strict)),redemption_stored_pc64=bits(float(stored))))
    cases=[]
    for item in discriminants:
      q=Q(value(item['base_bits']));square=q*q
      for frequency in(1,2,4):
        months=2020*12+2*(12//frequency);year,month=divmod(months,12)
        maturity=(dt.date(year,month+1,1)-dt.date(1899,12,30)).days
        for exponent in(-52,-51,-20,-1,0,1,20,51,52):
            scale=Q(2)**exponent;coupon=square*scale;rate=coupon*frequency/100
            if BINARY64.round(rate)!=rate or Q(100)*rate/frequency!=coupon:raise ValueError('coupon prefix is not exact')
            cases.append(dict(item,settlement=43831,maturity=maturity,frequency=frequency,
                yield_bits=bits(float((q-1)*frequency)),coupon_rate_bits=bits(float(rate)),scale_exponent=exponent,
                coupon_denominator_prefix=bits(float((q+1)*scale)),coupon_numerator_prefix=bits(float((q+2)*scale))))
    source=Path(checkout)/'crates/calc-func/src/ported/bond_core_family.rs'
    bindings={str(source.resolve()):sha(source),str(Path(__file__).resolve()):sha(__file__)}
    snapshots=Path(output).resolve().parent/'producers';snapshots.mkdir(exist_ok=True)
    for path,digest in bindings.items():
        dest=snapshots/(digest+'-'+Path(path).name)
        if not dest.exists():dest.write_bytes(Path(path).read_bytes())
    write(str(output)+'.questions.json',dict(format='dnacalc-duration-two-coupon-questions/1',complete=True,cases=cases,
        trials=trials,distinct_bases=len(seen),discriminants=discriminants,source_bindings=bindings,
        construction='q=k/2^25 with k divisible by5, so q^2 and q^2*2^e*f/100 are exact binary64. '
        'Two coupons have exact PV sum(q+1)*2^e and exact weighted sum(q+2)*2^e. '
        'Select only100/q^2 strict/PC64 division disagreements, before asking any target. '
        'Redemption quotient, final denominator sum, numerator sum and ratio each have two fixed stored arithmetic words.'))
    if not cases:raise ValueError('no primitive discriminants in declared mining scope')
    powers=service.ask_rows('POWER',[[client.number(c['base_bits']),client.number(2.)]for c in cases])
    if any(r!=dict(type='number',bits=c['square_bits'])for c,r in zip(cases,powers)):raise ValueError('exact square public POWER changed')
    arguments=[[client.number(c['settlement']),client.number(c['maturity']),client.number(c['coupon_rate_bits']),
                client.number(c['yield_bits']),client.number(c['frequency']),client.number(0.)]for c in cases]
    prices=service.ask_rows('PRICE',[a[:4]+[client.number(100.)]+a[4:]for a in arguments])
    answers=service.ask_rows('DURATION',arguments);modified=service.ask_rows('MDURATION',arguments)
    programs=[dict(words=list(words),matches=0,first_failure=None)for words in itertools.product(('strict','stored_pc64'),repeat=4)]
    price_programs=[dict(words=list(words),matches=0,first_failure=None)for words in itertools.product(('strict','stored_pc64'),repeat=2)]
    records=[];counts=Counter()
    for index,(case,answer,price,md,power)in enumerate(zip(cases,answers,prices,modified,powers)):
        den0=Q(value(case['coupon_denominator_prefix']));num0=Q(value(case['coupon_numerator_prefix']));matched=[]
        for pi,program in enumerate(programs):
            word,dw,nw,fw=program['words'];red=Q(value(case['redemption_'+word]))
            den=operation(den0,red,'add',dw);num=operation(num0,2*red,'add',nw)
            ratio=operation(num,den,'div',fw);predicted=bits(float(ratio/case['frequency']))
            good=answer.get('type')=='number'and answer['bits']==predicted;program['matches']+=good
            if good:matched.append(pi)
            elif program['first_failure']is None:program['first_failure']=dict(index=index,answer=answer,prediction=predicted)
        for program in price_programs:
            word,finish=program['words'];red=Q(value(case['redemption_'+word]))
            predicted=bits(float(operation(den0,red,'add',finish)))
            good=price.get('type')=='number'and price['bits']==predicted;program['matches']+=good
            if not good and program['first_failure']is None:program['first_failure']=dict(index=index,answer=price,prediction=predicted)
        for word in('strict','stored_pc64'):
            if answer.get('type')=='number':
                predicted=bits(float(operation(Q(value(answer['bits'])),Q(value(case['base_bits'])),'div',word)))
                counts['mduration_'+word+'_matches']+=md.get('type')=='number'and md['bits']==predicted
        counts['rows']+=1;counts['all_programs_refuted']+=not matched
        records.append(dict(case,answer=answer,price=price,modified_duration=md,power=power,matching_programs=matched))
    result=dict(format='dnacalc-duration-two-coupon-transfer/1',complete=True,counts=counts,programs=programs,
        complete_programs=[p for p in programs if p['matches']==len(cases)],price_programs=price_programs,
        complete_price_programs=[p for p in price_programs if p['matches']==len(cases)],records=records,
        receipts=service.receipts,source_bindings=bindings,questions_sha256=sha(str(output)+'.questions.json'),
        limits='Exact controlled nonzero coupon contributions on two-period on-coupon bonds only. '
        'No arbitrary coupon sum, fractional POWER, unobserved helper or complete duration kernel is inferred.')
    write(output,result);print(json.dumps(dict(counts=counts,complete=result['complete_programs'],price=result['complete_price_programs'])),flush=True)
    return result


def duration_off_coupon(checkout,output,cache,library,offline=False,count=16,direct_days=False):
    """Test coupon-term words with independently observed fractional powers."""
    import itertools,datetime as dt
    if not 8<=count<=256:raise ValueError('bounded base window required')
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    native=Native(library);states=[];cases=[];ladders=[];power_requests=[]
    for frequency in(1,2,4):
      for n in(2,3,5,9,17,31):
        months=2020*12+n*(12//frequency);year,month=divmod(months,12)
        maturity=(dt.date(year,month+1,1)-dt.date(1899,12,30)).days
        e=360./frequency;off=(e-41.)/e
        state=dict(settlement=43873,maturity=maturity,frequency=frequency,n=n,e=e,a=41.,off_bits=bits(off));states.append(state)
        for center in(1.001,1.03125,1.125,1.5,1.875):
          anchor=int(bits(center),16)
          for step in range(-count//2,count-count//2):
            q=value(f'{anchor+step:016x}');yld=(q-1.)*frequency
            if 1.+yld/frequency!=q:raise ValueError('exact base construction failed')
            slots=[]
            for k in range(n):
                exponent=(e-41.+k*e)/e if direct_days else off+k
                if exponent<=0 or exponent==.5 or exponent.is_integer():raise ValueError('fractional-power scope changed')
                slots.append(len(power_requests));power_requests.append([client.number(q),client.number(exponent)])
            ladder_index=len(ladders);ladders.append(dict(state,base_bits=bits(q),yield_bits=bits(yld),slots=slots))
            for rate in(.03,.06,.123456789):cases.append(dict(state,base_bits=bits(q),yield_bits=bits(yld),coupon_rate_bits=bits(rate),ladder=ladder_index))
    source=Path(checkout)/'crates/calc-func/src/ported/bond_core_family.rs'
    bindings={str(source.resolve()):sha(source),str(Path(__file__).resolve()):sha(__file__),
              str(ROOT/'racer/x87_programs.py'):sha(ROOT/'racer/x87_programs.py'),str(ROOT/'racer/x87_native.c'):sha(ROOT/'racer/x87_native.c')}
    snapshots=Path(output).resolve().parent/'producers';snapshots.mkdir(exist_ok=True)
    for path,digest in bindings.items():
        dest=snapshots/(digest+'-'+Path(path).name)
        if not dest.exists():dest.write_bytes(Path(path).read_bytes())
    word_profiles=list(itertools.product(('strict','stored_pc64'),repeat=5))
    write(str(output)+'.questions.json',dict(format='dnacalc-duration-off-coupon-questions/1',complete=True,cases=cases,
        states=states,ladders=ladders,source_bindings=bindings,word_profiles=word_profiles,direct_days=direct_days,
        sites=['weighted coupon product','coupon quotient','numerator sum','denominator sum','final ratio'],
        selection='Fixed18 schedules, five adjacent-double base windows and three coupons. Positive noninteger '
        'exponents avoid the public square-root/integer special cases. Every discount POWER is observed separately. '
        'Redemption divisions are fixed to the independently proven stored-PC64 word; other words are independent. '
        'The optional direct-day time is RN53((E-A+k*E)/E), separately passed to POWER and the weighted product.'))
    factor_args=[[client.number(s[k])for k in('settlement','maturity','frequency')]+[client.number(0.)]for s in states]
    for name in('COUPNUM','COUPDAYS','COUPDAYBS'):
        readings=service.ask(name,factor_args)['rows']
        field={'COUPNUM':'n','COUPDAYS':'e','COUPDAYBS':'a'}[name]
        if any(r!=dict(type='number',bits=bits(s[field]))for s,r in zip(states,readings)):raise ValueError('off-coupon public factors changed')
    powers=service.ask_rows('POWER',power_requests)
    if any(r['type']!='number'for r in powers):raise ValueError('nonnumeric fractional discount control')
    args=[[client.number(c['settlement']),client.number(c['maturity']),client.number(c['coupon_rate_bits']),
           client.number(c['yield_bits']),client.number(c['frequency']),client.number(0.)]for c in cases]
    answers=service.ask_rows('DURATION',args);modified=service.ask_rows('MDURATION',args)
    # Continue only at explicitly stored binary64 accumulators. This keeps
    # every graph within the existing128-node interpreter bound without
    # changing that bound or introducing a store into a resident program.
    graphs={};finishes={}
    for pi,words in enumerate(word_profiles):
      wm,qd,na,da,fd=words
      for width in(1,2,3,5,7,8):
       for last in(False,True):
        g=Builder(3+2*width);hundred=g.constant(bits(100.))
        def op(kind,a,b,word):
            at=g.op(kind,a,b,store=False);g.nodes[at]['site']=0
            if word=='strict':g.nodes[at]['phase']=1
            return at
        coupon=0;num=1;den=2
        for j in range(width):
            disc=3+2*j;diff=4+2*j
            num=op('add',num,op('div',op('mul',diff,coupon,wm),disc,qd),na)
            den=op('add',den,op('div',coupon,disc,qd),da)
            if last and j+1==width:
                num=op('add',num,op('div',op('mul',diff,hundred,'strict'),disc,'stored_pc64'),na)
                den=op('add',den,op('div',hundred,disc,'stored_pc64'),da)
        g.sites=1;graphs[pi,width,last]=g.graph(num,den)
      g=Builder(3)
      quotient=g.op('div',0,1,store=False);g.nodes[quotient].update(site=0,phase=int(fd=='strict'))
      result=g.op('div',quotient,2,store=False);g.nodes[result].update(site=0,phase=1)
      g.sites=1;finishes[pi]=g.graph(result,result)
    programs=[dict(words=list(w),matches=0,first_failure=None)for w in word_profiles];records=[];counts=Counter()
    for index,(case,answer,md)in enumerate(zip(cases,answers,modified)):
        ladder=ladders[case['ladder']];discounts=[powers[i]for i in ladder['slots']]
        matched=[]
        for pi,program in enumerate(programs):
            coupon=int(bits(100.*value(case['coupon_rate_bits'])/case['frequency']),16);num=den=0
            for start in range(0,case['n'],8):
                width=min(8,case['n']-start);last=start+width==case['n'];operands=[coupon,num,den]
                for k in range(start,start+width):
                    diff=(case['e']-case['a']+k*case['e'])/case['e']if direct_days else value(case['off_bits'])+k
                    operands.extend((int(discounts[k]['bits'],16),int(bits(diff),16)))
                outputs=native.evaluate(graphs[pi,width,last],operands,0x133f,1,inner_cw=0x123f)
                num,den=native.publish(outputs[:10]),native.publish(outputs[10:])
            outputs=native.evaluate(finishes[pi],[num,den,int(bits(case['frequency']),16)],0x133f,1,inner_cw=0x123f)
            prediction=f'{native.publish(outputs[:10]):016x}'
            good=answer.get('type')=='number'and answer['bits']==prediction;program['matches']+=good
            if good:matched.append(pi)
            elif program['first_failure']is None:program['first_failure']=dict(index=index,answer=answer,prediction=prediction)
        counts['rows']+=1;counts['all_declared_programs_refuted']+=not matched
        records.append(dict(case,discounts=discounts,answer=answer,modified_duration=md,matching_programs=matched))
    result=dict(format='dnacalc-duration-off-coupon/1',complete=True,counts=counts,programs=programs,
        complete_programs=[p for p in programs if p['matches']==len(cases)],records=records,receipts=service.receipts,
        source_bindings=bindings,library_sha256=sha(library),questions_sha256=sha(str(output)+'.questions.json'),direct_days=direct_days,
        limits='Named binary64 coupon and observed positive fractional POWER leaves are declared premises. '
        'PC53 and stored PC64 words vary independently at five sites; all operands at these sites are stored '
        'binary64. Does not identify private helpers or nonpositive/half-integer discount branches.')
    write(output,result);print(json.dumps(dict(counts=counts,complete=result['complete_programs'])),flush=True)
    return result


def sinh_small_stages(checkout,output,cache,library,offline=False,bank_trace=None,count=128):
    """Constrain the paired Kahan terms with public EXP/LN leaves."""
    from exact import operation
    if not 8<=count<=1024:raise ValueError('bounded small-SINH window required')
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline);native=Native(library)
    inputs={}
    for center in[math.ldexp(1.,e)for e in(-60,-54,-53,-52,-40,-27,-26,-20,-10)]+[.01,.125,.5,.9]:
        anchor=int(bits(center),16)
        for step in range(-count//2,count-count//2):
            for sign in(1.,-1.):inputs.setdefault(bits(sign*value(f'{anchor+step:016x}')),[]).append('fixed adjacent/exponent window')
    if bank_trace:
        for line in Path(bank_trace).read_text().splitlines():
            row=json.loads(line)
            if row['function']=='SINH'and row['actual']!=row['expected']:
                raw=row['args'][0].removeprefix('0x');x=value(raw)
                if 0<abs(x)<1:
                    inputs.setdefault(raw,[]).append('original bank failure input')
                    inputs.setdefault(bits(-x),[]).append('paired sign of original bank failure input')
    raw_inputs=list(inputs)
    source=Path(checkout)/'crates/calc-func/src/ported/sinh.rs'
    bindings={str(source.resolve()):sha(source),str(Path(__file__).resolve()):sha(__file__),
              str(ROOT/'racer/x87_programs.py'):sha(ROOT/'racer/x87_programs.py'),str(ROOT/'racer/x87_native.c'):sha(ROOT/'racer/x87_native.c')}
    snapshots=Path(output).resolve().parent/'producers';snapshots.mkdir(exist_ok=True)
    for path,digest in bindings.items():
        dest=snapshots/(digest+'-'+Path(path).name)
        if not dest.exists():dest.write_bytes(Path(path).read_bytes())
    write(str(output)+'.questions.json',dict(format='dnacalc-sinh-small-questions/1',complete=True,
        inputs=[dict(bits=raw,roles=inputs[raw])for raw in raw_inputs],source_bindings=bindings,
        bank_trace_sha256=sha(bank_trace)if bank_trace else None,
        selection='Fixed exponent/adjacent windows plus explicitly labeled old failure inputs and their paired signs. '
        'EXP(+x),EXP(-x) and LN of each published EXP are independent observations. Kahan numerator/division '
        'stores and the outer subtraction are varied without changing any logarithm producer or coefficient.'))
    exp_rows=service.ask_rows('EXP',[[client.number(raw)]for raw in raw_inputs]+[[client.number(-value(raw))]for raw in raw_inputs])
    if any(r['type']!='number'for r in exp_rows):raise ValueError('nonnumeric small exponential control')
    log_rows=service.ask_rows('LN',[[client.number(r['bits'])]for r in exp_rows])
    if any(r['type']!='number'for r in log_rows):raise ValueError('nonnumeric small logarithm control')
    args=[[client.number(raw)]for raw in raw_inputs]
    answers=service.ask_rows('SINH',args);coshes=service.ask_rows('COSH',args)
    cosecants=service.ask_rows('CSCH',args);tangents=service.ask_rows('TANH',args)
    graphs={};programs=[]
    for association in('product_first','ratio_first','input_ratio_first'):
      for identities in((False,False),(False,True),(True,False),(True,True)):
        g=Builder(5);negative=g.op('neg',0,store=False)
        def op(kind,a,b,site):
            at=g.op(kind,a,b,store=False);g.nodes[at]['site']=site;return at
        def term(t,u,log,identity):
            if identity:return t
            difference=op('sub',u,g.one,0)
            if association=='product_first':return op('div',op('mul',difference,t,1),log,2)
            if association=='ratio_first':return op('mul',op('div',difference,log,2),t,1)
            return op('mul',difference,op('div',t,log,2),1)
        plus=term(0,1,3,identities[0]);minus=term(negative,2,4,identities[1])
        result=op('div',op('sub',plus,minus,3),g.constant(bits(2.)),4);g.sites=5
        graphs[association,identities]=g.graph(result,result)
      for cw in(0x123f,0x133f):
       for mask in range(32):programs.append(dict(association=association,cw=cw,mask=mask,matches=0,first_failure=None))
    records=[];counts=Counter();size=len(raw_inputs)
    for i,(raw,answer,cosh,csch,tanh)in enumerate(zip(raw_inputs,answers,coshes,cosecants,tangents)):
        plus,minus=exp_rows[i],exp_rows[size+i];plog,mlog=log_rows[i],log_rows[size+i]
        identities=(value(plus['bits'])==1.,value(minus['bits'])==1.)
        operands=list(map(lambda x:int(x,16),(raw,plus['bits'],minus['bits'],plog['bits'],mlog['bits'])));predictions={}
        for association in('product_first','ratio_first','input_ratio_first'):
          for cw in(0x123f,0x133f):
            values=native.masks(graphs[association,identities],operands,cw)
            for mask in range(32):predictions[association,cw,mask]=f'{native.publish(values[20*mask:20*mask+10]):016x}'
        matched=[]
        for pi,program in enumerate(programs):
            prediction=predictions[program['association'],program['cw'],program['mask']]
            good=answer.get('type')=='number'and answer['bits']==prediction;program['matches']+=good
            if good:matched.append(pi)
            elif program['first_failure']is None:program['first_failure']=dict(index=i,answer=answer,prediction=prediction)
        controls={}
        for word in('strict','stored_pc64'):
            if answer.get('type')=='number':
                s=Q(value(answer['bits']));controls['csch_'+word]=bits(float(operation(Q(1),s,'div',word)))
                counts['csch_'+word+'_matches']+=csch.get('type')=='number'and csch['bits']==controls['csch_'+word]
                if cosh.get('type')=='number':
                    controls['tanh_'+word]=bits(float(operation(s,Q(value(cosh['bits'])),'div',word)))
                    counts['tanh_'+word+'_matches']+=tanh.get('type')=='number'and tanh['bits']==controls['tanh_'+word]
        counts['rows']+=1;counts['all_declared_programs_refuted']+=not matched
        counts['identity_matches']+=answer.get('type')=='number'and answer['bits']==raw
        records.append(dict(input=raw,roles=inputs[raw],exp_plus=plus,exp_minus=minus,ln_plus=plog,ln_minus=mlog,
                            answer=answer,cosh=cosh,csch=csch,tanh=tanh,controls=controls,matching_programs=matched))
    result=dict(format='dnacalc-sinh-small-stages/1',complete=True,counts=counts,records=records,programs=programs,
        complete_programs=[p for p in programs if p['matches']==size],receipts=service.receipts,source_bindings=bindings,
        library_sha256=sha(library),questions_sha256=sha(str(output)+'.questions.json'),
        limits='Public EXP/LN leaves are fixed, not fitted or replaced by a private alnrel producer. '
        'Arithmetic masks identify only declared functional compositions. No tiny-input threshold is inferred from '
        'partial matches; original failures and fixed windows remain separately labeled.')
    write(output,result);print(json.dumps(dict(counts=counts,complete=result['complete_programs'])),flush=True)
    return result


def weeknum_controls(checkout,output,cache,bank_trace,offline=False,repeat_count=1):
    """Check unsupported return-type readings against new exact scalar calls."""
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    source=Path(bank_trace);rows=[json.loads(line)for line in source.read_text().splitlines()]
    requests={};old=[]
    def add(arguments,role):
        key=tuple(arguments);requests.setdefault(key,[]).append(role)
    for serial in(0.,1.,2.,59.,60.,61.,43831.,44013.,2958465.):
        for kind in range(65):add((bits(serial),bits(kind)),'fixed serial/return-type grid')
    for kind in(1.,2.,11.,17.,21.,36.,54.,float(1<<31),float(1<<32)+21.):
        for value_ in(math.nextafter(kind,-math.inf),kind,math.nextafter(kind,math.inf)):
            add((bits(43831.),bits(value_)),'adjacent return-type boundary')
    for row in rows:
        if row['function']=='WEEKNUM'and row['actual']!=row['expected']:
            if all(isinstance(a,str)and re.fullmatch('0x[0-9a-f]{16}',a)for a in row['args']):
                key=tuple(a.removeprefix('0x')for a in row['args']);add(key,'original mismatched input');old.append(row)
    if not 1<=repeat_count<=8:raise ValueError('bounded repeat count required')
    keys=list(requests)*repeat_count
    if repeat_count>1:old=[]
    write(str(output)+'.questions.json',dict(format='dnacalc-weeknum-controls-questions/1',complete=True,
        arguments=[dict(args=list(k),roles=requests[k])for k in keys],bank_trace_sha256=sha(source),repeat_count=repeat_count,
        selection='Fixed serial and return-type grid, adjacent return-type edges, and original mismatched numeric inputs. '
        'Expected answers are absent from service questions. Repeated historical input conflicts are retained, not resolved by a row selector.'))
    answers=service.ask_rows('WEEKNUM',[[client.number(x)for x in key]for key in keys]);by_key=dict(zip(keys,answers));counts=Counter()
    def rendered(answer):
        if answer['type']=='number':return'0x'+answer['bits']
        if answer['type']=='error':return'error:'+{'#NUM!':'Num','#VALUE!':'Value','#DIV/0!':'Div0','#N/A':'NA'}.get(answer['code'],answer['code'])
        return answer['type']
    historical=[]
    for row in old:
        key=tuple(a.removeprefix('0x')for a in row['args']);fresh=by_key[key];counts['original_rows']+=1
        counts['original_expected_repeated_live']+=rendered(fresh)==row['expected']
        counts['current_calculation_repeated_live']+=rendered(fresh)==row['actual']
        historical.append(dict(original=row,fresh=fresh))
    unsupported=[]
    for key,answer in zip(keys,answers):
        kind=math.trunc(value(key[1]));counts['calls']+=1
        if kind not in(1,2,11,12,13,14,15,16,17,21)and answer['type']=='number':unsupported.append(dict(args=key,answer=answer))
    repeated={}
    for key,answer in zip(keys,answers):repeated.setdefault(key,set()).add(rendered(answer))
    counts['unique_inputs']=len(repeated);counts['inputs_with_changed_answers']=sum(len(v)>1 for v in repeated.values())
    result=dict(format='dnacalc-weeknum-controls/1',complete=True,counts=counts,historical=historical,
        repeat_count=repeat_count,repeated_conflicts=[dict(args=key,answers=sorted(values))for key,values in repeated.items()if len(values)>1],
        unsupported_numeric_controls=unsupported,records=[dict(args=key,answer=answer,roles=requests[key])for key,answer in zip(keys,answers)],
        receipts=service.receipts,bank_trace_sha256=sha(source),questions_sha256=sha(str(output)+'.questions.json'),
        limits='Fresh scalar inputs and original readings only. A historical contradiction is not repaired by changing the evaluator '
        'or discarding an original failure. Visible input metadata does not establish every effective callee input.')
    write(output,result);print(json.dumps(dict(counts=counts,unsupported_numeric_controls=len(unsupported))),flush=True);return result


def yield_multi_coupon(checkout, output, cache, offline=False):
    """On-coupon zero-coupon inversions and paired public PRICE neighborhoods."""
    from datetime import date
    import itertools
    client, _ = load_service(checkout)
    service = CachedService(client.Client(timeout=300), cache, offline)
    output = Path(output)
    output.mkdir(exist_ok=False)
    settlement = date(2024, 1, 15)
    serial = lambda d: float((d-date(1899, 12, 30)).days)
    cases = []
    for count, frequency, redemption, ratio in itertools.product(
            (2, 3, 4, 8, 16, 32), (1, 2, 4), (100., 128.), (.5, .75, 1., 1.25, 2.)):
        months = count*12//frequency
        index = settlement.year*12+settlement.month-1+months
        maturity = date(index//12, index%12+1, settlement.day)
        arguments = [serial(settlement), serial(maturity), 0., redemption*ratio,
                     redemption, float(frequency), 0.]
        cases.append(dict(id=str(len(cases)), count=count, frequency=frequency,
                          redemption=redemption, price=redemption*ratio,
                          arguments=[client.number(v) for v in arguments]))
    write(output/'questions.json', cases)
    yields = service.ask_rows('YIELD', [c['arguments'] for c in cases])
    prices, mapping, untransportable = [], [], []
    for index, (case, answer) in enumerate(zip(cases, yields)):
        if answer['type'] != 'number': continue
        y = value(answer['bits'])
        for offset in range(-4, 5):
            neighbor = y
            for _ in range(abs(offset)):
                neighbor = math.nextafter(neighbor, -math.inf if offset<0 else math.inf)
            if neighbor != 0. and abs(neighbor) < sys.float_info.min:
                untransportable.append(dict(case=index, offset=offset, yield_bits=bits(neighbor),
                    reason='Subnormal neighbor is changed by the worksheet input transport; no reading requested.'))
                continue
            args = list(case['arguments'])
            args[3] = client.number(neighbor)
            prices.append(args)
            mapping.append(dict(case=index, offset=offset, yield_bits=bits(neighbor)))
    price_answers = service.ask_rows('PRICE', prices, chunk_size=2000)
    powers = service.ask_rows('POWER', [[client.number(c['redemption']/c['price']),
        client.number(1./c['count'])] for c in cases])
    observations = dict(yields=yields, price_mapping=mapping, prices=price_answers, powers=powers,
                        untransportable_neighbors=untransportable)
    write(output/'answers.json', observations)
    counts = Counter()
    for c, y, p in zip(cases, yields, powers):
        counts['yield_'+y['type']] += 1
        if y['type'] == p['type'] == 'number':
            analytic = c['frequency']*(value(p['bits'])-1.)
            counts['public_power_composition_equal'] += bits(analytic) == y['bits']
    for mapped, answer in zip(mapping, price_answers):
        if mapped['offset'] == 0:
            c = cases[mapped['case']]
            counts['public_price_at_returned_yield_equal'] += answer.get('bits') == bits(c['price'])
    result = dict(complete=True, cases=len(cases), calls=len(cases)*2+len(prices), counts=counts,
        untransportable_neighbors=len(untransportable),
        receipts=service.receipts, source_sha256=sha(Path(__file__)),
        limits='Controlled on-coupon zero-coupon bonds; public PRICE can differ from a solver private leaf. '
               'Paired PRICE and POWER calls are diagnostics, never substituted Excel answers or an identified solver.')
    write(output/'summary.json', result)
    print(json.dumps(dict(cases=len(cases), counts=counts)))


def yield_domain_raster(checkout,output,cache,offline=False):
    """Read the negative-root numerical domain at exact dyadic price ratios."""
    from datetime import date
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    serial=lambda d:float((d-date(1899,12,30)).days);cases=[]
    for frequency in(1,2,4):
        month=1+24//frequency;year=2024+(month-1)//12;month=(month-1)%12+1
        for scale in(-20,0,20):
            price=math.ldexp(128.,scale)
            for numerator in range(1,128):
                cases.append(dict(frequency=frequency,scale=scale,numerator=numerator,
                    args=[client.number(x)for x in(serial(date(2024,1,15)),serial(date(year,month,15)),0.,price,math.ldexp(float(numerator),scale),float(frequency),0.)]))
    write(str(output)+'.questions.json',dict(cases=cases,source_sha256=sha(Path(__file__)),
        selection='All127 fixed dyadic redemption/price ratios below par, three exact money scales and three frequencies; two coupons on-calendar. All mathematical roots are inside the legal negative-yield interval.'))
    answers=service.ask_rows('YIELD',[c['args']for c in cases]);groups=[]
    for frequency in(1,2,4):
        for scale in(-20,0,20):
            rows=[(c,a)for c,a in zip(cases,answers)if c['frequency']==frequency and c['scale']==scale]
            runs=[]
            for c,a in rows:
                kind=a['type']if a['type']=='number'else a['code']
                if runs and runs[-1]['kind']==kind:runs[-1]['last']=c['numerator']
                else:runs.append(dict(kind=kind,first=c['numerator'],last=c['numerator']))
            groups.append(dict(frequency=frequency,scale=scale,runs=runs))
    write(output,dict(complete=True,cases=cases,answers=answers,groups=groups,receipts=service.receipts,
        source_sha256=sha(Path(__file__)),limits='Numerical admission raster only. NUM can arise at any iteration; an edge alone does not identify an initial guess, derivative or fallback. No tolerance or seed fitted.'))
    print(json.dumps(groups),flush=True)


def yield_exact_roots(checkout, output, cache, offline=False):
    """Dyadic zero-coupon roots certify prices before solver observations."""
    from datetime import date
    import itertools
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    serial=lambda d:float((d-date(1899,12,30)).days);settlement=date(2024,1,15);cases=[]
    for n,f,scale,q in itertools.product((2,3,4,8),(1,2,4),(-20,0,20),(.5,.75,1.,1.25,1.5,2.,4.)):
        months=n*12//f;index=settlement.year*12+settlement.month-1+months;maturity=date(index//12,index%12+1,15)
        price=math.ldexp(128.,scale);red=price*q**n
        if Q(red)!=Q(price)*Q(q)**n:raise ValueError('uncertified dyadic price')
        args=[serial(settlement),serial(maturity),0.,price,red,float(f),0.]
        first_index=settlement.year*12+settlement.month-1+12//f;first=date(first_index//12,first_index%12+1,15)
        odd=[args[0],args[1],serial(date(2023,1,15)),serial(first),0.,price,red,float(f),0.]
        cases.append(dict(n=n,frequency=f,scale=scale,q=q,exact_root=bits(float(f*(q-1))),args=[client.number(x)for x in args],odd=[client.number(x)for x in odd]))
    write(str(output)+'.questions.json',dict(cases=cases,source_sha256=sha(Path(__file__)),
        selection='Fixed dyadic base, coupon count, frequency and price powers; exact rational equality certifies each price before asking either solver.'))
    answers={f:service.ask_rows(f,[c['args'if f=='YIELD'else'odd']for c in cases])for f in ('YIELD','ODDFYIELD')}
    counts=Counter()
    for i,c in enumerate(cases):
        for f in answers:
            a=answers[f][i];counts[f+'/'+a['type']]+=1;counts[f+'/exact_root']+=a==dict(type='number',bits=c['exact_root'])
        counts['paired_equal']+=answers['YIELD'][i]==answers['ODDFYIELD'][i]
    write(output,dict(complete=True,cases=cases,answers=answers,counts=counts,receipts=service.receipts,
        source_sha256=sha(Path(__file__)),limits='Exact mathematical roots and paired public solver observations only. A close numerical return is not scored as exact. No private iteration, initial value or tolerance identified.'))
    print(json.dumps(counts),flush=True)


def yield_date_domain(checkout, output, cache, offline=False):
    """Numeric date-error readings through the same scalar question service."""
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    cases=[];base=[39493.,42689.,.0575,95.04287,100.,2.,0.]
    for slot in (0,1):
        for date in (-2147483648.,-1.,0.,.5,1.,1.5,59.,59.5,60.,60.5,61.,
                     2958465.,2958465.5,2958466.,4294967295.,4294967296.):
            args=base.copy();args[slot]=date
            cases.append(dict(slot=slot,date=date,args=[client.number(x)for x in args]))
    write(str(output)+'.questions.json',dict(cases=cases,source_sha256=sha(Path(__file__)),
        selection='Fixed numeric date boundaries and halves; no subnormal input whose cell transport changes it.'))
    answers=service.ask_rows('YIELD',[r['args']for r in cases])
    write(output,dict(complete=True,cases=cases,answers=answers,receipts=service.receipts,
        source_sha256=sha(Path(__file__)),limits='Date admission only; numeric solver outputs are retained observations and do not identify a solver.'))
    print(json.dumps(Counter(json.dumps(a,sort_keys=True)if a['type']!='number'else'number'for a in answers)))


def yield_scale_questions(checkout, output, cache, offline=False):
    """Power-of-two price/redemption scales isolate solver stopping premises."""
    from datetime import date
    import itertools
    client,_=load_service(checkout);service=CachedService(client.Client(timeout=300),cache,offline)
    serial=lambda d:float((d-date(1899,12,30)).days)
    rows=[]
    for count,scale,ratio in itertools.product((2,3,8),(-40,-10,0,10,40),(.5,.75,1.,1.25,2.)):
        redemption=math.ldexp(128.,scale)
        rows.append(dict(count=count,scale=scale,ratio=ratio,args=[client.number(x) for x in
            (serial(date(2024,1,15)),serial(date(2024+count,1,15)),0.,redemption*ratio,redemption,1.,0.)]))
    answers=service.ask_rows('YIELD',[r['args']for r in rows]);groups={}
    for r,a in zip(rows,answers):groups.setdefault((r['count'],r['ratio']),[]).append(dict(scale=r['scale'],answer=a))
    write(output,dict(complete=True,cases=rows,answers=answers,groups=[dict(count=k[0],ratio=k[1],scaled=v,
        distinct_answers=len({json.dumps(x['answer'],sort_keys=True)for x in v}))for k,v in groups.items()],
        receipts=service.receipts,source_sha256=sha(Path(__file__)),
        limits='Dyadic scaling is an exact necessary invariance of homogeneous zero-coupon price, relative stopping and coherent operations away from exponent edges. Different roots refute that combined premise, not a uniquely identified stopping rule.'))
    print(json.dumps(Counter(len({json.dumps(x['answer'],sort_keys=True)for x in v})for v in groups.values())))


def yield_solver_mapping(checkout, output, cache, offline=False):
    """Map around par, by money scale and adjacent prices before solver search."""
    from datetime import date
    client,_=load_service(checkout); service=CachedService(client.Client(timeout=180),cache,offline)
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    serial=lambda d:float((d-date(1899,12,30)).days)
    settlement=date(2024,1,15);cases=[]
    for n in (2,4,8):
        for frequency in (1,2,4):
            index=settlement.year*12+settlement.month-1+n*12//frequency
            maturity=date(index//12,index%12+1,15)
            first_index=settlement.year*12+settlement.month-1+12//frequency
            first=date(first_index//12,first_index%12+1,15)
            for coupon in (0.,.03125,.0625,.125):
                prices=(50.,75.,90.,95.,99.,99.9,100.,100.1,101.,105.,110.,125.,150.,200.)
                if n==4 and frequency==2 and coupon==.0625:
                    prices=tuple(sorted(set(prices+(tuple(100.+sign*2.**e for sign in (-1,1) for e in range(-24,5))))))
                for price in prices:
                    args=[serial(settlement),serial(maturity),coupon,price,100.,float(frequency),0.]
                    odd=[args[0],args[1],serial(date(2023,1,15)),serial(first),coupon,price,100.,float(frequency),0.]
                    cases.append(dict(id=str(len(cases)),n=n,frequency=frequency,coupon=coupon,price=price,
                        args=[client.number(x)for x in args],odd=[client.number(x)for x in odd]))
    write(output/'questions.json',dict(cases=cases,source_sha256=sha(Path(__file__)),
        selection='Fixed on-coupon mapping across coupon, frequency and count; one power-of-two ladder around par. No prediction/match-guided selection.'))
    answers={f:service.ask_rows(f,[c['args'if f=='YIELD'else'odd']for c in cases])for f in ('YIELD','ODDFYIELD')}
    rows=[]
    for i,c in enumerate(cases):
        a=answers['YIELD'][i]
        expected='0x'+a['bits']if a['type']=='number'else'error:Num'if a.get('code')=='#NUM!'else None
        if expected is None:raise ValueError('unclassified yield result')
        rows.append(dict(id=c['id'],function='YIELD',args=['0x'+x['bits']for x in c['args']],expected=expected,
            original_case=c,odd_answer=answers['ODDFYIELD'][i]))
    with(output/'solver-inputs.jsonl').open('x')as stream:
        for r in rows:stream.write(json.dumps(r,separators=(',',':'))+'\n')
    result=dict(format='dnacalc-yield-solver-mapping/1',complete=True,cases=cases,answers=answers,
        paired_equal=sum(a==b for a,b in zip(*answers.values())),receipts=service.receipts,
        source_sha256=sha(Path(__file__)),limits='Original public answers and readbacks; no solver, stop rule or iteration count identified.')
    write(output/'result.json',result);print(json.dumps(dict(cases=len(cases),paired_equal=result['paired_equal'])),flush=True)


if __name__ == '__main__':
    if len(sys.argv)>1 and sys.argv[1]=='yield-solver-map':
        yield_solver_mapping(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='yield-domain-raster':
        yield_domain_raster(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='yield-exact-roots':
        yield_exact_roots(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='yield-dates':
        yield_date_domain(*sys.argv[2:]);raise SystemExit(0)
    if len(sys.argv)>1 and sys.argv[1]=='yield-scales':
        yield_scale_questions(*sys.argv[2:]);raise SystemExit(0)


def render(result, path):
    lines = ['# Numerical fingerprint map', '',
        'Regenerate from DnaCalc-Research, using retained service replies (or asking Excel for missing replies):', '',
        '```sh',
        'python3 racer/fingerprint.py map ../DnaCalc ../DnaCalc/docs/generated/numeric-fingerprints --cache evidence/round26/fingerprints',
        '```', '',
        'Add `--offline` to require existing answers. The JSON retains exact windows, spectra, request identities and worker provenance.', '',
        *[s+'\n' for s in result['limitations']]]
    isolated=[r for r in result['functions']if 'isolated_stage'in r]
    if isolated:
        lines.extend(['R28 isolates the public LN/EXP reduction constants within the declared native instruction models, '
            'and GEOMEAN\'s exact integer exponent normalization followed by stored LN, reciprocal and EXP-argument stages. '
            'The sampled normalized composition is identified; private implementation internals and the product fallback remain open. '
            'The JSON binds the independent proofs and retains each stage\'s scope.', '',
            '| Function | Isolated stage |', '|---|---|'])
        for row in isolated:
            stage=row['isolated_stage']
            label=stage.get('model',stage.get('normalization','Conditional return/seed constraints; no algorithm adopted'))
            if row['function']=='GEOMEAN':label+='; '+stage['bounded_core']
            lines.append('| '+row['function']+' | '+label+' |')
        lines.append('')
    lines.extend([
        '| Function | Status | Rounding | Stored stages | Window peaks | Route edges |',
        '|---|---|---|---|---|---|'])
    for row in result['functions']:
        peaks = ', '.join(str(len(w['spectrum']['lines'])) for w in row['windows']) or '—'
        lines.append('| '+' | '.join([row['function'], row['status'], row['rounding'],
            str(row['stored_stages']) if row['stored_stages'] is not None else 'unresolved', peaks,
            str(len(row['band_edges'])) if row['band_edges'] else 'unmeasured'])+' |')
    Path(path).write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    multi_yield_parser = commands.add_parser('yield-multi-coupon')
    multi_yield_parser.add_argument('checkout')
    multi_yield_parser.add_argument('output')
    multi_yield_parser.add_argument('--cache', required=True)
    multi_yield_parser.add_argument('--offline', action='store_true')
    duration_parser=commands.add_parser('duration-zero-coupon')
    duration_parser.add_argument('checkout')
    duration_parser.add_argument('output')
    duration_parser.add_argument('--cache',required=True)
    duration_parser.add_argument('--library',required=True)
    duration_parser.add_argument('--count',type=int,default=64)
    duration_parser.add_argument('--holdout',action='store_true')
    duration_parser.add_argument('--offline',action='store_true')
    duration_transfer_parser=commands.add_parser('duration-two-coupon-transfer')
    duration_transfer_parser.add_argument('checkout')
    duration_transfer_parser.add_argument('output')
    duration_transfer_parser.add_argument('--cache',required=True)
    duration_transfer_parser.add_argument('--offline',action='store_true')
    duration_transfer_parser.add_argument('--trials',type=int,default=65536)
    duration_off_parser=commands.add_parser('duration-off-coupon')
    duration_off_parser.add_argument('checkout')
    duration_off_parser.add_argument('output')
    duration_off_parser.add_argument('--cache',required=True)
    duration_off_parser.add_argument('--library',required=True)
    duration_off_parser.add_argument('--offline',action='store_true')
    duration_off_parser.add_argument('--count',type=int,default=16)
    duration_off_parser.add_argument('--direct-days',action='store_true')
    sinh_parser=commands.add_parser('sinh-small-stages')
    sinh_parser.add_argument('checkout')
    sinh_parser.add_argument('output')
    sinh_parser.add_argument('--cache',required=True)
    sinh_parser.add_argument('--library',required=True)
    sinh_parser.add_argument('--offline',action='store_true')
    sinh_parser.add_argument('--bank-trace')
    sinh_parser.add_argument('--count',type=int,default=128)
    weeknum_parser=commands.add_parser('weeknum-controls')
    weeknum_parser.add_argument('checkout')
    weeknum_parser.add_argument('output')
    weeknum_parser.add_argument('--cache',required=True)
    weeknum_parser.add_argument('--bank-trace',required=True)
    weeknum_parser.add_argument('--offline',action='store_true')
    weeknum_parser.add_argument('--repeat',type=int,default=1)
    yield_parser=commands.add_parser('yield-single-coupon')
    yield_parser.add_argument('checkout')
    yield_parser.add_argument('output')
    yield_parser.add_argument('--cache',required=True)
    yield_parser.add_argument('--offline',action='store_true')
    yield_parser.add_argument('--month-ends',action='store_true')
    yield_parser.add_argument('--date-holdout',action='store_true')
    yield_programs_parser=commands.add_parser('yield-coupon-programs')
    yield_programs_parser.add_argument('original')
    yield_programs_parser.add_argument('output')
    yield_programs_parser.add_argument('--library',required=True)
    yield_programs_parser.add_argument('--admitted')
    yield_programs_parser.add_argument('--day-boundaries',action='store_true')
    yield_programs_parser.add_argument('--nasd-days',action='store_true')
    yield_programs_parser.add_argument('--numerator-days',action='store_true')
    yield_days_parser=commands.add_parser('yield-coupon-days')
    yield_days_parser.add_argument('checkout')
    yield_days_parser.add_argument('original')
    yield_days_parser.add_argument('output')
    yield_days_parser.add_argument('--cache',required=True)
    yield_days_parser.add_argument('--offline',action='store_true')
    mapping = commands.add_parser('map')
    mapping.add_argument('checkout')
    mapping.add_argument('output', type=Path)
    mapping.add_argument('--cache', required=True)
    mapping.add_argument('--count', type=int, default=512)
    mapping.add_argument('--offline', action='store_true')
    mapping.add_argument('--only')
    edge_parser=commands.add_parser('edges')
    edge_parser.add_argument('checkout')
    edge_parser.add_argument('output',type=Path)
    edge_parser.add_argument('--cache',required=True)
    edge_parser.add_argument('--offline',action='store_true')
    band_parser=commands.add_parser('bessel-bands')
    band_parser.add_argument('checkout')
    band_parser.add_argument('output',type=Path)
    band_parser.add_argument('--cache',required=True)
    band_parser.add_argument('--offline',action='store_true')
    root_parser=commands.add_parser('roots')
    root_parser.add_argument('checkout')
    root_parser.add_argument('output',type=Path)
    root_parser.add_argument('--cache',required=True)
    root_parser.add_argument('--offline',action='store_true')
    product_parser=commands.add_parser('products')
    product_parser.add_argument('checkout')
    product_parser.add_argument('output',type=Path)
    product_parser.add_argument('--cache',required=True)
    product_parser.add_argument('--offline',action='store_true')
    store_parser=commands.add_parser('product-stores')
    store_parser.add_argument('checkout')
    store_parser.add_argument('output',type=Path)
    store_parser.add_argument('--cache',required=True)
    store_parser.add_argument('--offline',action='store_true')
    stages_parser=commands.add_parser('product-stages')
    stages_parser.add_argument('checkout')
    stages_parser.add_argument('output',type=Path)
    stages_parser.add_argument('--cache',required=True)
    stages_parser.add_argument('--offline',action='store_true')
    poisson_parser=commands.add_parser('poisson')
    poisson_parser.add_argument('checkout')
    poisson_parser.add_argument('output',type=Path)
    poisson_parser.add_argument('--cache',required=True)
    poisson_parser.add_argument('--offline',action='store_true')
    binomial_parser=commands.add_parser('binomial')
    binomial_parser.add_argument('checkout')
    binomial_parser.add_argument('output',type=Path)
    binomial_parser.add_argument('--cache',required=True)
    binomial_parser.add_argument('--offline',action='store_true')
    primitive_parser=commands.add_parser('primitives')
    primitive_parser.add_argument('checkout')
    primitive_parser.add_argument('output',type=Path)
    primitive_parser.add_argument('--cache',required=True)
    primitive_parser.add_argument('--offline',action='store_true')
    power_parser=commands.add_parser('integer-powers')
    power_parser.add_argument('checkout')
    power_parser.add_argument('output',type=Path)
    power_parser.add_argument('--cache',required=True)
    power_parser.add_argument('--offline',action='store_true')
    admission_parser=commands.add_parser('binomial-admission')
    admission_parser.add_argument('checkout')
    admission_parser.add_argument('output',type=Path)
    admission_parser.add_argument('--cache',required=True)
    admission_parser.add_argument('--offline',action='store_true')
    resident_parser=commands.add_parser('resident-stages')
    resident_parser.add_argument('checkout')
    resident_parser.add_argument('output',type=Path)
    resident_parser.add_argument('--cache',required=True)
    resident_parser.add_argument('--library',required=True)
    resident_parser.add_argument('--offline',action='store_true')
    constant_parser=commands.add_parser('reduction-constants')
    constant_parser.add_argument('checkout')
    constant_parser.add_argument('output',type=Path)
    constant_parser.add_argument('--cache',required=True)
    constant_parser.add_argument('--library',required=True)
    constant_parser.add_argument('--offline',action='store_true')
    refine_parser=commands.add_parser('refine-log-constant')
    refine_parser.add_argument('checkout')
    refine_parser.add_argument('prior')
    refine_parser.add_argument('output',type=Path)
    refine_parser.add_argument('--cache',required=True)
    refine_parser.add_argument('--library',required=True)
    refine_parser.add_argument('--offline',action='store_true')
    refine_parser.add_argument('--spread',action='store_true')
    anchor_parser=commands.add_parser('root-anchors')
    anchor_parser.add_argument('checkout')
    anchor_parser.add_argument('output',type=Path)
    anchor_parser.add_argument('--cache',required=True)
    anchor_parser.add_argument('--offline',action='store_true')
    spill_parser=commands.add_parser('root-primitive-stores')
    spill_parser.add_argument('checkout')
    spill_parser.add_argument('output',type=Path)
    spill_parser.add_argument('--cache',required=True)
    spill_parser.add_argument('--library',required=True)
    hyper_return_parser=commands.add_parser('hypergeom-returns')
    hyper_return_parser.add_argument('checkout')
    hyper_return_parser.add_argument('output',type=Path)
    hyper_return_parser.add_argument('--cache',required=True)
    hyper_return_parser.add_argument('--offline',action='store_true')
    poisson_seed_parser=commands.add_parser('poisson-seeds')
    poisson_seed_parser.add_argument('checkout')
    poisson_seed_parser.add_argument('output',type=Path)
    poisson_seed_parser.add_argument('--cache',required=True)
    poisson_seed_parser.add_argument('--offline',action='store_true')
    equivalence_parser=commands.add_parser('root-equivalence')
    equivalence_parser.add_argument('checkout')
    equivalence_parser.add_argument('output',type=Path)
    equivalence_parser.add_argument('--cache',required=True)
    equivalence_parser.add_argument('--offline',action='store_true')
    scale_parser=commands.add_parser('root-scales')
    scale_parser.add_argument('checkout')
    scale_parser.add_argument('output',type=Path)
    scale_parser.add_argument('--cache',required=True)
    scale_parser.add_argument('--offline',action='store_true')
    output_edge_parser=commands.add_parser('root-output-edges')
    output_edge_parser.add_argument('checkout')
    output_edge_parser.add_argument('output',type=Path)
    output_edge_parser.add_argument('--cache',required=True)
    output_edge_parser.add_argument('--offline',action='store_true')
    normalization_parser=commands.add_parser('root-normalization')
    normalization_parser.add_argument('checkout')
    normalization_parser.add_argument('original')
    normalization_parser.add_argument('output',type=Path)
    normalization_parser.add_argument('--cache',required=True)
    normalization_parser.add_argument('--offline',action='store_true')
    log_stage_parser=commands.add_parser('root-log-stage')
    log_stage_parser.add_argument('checkout')
    log_stage_parser.add_argument('original')
    log_stage_parser.add_argument('output',type=Path)
    log_stage_parser.add_argument('--cache',required=True)
    log_stage_parser.add_argument('--library',required=True)
    log_stage_parser.add_argument('--offline',action='store_true')
    public_stage_parser=commands.add_parser('root-public-stages')
    public_stage_parser.add_argument('checkout')
    public_stage_parser.add_argument('original')
    public_stage_parser.add_argument('output',type=Path)
    public_stage_parser.add_argument('--cache',required=True)
    public_stage_parser.add_argument('--offline',action='store_true')
    store_contrast_parser=commands.add_parser('root-store-contrasts')
    store_contrast_parser.add_argument('original')
    store_contrast_parser.add_argument('output',type=Path)
    store_contrast_parser.add_argument('--library',required=True)
    poisson_argument_parser=commands.add_parser('poisson-exp-arguments')
    poisson_argument_parser.add_argument('checkout')
    poisson_argument_parser.add_argument('original')
    poisson_argument_parser.add_argument('output',type=Path)
    poisson_argument_parser.add_argument('--cache',required=True)
    poisson_argument_parser.add_argument('--library',required=True)
    poisson_argument_parser.add_argument('--offline',action='store_true')
    poisson_log_parser=commands.add_parser('poisson-log-difference')
    poisson_log_parser.add_argument('checkout')
    poisson_log_parser.add_argument('original')
    poisson_log_parser.add_argument('output',type=Path)
    poisson_log_parser.add_argument('--cache',required=True)
    poisson_log_parser.add_argument('--library',required=True)
    poisson_log_parser.add_argument('--offline',action='store_true')
    poisson_contrast_parser=commands.add_parser('poisson-log-contrasts')
    poisson_contrast_parser.add_argument('original')
    poisson_contrast_parser.add_argument('output',type=Path)
    poisson_contrast_parser.add_argument('--library',required=True)
    poisson_one_parser=commands.add_parser('poisson-one-primitives')
    poisson_one_parser.add_argument('original')
    poisson_one_parser.add_argument('output',type=Path)
    poisson_one_parser.add_argument('--library',required=True)
    hypergeom_parser=commands.add_parser('hypergeom-products')
    diagonal_parser=commands.add_parser('hypergeom-diagonal')
    diagonal_parser.add_argument('checkout')
    diagonal_parser.add_argument('output',type=Path)
    diagonal_parser.add_argument('--cache',required=True)
    diagonal_parser.add_argument('--offline',action='store_true')
    bessel_admission_parser=commands.add_parser('bessel-admission')
    bessel_admission_parser.add_argument('checkout')
    bessel_admission_parser.add_argument('output',type=Path)
    bessel_admission_parser.add_argument('--cache',required=True)
    bessel_admission_parser.add_argument('--offline',action='store_true')
    hypergeom_parser.add_argument('input')
    hypergeom_parser.add_argument('output')
    log_parser=commands.add_parser('log-consumers')
    log_parser.add_argument('checkout')
    log_parser.add_argument('output',type=Path)
    log_parser.add_argument('--cache',required=True)
    log_parser.add_argument('--offline',action='store_true')
    book_parser=commands.add_parser('workbooks')
    for name in ('checkout','output','executable','powershell'):book_parser.add_argument(name)
    book_parser.add_argument('--root-transfer',action='store_true')
    book_parser.add_argument('--rank-dependencies',action='store_true')
    book_parser.add_argument('--rank-edges')
    args = parser.parse_args()
    if args.command == 'yield-multi-coupon':
        yield_multi_coupon(args.checkout, args.output, args.cache, args.offline)
        raise SystemExit
    if args.command=='weeknum-controls':
        weeknum_controls(args.checkout,args.output,args.cache,args.bank_trace,args.offline,args.repeat)
        raise SystemExit
    if args.command=='sinh-small-stages':
        sinh_small_stages(args.checkout,args.output,args.cache,args.library,args.offline,args.bank_trace,args.count)
        raise SystemExit
    if args.command=='duration-off-coupon':
        duration_off_coupon(args.checkout,args.output,args.cache,args.library,args.offline,args.count,args.direct_days)
        raise SystemExit
    if args.command=='duration-two-coupon-transfer':
        duration_two_coupon_transfer(args.checkout,args.output,args.cache,args.offline,args.trials)
        raise SystemExit
    if args.command=='duration-zero-coupon':
        duration_zero_coupon(args.checkout,args.output,args.cache,args.library,args.offline,args.count,args.holdout)
        raise SystemExit
    if args.command=='yield-coupon-days':
        yield_coupon_days(args.checkout,args.original,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='yield-coupon-programs':
        yield_coupon_programs(args.original,args.output,args.library,args.admitted,args.day_boundaries,args.nasd_days,args.numerator_days)
        raise SystemExit
    if args.command=='yield-single-coupon':
        yield_single_coupon(args.checkout,args.output,args.cache,args.offline,args.month_ends,args.date_holdout)
        raise SystemExit
    if args.command=='edges':
        edges(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='bessel-bands':
        bessel_bands(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='roots':
        root_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='products':
        product_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='product-stores':
        product_store_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='product-stages':
        product_stage_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='poisson':
        poisson_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='binomial':
        binomial_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='primitives':
        primitive_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='integer-powers':
        integer_power_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='binomial-admission':
        binomial_admission_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='resident-stages':
        resident_stage_study(args.checkout,args.output,args.cache,args.library,args.offline)
        raise SystemExit
    if args.command=='reduction-constants':
        reduction_constant_study(args.checkout,args.output,args.cache,args.library,args.offline)
        raise SystemExit
    if args.command=='refine-log-constant':
        refine_log_constant(args.checkout,args.prior,args.output,args.cache,args.library,args.offline,args.spread)
        raise SystemExit
    if args.command=='root-anchors':
        root_anchor_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-primitive-stores':
        root_primitive_spill_study(args.checkout,args.output,args.cache,args.library)
        raise SystemExit
    if args.command=='hypergeom-returns':
        hypergeom_return_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='poisson-seeds':
        poisson_seed_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-equivalence':
        root_equivalence_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-scales':
        root_scale_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-output-edges':
        root_output_edges(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-normalization':
        root_normalization_study(args.checkout,args.original,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-log-stage':
        root_log_stage(args.checkout,args.original,args.output,args.cache,args.library,args.offline)
        raise SystemExit
    if args.command=='root-public-stages':
        root_public_stages(args.checkout,args.original,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='root-store-contrasts':
        root_store_contrasts(args.original,args.output,args.library)
        raise SystemExit
    if args.command=='poisson-exp-arguments':
        poisson_exp_argument_study(args.checkout,args.original,args.output,args.cache,args.library,args.offline)
        raise SystemExit
    if args.command=='poisson-log-difference':
        poisson_log_difference_study(args.checkout,args.original,args.output,args.cache,args.library,args.offline)
        raise SystemExit
    if args.command=='poisson-log-contrasts':
        poisson_log_primitive_contrasts(args.original,args.output,args.library)
        raise SystemExit
    if args.command=='poisson-one-primitives':
        poisson_one_primitive_contrasts(args.original,args.output,args.library)
        raise SystemExit
    if args.command=='bessel-admission':
        bessel_admission_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='hypergeom-diagonal':
        hypergeom_diagonal_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit(0)
    if args.command=='hypergeom-products':
        hypergeom_product_constraints(args.input,args.output)
        raise SystemExit
    if args.command=='log-consumers':
        log_consumer_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='workbooks':
        make_workbooks(args.checkout,args.output,args.executable,args.powershell,args.root_transfer,args.rank_dependencies,args.rank_edges)
        raise SystemExit
    if args.count < 32 or args.count > 4096 or args.count & (args.count-1):
        parser.error('--count must be a power of two from32 through4096')
    scan(args.checkout, args.output, args.cache, args.count, args.offline, args.only)
