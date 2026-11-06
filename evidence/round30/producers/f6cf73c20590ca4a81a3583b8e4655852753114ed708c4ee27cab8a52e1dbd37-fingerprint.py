"""Consecutive-double Excel fingerprints, using the existing scalar service.

Original MIT tooling. Inputs, typed answers and worker provenance are retained
in small gzip files. Spectral peaks are observations, not a count of hidden
stores: only a discriminated, fully matching computation tree can supply that.
No coefficient, threshold or tree is ranked by the number of matching answers.
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
                      analysis_model='gpt-6.1-sol (session launcher)',
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

    def ask(self, function, rows):
        request = dict(function=function, rows=rows)
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
                reply = self.client.ask(function, rows)
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


def poisson_study(checkout,output,cache,offline=False,count=512):
    """Fixed zero/one-count textbook products; independent EXP readings.

    Measured public EXP is an explicitly conditional binary64 leaf. Surviving
    products identify only this stated consumer on the sampled windows.
    """
    client,_=load_service(checkout)
    service=CachedService(client.Client(timeout=900),cache,offline)
    centers=(.13,1.3,3.7,10.,50.)
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
                tally=Counter();first=None
                for x,e,a in zip(xs,exponential,answers[j*len(xs):(j+1)*len(xs)]):
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
                    if key=='mismatched'and first is None:first=dict(input=x,expo=e['bits'],excel=a['bits'],candidate=predicted)
                profiles.append(dict(shape=shape,model=model,counts=tally,first_contradiction=first,
                                     consistent=tally['matched']==len(xs)))
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
        limitations='Positive finite small means and counts0..2 only. EXP readings are public scalar binary64 leaves, not recovered private Ext80 seeds. A mismatch rejects only the specified consumer. No threshold, polynomial, recurrence selector or kernel is fit to counts.')
    write(output,result)
    return result


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


def positive_extended_cell(native, graph, answer, center,decreasing=False):
    """Conditional inverse of one fixed positive extended-input primitive.

    Find both ends of an output's preimage on the PC64 input lattice. The
    local monotonicity of this fixed instruction chain is a declared premise;
    Excel is never queried with an inferred extended operand. No logarithm
    producer is assumed. Empty preimages are kept as contradictions.
    """
    from x87 import encode,decode
    significand=EXTENDED64.floor_coordinate(EXTENDED64.round(center))
    radius=1<<18
    lo,hi=significand-radius,significand+radius
    if not 0<lo<hi<EXTENDED64.last:
        raise ValueError('inverse bracket leaves positive finite extended values')
    def pack(sig):return encode(EXTENDED64.value(sig))
    def evaluate(sig):
        return native.publish(native.evaluate(graph,[],0x133f,extended=pack(sig))[:10])
    direction=-1 if decreasing else 1
    target=direction*int(answer,16)
    if not direction*evaluate(lo)<target<direction*evaluate(hi):
        raise ValueError('inverse bracket does not enclose requested output')
    def first(strict):
        a,b=lo,hi
        while b-a>1:
            m=(a+b)//2;got=direction*evaluate(m)
            if got>target if strict else got>=target:b=m
            else:a=m
        return b
    lower,upper=first(False),first(True)-1
    if lower>upper:return None
    if direction*evaluate(lower)!=target or direction*evaluate(upper)!=target:
        raise ValueError('inverse endpoint does not publish requested output')
    return dict(low=pack(lower).hex(),high=pack(upper).hex(),
                predecessor_output=f'{evaluate(lower-1):016x}',
                successor_output=f'{evaluate(upper+1):016x}')


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
