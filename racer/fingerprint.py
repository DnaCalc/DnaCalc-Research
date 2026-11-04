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
        if PRODUCER_BINDINGS!={p:sha(p)for p in PRODUCER_BINDINGS}:
            raise ValueError('analysis source changed; raw replies remain retained, regenerate this analysis')
        document=dict(document,producer_bindings=PRODUCER_BINDINGS,
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


def append_question(book,name,formula,role,**fields):
    row=len(book)+1
    book.append(dict(id=f'r26-{name}-{row:05}',row=row,column=4,formula=formula,role=role,**fields))
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


def make_workbooks(checkout,output,executable,powershell,root_transfer=False):
    checkout,output,executable,powershell=map(lambda p:Path(p).resolve(),(checkout,output,executable,powershell))
    output.mkdir(parents=True,exist_ok=False)
    books,requests=root_transfer_questions()if root_transfer else workbook_questions(checkout)
    generator=checkout/'crates/calc-func/tests/fixtures/make-catalogue-fixtures.ps1'
    provenance=subprocess.check_output([str(executable),'provenance'],timeout=60)
    (output/'prediction-binding.json').write_bytes(provenance)
    inventory=[]
    for name,questions in books.items():
        directory=output/name;directory.mkdir()
        case='dcj9k_r26_'+name
        doc=directory/'workbook-questions.json';write(doc,dict(sheet=case,allow_refusals=True,questions=questions))
        section=directory/('sheet-'+case+'.txt')
        section.write_text('# R26 exact operands and unanswered fingerprint questions.\n'+'\n'.join(q['formula']for q in questions)+'\n')
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
    manifest=dict(format='dnacalc-r26-fingerprint-questions/1',new_excel_answers=0,workbooks=inventory,
                  scalar_requests=requests,generator_sha256=sha(generator),producer_sha256=sha(__file__),
                  executable_sha256=sha(executable),
                  files_sha256={str(p.relative_to(output)):sha(p)for p in output.rglob('*')if p.is_file()},
                  limits=('New workbook only. Reference-cell inputs use exact integer limbs and echoes. Decimal arrays have their own INDEX operand echo, without an assumed binary64 value. LN, reciprocal, mean and EXP occupy separate cells. Every workbook outcome awaits capture; scalar observations do not establish input-route transfer. Predictions are application calculations, not Excel readings.'if root_transfer else
                          'New workbooks only. Every numerical operand uses exact small integer limbs and is echoed. Predictions are application calculations, not Excel readings. QR rank thresholds and complete arithmetic remain unknown; no tolerance or coefficient fitted.'))
    write(output/'manifest.json',manifest)
    if root_transfer:
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
        *[s+'\n' for s in result['limitations']],
        '| Function | Status | Rounding | Stored stages | Window peaks | Route edges |',
        '|---|---|---|---|---|---|']
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
    log_parser=commands.add_parser('log-consumers')
    log_parser.add_argument('checkout')
    log_parser.add_argument('output',type=Path)
    log_parser.add_argument('--cache',required=True)
    log_parser.add_argument('--offline',action='store_true')
    book_parser=commands.add_parser('workbooks')
    for name in ('checkout','output','executable','powershell'):book_parser.add_argument(name)
    book_parser.add_argument('--root-transfer',action='store_true')
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
    if args.command=='log-consumers':
        log_consumer_study(args.checkout,args.output,args.cache,args.offline)
        raise SystemExit
    if args.command=='workbooks':
        make_workbooks(args.checkout,args.output,args.executable,args.powershell,args.root_transfer)
        raise SystemExit
    if args.count < 32 or args.count > 4096 or args.count & (args.count-1):
        parser.error('--count must be a power of two from32 through4096')
    scan(args.checkout, args.output, args.cache, args.count, args.offline, args.only)
