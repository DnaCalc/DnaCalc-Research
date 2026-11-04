"""Independent rational replay of normal-CDF sketches and their complete bounds.

No racer/Fpu imports. Transcendentals are explicitly conditional concrete leaves,
bound by full input/output encodings and provider; they are not proved emulations.
Every arithmetic, square root, store, constant load and published contradiction
is checked without the search numeric substrate. No tolerance or fitting.
"""
from functools import lru_cache
from fractions import Fraction as R
import itertools
import json
from pathlib import Path
import time

from exact import (BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P53,
                   binary64, operation, two, Span, inverse)
from check import require, records, digest, tagged_value, write_result

MODELS = ('strict', 'stored_pc64', 'resident_pc64', 'stored_pc53', 'resident_pc53')
FIELDS = ('sketch', 'extended_constants', 'normalization', 'reduction', 'polynomial',
          'quotient', 'exponent', 'reconstruction', 'exponential', 'split', 'association')
SOURCE_HASH = digest(Path(__file__))


def grammar(wrappers=False, field=False, hart=False):
    for sketch in ((2,) if hart else (0, 1)):
        for ext, norm, reduce, poly, reconstruct, provider, split, assoc in itertools.product(
                (False, True), range(3 if sketch == 1 else 1),
                MODELS[:3] if sketch == 1 else ('strict',), MODELS, MODELS,
                range(4), range(3), range(3)):
            for reflection in range(4 if wrappers else 1):
                values = (sketch, ext, norm, reduce, poly, poly, reconstruct,
                          reconstruct, provider, split, assoc)
                yield values + ((reflection,) if field else ())


@lru_cache(maxsize=64)
def sqrt_round(x, model):
    """Nearest lattice root via exact squared comparisons and midpoint ties."""
    fmt = B if model == 'strict' else P53 if model.endswith('pc53') else E
    if model == 'strict':
        x = B.round(x)
    require(x >= 0, 'negative root')
    lo, hi = 0, fmt.last
    while lo < hi:
        mid = (lo+hi+1)//2
        if fmt.value(mid)**2 <= x:
            lo = mid
        else:
            hi = mid-1
    low, high = fmt.value(lo), fmt.value(lo+1)
    midpoint = (low+high)/2
    value = low if x < midpoint**2 or (x == midpoint**2 and lo % 2 == 0) else high
    return B.round(value) if model.startswith('stored') else value


def encoding(value):
    q = B.round(value)
    require(q is not None, 'nonfinite publication outside checked scope')
    if abs(q) < two(-1022):
        return '0x0000000000000000'
    coordinate = B.floor_coordinate(abs(q))
    return f'0x{((1 << 63) if q < 0 else 0) | coordinate:016x}'


class Replay:
    def __init__(self, nodes, constants, helpers):
        self.nodes, self.constants, self.helpers = nodes, constants, helpers
        self.values = []
        self.helper_uses = 0

    def node(self, op, value, **attributes):
        i = len(self.values)
        require(i < len(self.nodes), 'missing trace node')
        node = self.nodes[i]
        require({k:v for k,v in node.items() if k != 'value'} == dict(op=op, **attributes),
                f'node {i} differs from independently generated sketch')
        require(tagged_value(node['value']) == value, f'node {i} exact value mismatch')
        binary = op in ('literal', 'input') or (op == 'constant' and not attributes['extended'])
        if op in ('add', 'mul', 'div', 'sqrt'):
            binary = attributes['model'] == 'strict' or attributes['model'].startswith('stored')
        if op == 'exp':
            binary = attributes['provider'] != 0
        if op not in ('neg',):
            require(node['value'].startswith('binary64:') == binary, f'node {i} wrong value representation')
        self.values.append(value)
        return i

    def literal(self, bits):
        return self.node('literal', binary64(bits), bits=bits)

    def const(self, table, index, extended):
        c = self.constants['tables'][table][index]
        value = R(c['decimal'])
        value = E.round(value) if extended else B.round(value)
        return self.node('constant', value, table=table, index=index, extended=extended)

    def op(self, kind, model, a, b):
        value = operation(self.values[a], self.values[b], kind, model)
        require(value is not None, 'nonfinite arithmetic not certified')
        return self.node(kind, value, model=model, args=[a,b])

    def neg(self, a):
        return self.node('neg', -self.values[a], args=[a])

    def sub(self, model, a, b):
        return self.op('add', model, a, self.neg(b))

    def sqrt(self, a, model):
        return self.node('sqrt', sqrt_round(self.values[a], model), model=model, args=[a])

    def floor(self, a):
        return self.node('floor', R(self.values[a]//1), args=[a])

    def exp(self, a, provider):
        key = f'{provider}:'+self.nodes[a]['value']
        require(key in self.helpers, 'unfrozen transcendental leaf')
        leaf = self.helpers[key]
        require(leaf['provider'] == provider and leaf['input'] == self.nodes[a]['value'], 'helper input binding')
        self.helper_uses += 1
        return self.node('exp', tagged_value(leaf['output']), args=[a], provider=provider, concrete_leaf=key)

    def numerator(self, table, z, ext, model):
        coefficients = self.constants['tables'][table]
        result = self.const(table, len(coefficients)-1, ext)
        for i in range(len(coefficients)-1):
            result = self.op('mul', model, result, z)
            result = self.op('add', model, result, self.const(table, i, ext))
        return result

    def denominator(self, table, z, ext, model):
        count = len(self.constants['tables'][table])
        result = z
        for i in range(count):
            result = self.op('add', model, result, self.const(table, i, ext))
            if i+1 < count:
                result = self.op('mul', model, result, z)
        return result

    def ordinary_polynomial(self, table, z, ext, model):
        count = len(self.constants['tables'][table])
        result = self.const(table,count-1,ext)
        for index in reversed(range(count-1)):
            result = self.op('add',model,self.op('mul',model,result,z),self.const(table,index,ext))
        return result

    def triple(self, a, b, c, model, association):
        if association == 0:
            return self.op('mul', model, self.op('mul', model, a, b), c)
        if association == 1:
            return self.op('mul', model, a, self.op('mul', model, b, c))
        return self.op('mul', model, c, self.op('mul', model, a, b))


def replay(program, witness, nodes, constants, helpers):
    p, a = program, witness['args']
    require(witness['function'] == 'NORM.S.DIST' and witness['result_index'] is None and len(a) == 2,
            'unsupported original call')
    x = binary64(a[0])
    require(x == 0 or abs(x) >= two(-1022), 'subnormal input control')
    flag = a[1] if isinstance(a[1], bool) else binary64(a[1])
    require(flag != 0, 'density outside scope')
    m = Replay(nodes, constants, helpers)
    m.node('input', x, bits=a[0][2:].lower())
    half = m.literal('3fe0000000000000')
    one = m.literal('3ff0000000000000')
    positive = m.literal(f'{int(a[0],16)&((1<<63)-1):016x}')
    y = positive
    if p['sketch'] == 1:
        if p['normalization'] == 2:
            c = m.sqrt(m.literal('4000000000000000'), p['reduction'])
        else:
            c = m.const('sqrt2' if p['normalization'] == 0 else 'invsqrt2', 0, p['extended_constants'])
        y = m.op('mul' if p['normalization'] == 1 else 'div', p['reduction'], positive, c)
    # The decimal0.66291 is loaded as binary64; derive/check its encoding independently.
    central = p['sketch'] != 2 and m.values[y] <= (B.round(R('0.66291')) if p['sketch'] == 0 else R(15,32))
    middle = p['sketch'] == 2 or m.values[y] <= (binary64('4016a09e667f3bcd') if p['sketch'] == 0 else R(4))
    region = 'central' if central else 'middle' if middle else 'tail'
    prefix = 'direct_' if p['sketch'] == 0 else 'erfc_' if p['sketch'] == 1 else 'hart_'
    ext, poly, div, reconstruct = (p[k] for k in ('extended_constants','polynomial','quotient','reconstruction'))
    if central:
        z = m.op('mul', poly, y, y)
        n = m.numerator(prefix+'a', z, ext, poly)
        d = m.denominator(prefix+'b', z, ext, poly)
        if p['association'] == 0:
            delta = m.op('div', div, m.op('mul', div, y, n), d)
        elif p['association'] == 1:
            delta = m.op('mul', div, y, m.op('div', div, n, d))
        else:
            delta = m.op('mul', div, m.op('div', div, y, d), n)
        if p['sketch'] == 0:
            answer = m.sub(reconstruct, half, delta) if x < 0 else m.op('add', reconstruct, half, delta)
        else:
            wrapper = m.sub(reconstruct, one, delta) if x < 0 else m.op('add', reconstruct, one, delta)
            answer = m.op('mul', reconstruct, half, wrapper)
    elif abs(x) > 40:
        answer = m.literal('0000000000000000' if x < 0 else '3ff0000000000000')
    else:
        if p['sketch'] == 2:
            require(abs(x)<B.round(R('7.07106781186547')),'outside published Hart body')
            n = m.ordinary_polynomial('hart_p',y,ext,poly)
            d = m.ordinary_polynomial('hart_q',y,ext,poly)
            ratio = m.op('div',div,n,d)
        elif middle:
            n = m.numerator(prefix+'c', y, ext, poly)
            d = m.denominator(prefix+'d', y, ext, poly)
            ratio = m.op('div', div, n, d)
        else:
            z = m.op('div', poly, one, m.op('mul', poly, y, y))
            n = m.numerator(prefix+'p', z, ext, poly)
            d = m.denominator(prefix+'q', z, ext, poly)
            ratio = m.op('div', div, m.op('mul', div, z, n), d)
            c = m.const(prefix+'sqrpi', 0, ext)
            ratio = m.op('div', div, m.sub(div, c, ratio), y)
        exponent, provider = p['exponent'], p['exponential']
        scale = one if p['sketch'] == 1 else half
        if p['split'] == 0:
            sq = m.op('mul', exponent, y, y)
            sq = m.op('mul', exponent, sq, scale)
            ea, eb = m.neg(sq), m.literal('0000000000000000')
        else:
            grid = m.literal('4030000000000000' if p['split'] == 1 else '3ff999999999999a')
            t = m.op('div', exponent, m.floor(m.op('mul', exponent, y, grid)), grid)
            lo = m.sub(exponent, y, t)
            hi = m.op('add', exponent, y, t)
            delta = m.op('mul', exponent, lo, hi)
            delta = m.op('mul', exponent, delta, scale)
            sq = m.op('mul', exponent, t, t)
            sq = m.op('mul', exponent, sq, scale)
            ea, eb = m.neg(sq), m.neg(delta)
        ea, eb = m.exp(ea, provider), m.exp(eb, provider)
        small = m.triple(ea, eb, ratio, reconstruct, p['association'])
        if p['sketch'] == 1:
            small = m.op('mul', reconstruct, small, half)
        answer = m.sub(reconstruct, one, small) if x > 0 else small
    if region == 'middle':
        reflection = p.get('reflection', 0)
        if reflection == 1:
            answer = m.sub(reconstruct, one, m.sub(reconstruct, one, answer))
        elif reflection == 2:
            answer = m.sub(reconstruct, half, m.sub(reconstruct, half, answer))
        elif reflection == 3:
            twice = m.literal('4000000000000000')
            doubled = m.op('mul', reconstruct, answer, twice)
            erf = m.sub(reconstruct, one, doubled)
            erfc = m.sub(reconstruct, one, erf)
            answer = m.op('mul', reconstruct, half, erfc)
    require(len(m.values) == len(nodes), 'extra trace nodes')
    return encoding(m.values[answer]), region, len(nodes), m.helper_uses


def check(directory, constant_path, trace, output):
    start, directory = time.monotonic(), Path(directory)
    constants = json.loads(Path(constant_path).read_text())
    for table in constants['tables'].values():
        for c in table:
            require(binary64(c['binary64']) == B.round(R(c['decimal'])), 'binary64 constant load')
            require(tagged_value('ext80-le:'+c['ext80_le']) == E.round(R(c['decimal'])), 'extended constant load')
    programs = list(records(directory/'programs.jsonl'))
    field = 'reflection' in programs[0]
    hart = programs[0]['sketch'] == 2
    wrappers = len(programs) == 72000 or hart
    declared = list(grammar(wrappers,field,hart))
    fields = FIELDS + (('reflection',) if field else ())
    require(len(programs) == len(declared) == (7200 if hart else 72000 if wrappers else 18000), 'incomplete normal grammar')
    for i, (p, values) in enumerate(zip(programs, declared)):
        require(p['id'] == i and tuple(p[k] for k in fields) == values, 'grammar member missing or reordered')
    original = list(records(trace))
    helpers = json.loads((directory/'concrete-helpers.json').read_text())
    summary = json.loads((directory/'summary.json').read_text())
    seen, nodes_checked, helper_uses, groups = set(), 0, 0, 0
    for kill in records(directory/'first-kills.jsonl'):
        p, witness = kill['program'], kill['witness']
        require(original[kill['trace_index']] == witness, 'original observation binding')
        require(witness['reading']['excel_build'] == summary['build'], 'mixed environment contradiction')
        require(programs[p['id']] == p and p['id'] in kill['candidate_ids'], 'representative program binding')
        got, region, count, calls = replay(p, witness, kill['nodes'], constants, helpers)
        require(region == kill['region'] and got == kill['got'] != kill['expected'] == witness['expected'],
                'not an exact original contradiction')
        for id in kill['candidate_ids']:
            require(id not in seen, 'candidate killed twice')
            seen.add(id)
            a, b = dict(programs[id]), dict(p)
            a.pop('id'); b.pop('id')
            if p['sketch'] == 0 and region == 'central':
                for k, value in [('exponent','strict'),('exponential',0),('split',0)]:
                    a[k] = b[k] = value
            require(a == b, 'group contains an unproved alternative')
        nodes_checked += count
        helper_uses += calls
        groups += 1
    survivors = json.loads((directory/'survivors.json').read_text())
    require(seen.isdisjoint(p['id'] for p in survivors), 'survivor was killed')
    require(seen | {p['id'] for p in survivors} == set(range(len(programs))), 'unchecked program')
    require(not survivors, 'survivor replay must be separately checked before acceptance')
    require(digest(Path(__file__)) == SOURCE_HASH, 'normal checker changed during run')  # ubs:ignore — public artifact integrity, no secret comparison
    result = dict(complete=True, programs=len(programs), exact_contradictions=len(seen),
        trace_groups=groups, nodes_checked=nodes_checked, conditional_concrete_helper_uses=helper_uses,
        seconds=time.monotonic()-start, normal_checker_sha256=SOURCE_HASH,
        inputs={str(p):digest(Path(p)) for p in [constant_path,trace,directory/'programs.jsonl',directory/'first-kills.jsonl',directory/'concrete-helpers.json']},
        claim='Independent exact arithmetic, full-width traces, original witness binding and complete finite grammar; transcendental leaves remain conditional, not proved target helpers.',
        global_equivalence=False, fresh_excel_validation=0)
    write_result(output, result)


def linear_check(trace, directory, output):
    """Check arbitrary shared-coefficient exclusions with separate lattice arithmetic."""
    from collections import Counter
    start = time.monotonic()
    directory = Path(directory)
    profiles = json.loads((directory/'profiles.json').read_text())
    declared = [(form,models) for form,count in [(0,2),(1,3),(2,3)]
                for models in itertools.product(range(5),repeat=count)]
    require(len(profiles) == len(declared) == 275,'incomplete linear family')
    for index,(profile,(form,models)) in enumerate(zip(profiles,declared)):
        require(profile == dict(id=index,form=form,models=list(models)),'undeclared linear profile')
    original = list(records(trace))
    seen, reasons, builds = set(),Counter(),set()
    for cert in records(directory/'certificates.jsonl'):
        p = cert['profile']
        require(p == profiles[p['id']] and p['id'] not in seen,'wrong or duplicate linear profile')
        seen.add(p['id']);current=None;first=True
        require(1<=len(cert['witnesses'])<=3,'not a small sufficient subset')
        for entry in cert['witnesses']:
            row = original[entry['trace_index']]
            require(row == entry['witness'],'changed original observation')
            require(dict(reading=row['reading'],source=row['source'],flag=row['args'][1]) == cert['group'], 'coefficient group identity')
            require(row['function']=='NORM.S.DIST' and row['result_index'] is None and len(row['args'])==2,'unsupported surface')
            flag = row['args'][1] if isinstance(row['args'][1],bool) else binary64(row['args'][1])
            require(flag != 0,'density row in CDF certificate')
            x,y = binary64(row['args'][0]),binary64(row['expected'])
            require(two(-1022)<=abs(x)<=two(-27) and 0<y<1,'linear certificate domain')
            models = [MODELS[m] for m in p['models']]
            span = B.pullback(Span(y,y))
            if p['form']==0:
                span = inverse(span,R(1,2),'add',True,models[1])
            elif p['form']==1:
                span = inverse(span,R(1,2),'mul',True,models[2])
                span = inverse(span,R(1),'add',True,models[1]) if span else None
            else:
                span = inverse(span,R(1,2),'add',True,models[2])
                span = inverse(span,R(1,2),'mul',True,models[1]) if span else None
            span = inverse(span,x,'mul',True,models[0]) if span else None
            current = span if first else current.meet(span) if current else None
            first=False;builds.add(row['reading']['excel_build'])
        require(current is None or E.coordinates(current) is None,'shared coefficient is not excluded')
        require((current.json() if current else None)==cert['remaining_interval'],'linear inverse endpoint disagreement')
        reason='empty_real_intersection' if current is None else 'empty_ext80_lattice'
        require(cert['reason']==reason,'linear exclusion reason')
        reasons[reason]+=1
    require(seen == set(range(275)),'linear profile missing')
    require(digest(Path(__file__))==SOURCE_HASH,'checker changed during linear proof')  # ubs:ignore — public artifact integrity, no secret comparison
    write_result(output,dict(complete=True,profiles_checked=len(seen),reasons=reasons,builds=sorted(builds),
        seconds=time.monotonic()-start,normal_checker_sha256=SOURCE_HASH,
        inputs={str(p):digest(Path(p)) for p in [trace,directory/'profiles.json',directory/'certificates.jsonl']},
        claim='Independent exact shared-coefficient exclusion, with full typed original witness binding; no coefficient value selected and no transcendental leaf assumption.',
        limits='Only three declared linear forms, fixed tiny domain, five arithmetic models and shared coefficient within a reading/source/flag group; no general normal-program or build20527 claim.',new_excel_answers=0))


def symmetry_span(entry, profile, negative):
    y=binary64(entry['witness']['expected'])
    require(0<y<1,'symmetry endpoints outside this unit')
    span=B.pullback(Span(y,y));ms=[MODELS[m]for m in profile['models']]
    form=profile['form']
    if form in (0,1):
        if negative == (form==1):
            span=inverse(span,R(1),'add',True,ms[0])
            span=span.linear(R(-1))if span else None
    elif form==2:
        span=inverse(span,R(1,2),'add',True,ms[0 if negative else 1])
        span=span.linear(R(-1))if span and negative else span
    else:
        base=0 if negative else 2
        span=inverse(span,R(1,2),'mul',True,ms[base+1])
        span=inverse(span,R(1),'add',True,ms[base])if span else None
        span=span.linear(R(-1))if span and negative else span
    return span


def symmetry_meet(pair, profile):
    a,b=symmetry_span(pair['negative'],profile,True),symmetry_span(pair['positive'],profile,False)
    meet=a.meet(b)if a else None
    return meet.meet(Span(R(0),R(1,2)if profile['form']==2 else R(1)))if meet else None


def symmetry_check(trace,directory,output):
    start=time.monotonic();directory=Path(directory);original=list(records(trace))
    profiles=json.loads((directory/'profiles.json').read_text())
    declared=[(form,models)for form,count in [(0,1),(1,1),(2,2),(3,4)]
              for models in itertools.product(range(5),repeat=count)]
    require(len(profiles)==len(declared)==660,'incomplete sign/reconstruction grammar')
    for index,(p,(form,models))in enumerate(zip(profiles,declared)):
        require(p==dict(id=index,form=form,models=list(models)),'undeclared symmetry profile')
    pairs=json.loads((directory/'pairs.json').read_text())
    for pair in pairs:
        positive,negative=(pair[k]for k in ('positive','negative'))
        a,b=positive['witness'],negative['witness']
        require(original[positive['trace_index']]==a and original[negative['trace_index']]==b,'symmetry original binding')
        require(a['function']==b['function']=='NORM.S.DIST' and a['result_index'] is b['result_index'] is None,'symmetry surface')
        require(a['reading']==b['reading'] and a['source']==b['source'] and a['args'][1]==b['args'][1],'symmetry environment/source/flag')
        x,y=int(a['args'][0],16),int(b['args'][0],16)
        require(x>0 and x^(1<<63)==y,'not exact sign companions')
        require(0<binary64(b['expected'])<R(1,2)<binary64(a['expected'])<1,'ordinary sign targets')
    seen=set();nodes=0
    for cert in records(directory/'certificates.jsonl'):
        p=cert['profile'];cid=p['id'];require(p==profiles[cid]and cid not in seen,'duplicate/unknown sign profile')
        seen.add(cid);require(cert['pair']==pairs[cert['order']],'wrong killing pair')
        interval=symmetry_meet(cert['pair'],p)
        require(interval is None or E.coordinates(interval)is None,'sign relation is not excluded')
        require((interval.json()if interval else None)==cert['remaining_interval'],'sign inverse endpoints')
        require(cert['reason']==('empty_real_intersection'if interval is None else'empty_ext80_lattice'),'sign reason mismatch')
    survivors=json.loads((directory/'survivors.json').read_text())
    constraints={(r['profile_id'],r['positive_trace_index'],r['negative_trace_index']):r for r in records(directory/'survivor-constraints.jsonl')}
    require(len(constraints)==len(survivors)*len(pairs),'missing or duplicate survivor constraint')
    for p in survivors:
        require(p==profiles[p['id']]and p['id']not in seen,'invalid surviving sign relation')
        seen.add(p['id'])
        for pair in pairs:
            interval=symmetry_meet(pair,p);require(interval is not None and E.coordinates(interval)is not None,'survivor infeasible')
            key=(p['id'],pair['positive']['trace_index'],pair['negative']['trace_index'])
            require(interval.json()==constraints[key]['interval'],'surviving inverse endpoint mismatch')
    require(seen==set(range(660)),'unchecked sign profile')
    require(digest(Path(__file__))==SOURCE_HASH,'checker changed during sign proof')  # ubs:ignore — public artifact integrity, no secret comparison
    write_result(output,dict(complete=True,profiles=660,excluded=660-len(survivors),surviving_constraint_profiles=len(survivors),
        original_pairs=len(pairs),survivor_constraints=len(constraints),normal_checker_sha256=SOURCE_HASH,seconds=time.monotonic()-start,
        inputs={str(p):digest(Path(p))for p in [trace,directory/'profiles.json',directory/'pairs.json',directory/'certificates.jsonl',directory/'survivors.json',directory/'survivor-constraints.jsonl']},
        claim='Independent exact sign/reconstruction exclusions and feasible shared-latent constraints, with original paired observations bound by exact arguments/environment/source. Surviving relations are necessary constraints, not complete normal CDF programs.',
        new_excel_answers=0,fresh_validation=0))


def complement_sets(negative_targets):
    """Exact possible publications for every admissible latent, no ULP window."""
    latent=B.pullback(negative_targets);require(latent is not None,'no admissible negative targets')
    bounds=E.coordinates(latent);require(bounds is not None,'no latent lattice value')
    low,high=map(E.value,bounds);result=[]
    for model in MODELS:
        ends=[B.round(operation(R(1),-s,'add',model))for s in (low,high)]
        first,last=sorted(B.floor_coordinate(x)for x in ends)
        outcomes=[]
        for k in range(first,last+1):
            target=B.value(k);required=inverse(B.pullback(Span(target,target)),R(1),'add',True,model)
            required=required.linear(R(-1))if required else None
            intersection=latent.meet(required)
            if intersection is not None and E.coordinates(intersection)is not None:
                outcomes.append(f'0x{k:016x}')
        require(outcomes,'no complementary outcome')
        result.append(outcomes)
    return result


def question_check(trace,packet,output):
    packet=Path(packet);original=list(records(trace))
    base=json.loads((packet/'requests-NORM.S.DIST.json').read_text())['witnesses']
    extra=json.loads((packet/'requests-NORM.S.DIST-lattice.json').read_text())['witnesses']
    combined=json.loads((packet/'requests-NORM.S.DIST-all.json').read_text())['witnesses']
    require(combined==base+extra and len(base)==16 and len(extra)==8,'combined request binding')
    require(len({q['id']for q in combined})==len({tuple(q['args'])for q in combined})==24,'duplicate normal question')
    for q in combined:
        require(set(q)=={'id','args','input_transport'} and q['input_transport']=='cell_value2_bulk','request includes answers or an unsupported transport')
        require(len(q['args'])==2 and q['args'][1]=='0x3ff0000000000000','CDF flag/types')
        for a in q['args']:binary64(a)
    predictions=list(records(packet/'NORM-predictions.jsonl'))
    require([r['id']for r in predictions]==[q['id']for q in base],'base prediction IDs')
    selected=json.loads((packet/'NORM-selection.json').read_text())['control_and_discriminator_selection']
    comparisons=0
    for index,entry in enumerate(selected):
        negative,positive=base[2*index:2*index+2]
        old=original[entry['trace_index']]
        a,b=predictions[2*index:2*index+2]
        require(a['historical_observation']==old and negative['args']==old['args']and a['predicted']==old['expected'],'historical anchor binding')
        require(positive['args'][0]==f"0x{int(negative['args'][0],16)^(1<<63):016x}"and positive['args'][1]==negative['args'][1],'question signs differ by more than sign bit')
        require(b['condition']['negative_request']==negative['id']and b['condition']['negative_outcome']==old['expected'],'conditional prediction prerequisite')
        expected=complement_sets(Span(binary64(old['expected']),binary64(old['expected'])))
        require(b['alternatives']==[dict(profile_id=i,possible_outcomes=v)for i,v in enumerate(expected)],'base set prediction differs from independent preimage')
        comparisons+=5
    document=json.loads((packet/'NORM-lattice-predictions.json').read_text())
    p0,w=two(-54),two(-66)
    cases=[Span(p0-w,p0,True,False),Span(p0,p0),Span(p0,p0+w,False,True)]
    require(len(document['cases'])==3,'missing frozen lattice case')
    for actual,span in zip(document['cases'],cases):
        require(actual['negative_outcome_interval']==span.json(),'altered lattice case endpoints')
        expected=complement_sets(span)
        require(actual['alternatives']==[dict(profile_id=i,possible_outcomes=v)for i,v in enumerate(expected)],'lattice prediction mismatch')
        comparisons+=5
    # These named profiles have disjoint predictions in the above-midpoint case.
    above=complement_sets(cases[2]);require(set(above[0]).isdisjoint(above[1]),'no exact rounding discriminator')
    links=document['pairs'];require(len(links)==4,'lattice pair count')
    for index,link in enumerate(links):
        negative,positive=extra[2*index:2*index+2]
        require(negative['id']==link['negative_request']and positive['id']==link['positive_request'],'lattice pair IDs')
        require(int(negative['args'][0],16)^(1<<63)==int(positive['args'][0],16),'lattice pair signs')
    magnitudes=[int(q['args'][0],16)for q in extra[1::2]]
    require(magnitudes[0]+1==magnitudes[1]and magnitudes[1]+1==magnitudes[2],'not adjacent binary64 controls')
    require(magnitudes[1]in document['mathematical_proposal']['root_bracket']and magnitudes[3]in document['mathematical_proposal']['robust_bracket'],'proposal bracket binding')
    for name in ('NORM-manifest.json','NORM-lattice-manifest.json'):
        manifest=json.loads((packet/name).read_text())
        for f,h in manifest['files_sha256'].items():require(digest(packet/f)==h,'frozen packet file changed')  # ubs:ignore — public artifact integrity, no secret comparison
    require(digest(Path(__file__))==SOURCE_HASH,'question checker changed during run')  # ubs:ignore — public artifact integrity, no secret comparison
    write_result(output,dict(complete=True,normal_questions=24,PMT_questions=64,new_excel_answers=0,
        independent_conditional_set_comparisons=comparisons,exact_above_midpoint_discriminator=True,
        complete_normal_programs=0,normal_checker_sha256=SOURCE_HASH,
        claim='Frozen requests contain no answers; all conditional output sets and open/closed probability cases checked independently. Mathematical proposal values are not treated as target helper outputs.',
        limitations='Profile/build and negative-anchor prerequisites remain unobserved on20527; the five survivors are reconstruction constraints, not full CDF programs.'))


def temme_check(trace, changes, output):
    """Exact replay of the final gamma correction, conditional on original leaves."""
    originals = list(records(trace))
    checked = []
    for row in records(changes):
        if row['program'] != 34:
            continue
        original = originals[row['trace_index']]
        require(all(row[k] == original[k] for k in ('id', 'expected')), 'Temme original binding')
        require(row['baseline'] == original['actual'], 'Temme original calculation binding')
        l, c, w, t, a, coefficient = map(binary64, row['temme_leaf_l_c_w_t_a_rt2pin'])
        root = sqrt_round(a, 'strict')
        def finish(model):
            correction = operation(operation(coefficient, t, 'mul', model), root, 'div', model)
            inner = operation(w, correction if l >= 1 else -correction, 'add', model)
            answer = B.round(operation(c, inner, 'mul', model))
            complement = operation(R(1, 2), operation(R(1, 2), -answer, 'add', model), 'add', model)
            p, q = (complement, answer) if l >= 1 else (answer, complement)
            return encoding(q if original['function'] in ('POISSON', 'POISSON.DIST') else p)
        require(finish('strict') == row['baseline'], 'original Temme arithmetic differs')
        require(finish('resident_pc64') == row['prediction'] == row['expected'], 'retained Temme arithmetic differs')
        checked.append(dict(id=row['id'], original=row['baseline'], retained=row['prediction']))
    require(checked, 'empty Temme proof')
    write_result(output, dict(complete=True, conditional_corrections=checked,
        inputs_sha256={str(p):digest(p) for p in (trace, changes)},
        normal_checker_sha256=SOURCE_HASH,
        limits='Independent rational PC64/PC53 rounding and exact squared-midpoint square root. Polynomial, EXP and scaled-erfc leaves are conditional original producer outputs; this does not identify their internal algorithms.'))


if __name__ == '__main__':
    import sys
    if sys.argv[1]=='linear':linear_check(*sys.argv[2:])
    elif sys.argv[1]=='symmetry':symmetry_check(*sys.argv[2:])
    elif sys.argv[1]=='questions':question_check(*sys.argv[2:])
    elif sys.argv[1]=='temme':temme_check(*sys.argv[2:])
    else:check(*sys.argv[1:])
