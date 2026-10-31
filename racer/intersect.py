"""Exact shared-rate/nper expm1 constraints for the R7 FV=0 grammar.

Catalogues describe research programs, never answers. Raw Excel witnesses come
only from the existing loader's frozen trace. Environment remains part of keys.
"""
from fractions import Fraction as Q
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse
import collections
import hashlib
import json
import sys
import time
from intervals import B64, E64, Interval, inverse_op, number, stored

MODELS = dict(strict=0, stored_pc64=1, resident_pc64=2)


def rnd(x, model):
    return B64.rn(x) if model == 0 else stored(x) if model == 1 else E64.rn(x)


def compile_expr(expr, models, site):
    if 'leaf' in expr:
        return expr['leaf'], site
    op = next(iter(expr))
    a, site = compile_expr(expr[op][0], models, site)
    b, site = compile_expr(expr[op][1], models, site)
    return (op, a, b, models[site]), site + 1


def has_em(expr):
    return expr == 2 if isinstance(expr, int) else has_em(expr[1]) or has_em(expr[2])


def const_value(expr, args):
    if isinstance(expr, int):
        assert expr != 2
        return args[expr]
    op, a, b, model = expr
    a, b = const_value(a, args), const_value(b, args)
    if model == 0:
        a, b = B64.rn(a), B64.rn(b)
    if a is None or b is None or (op == 'div' and b == 0):
        raise ArithmeticError('nonfinite/zero constant subtree')
    result = rnd(a*b if op == 'mul' else a/b, model)
    if result is None:
        raise ArithmeticError('nonfinite constant subtree')
    return result


def inverse_expr(expr, args, interval):
    if interval is None:
        return None
    if isinstance(expr, int):
        assert expr == 2
        return interval
    op, a, b, model = expr
    left = has_em(a)
    assert left != has_em(b), 'exactly one occurrence of em'
    c = const_value(b if left else a, args)
    interval = inverse_op(interval, model, op, c, left)
    return inverse_expr(a if left else b, args, interval)


def row_preimage(row, profile):
    c, expr, models = profile
    rate, pv, _, timing, _ = [Q(number(f'{b:016x}')) for b in row['leaves']]
    tf = rnd(1 + rate*timing, models[0])
    target = row['raw']['expected']
    if not target.startswith('0x'):
        raise ArithmeticError('nonnumeric expected outcome')
    value = number(target)
    if value == 0 or not (-float('inf') < value < float('inf')):
        raise ArithmeticError('nonfinite/zero target')
    if pv == 0 or tf is None or tf == 0:
        raise ArithmeticError('zero/nonfinite constant')
    target = B64.preimage(Interval(Q(value), Q(value)))
    return inverse_expr(expr, [pv, rate, None, tf], target)


def key(row):
    raw = row['raw']
    return (row['leaves'][0], raw['args'][1], json.dumps(raw['reading'], sort_keys=True))


ROWS = None
PROFILES = None


def process(item):
    gid, ids = item
    records = []
    counts = collections.Counter()
    for profile in PROFILES:
        c = profile[0]
        current = None
        lower_id = upper_id = None
        unsupported = 0
        visited = 0
        contradiction = None
        for idx in ids:
            row = ROWS[idx]
            try:
                preimage = row_preimage(row, profile)
            except (ArithmeticError, AssertionError) as e:
                unsupported += 1
                counts['unsupported_row_profile_constraints'] += 1
                continue
            visited += 1
            if preimage is None:
                contradiction = dict(new_row=idx, previous_interval=None,
                                     new_preimage=None, bound_rows=[])
                current = None
                break
            if current is None:
                current = preimage
                lower_id = upper_id = idx
                continue
            new = current.intersect(preimage)
            if new is None:
                contradiction = dict(new_row=idx, previous_interval=current.json(),
                                     new_preimage=preimage.json(),
                                     bound_rows=sorted(set([lower_id, upper_id])))
                current = None
                break
            if new.lo != current.lo or new.lc != current.lc:
                lower_id = idx
            if new.hi != current.hi or new.hc != current.hc:
                upper_id = idx
            current = new
        status = 'empty' if contradiction else 'nonempty' if current is not None else 'unsupported'
        counts[status + '_group_profiles'] += 1
        b64 = current is not None and B64.preimage(current) is not None
        e64 = current is not None and E64.preimage(current) is not None
        counts['group_profiles_with_binary64_em'] += b64
        counts['group_profiles_with_ext80_em'] += e64
        records.append(dict(group=gid, combine=c['id'], constraints_visited=visited,
                            unsupported_constraints=unsupported, status=status,
                            interval=current.json() if current is not None else None,
                            binary64_value_exists=b64, ext80_value_exists=e64,
                            bound_rows=sorted(set([lower_id,upper_id])) if current else [],
                            first_empty=contradiction))
    return gid, records, dict(counts)


def main():
    global ROWS, PROFILES
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('catalogue')
    parser.add_argument('output')
    parser.add_argument('--profiles', choices=['baseline', 'all', 'uniform'], default='all')
    parser.add_argument('--groups', choices=['all', 'discovery'], default='discovery')
    parser.add_argument('--group-limit', type=int, default=32)
    parser.add_argument('--workers', type=int, default=8)
    a = parser.parse_args()
    start = time.monotonic()
    # Original complete records stay in the frozen input. Fork only fields
    # needed by rational arithmetic; this keeps the eight-worker budget small.
    ROWS = []
    environments = {}
    for line in open(a.input):
        r = json.loads(line)
        if number(f'{r["leaves"][2]:016x}') != 0:
            continue
        raw = r['raw']
        env = json.dumps(raw['reading'], sort_keys=True)
        environment = environments.setdefault(env, raw['reading'])
        ROWS.append(dict(leaves=r['leaves'], trace_index=r['trace_index'],
            raw=dict(args=[None, raw['args'][1]], expected=raw['expected'], reading=environment)))
    groups = collections.defaultdict(list)
    for idx, r in enumerate(ROWS):
        groups[key(r)].append(idx)
    group_list = list(groups.values())
    if a.groups == 'discovery':
        group_list = group_list[:a.group_limit]
    cat = Path(a.catalogue)
    expressions = json.loads((cat/'expressions.json').read_text())
    combines = json.loads((cat/'combines.json').read_text())
    if a.profiles == 'baseline':
        combines = [c for c in combines if c['id'] == 36]
    if a.profiles == 'uniform':
        combines = [c for c in combines if len(set(c['models'])) == 1]
    PROFILES = []
    for c in combines:
        models = [MODELS[m] for m in c['models']]
        expr, site = compile_expr(expressions[c['expression']], models, 1)
        assert site == 4
        PROFILES.append((c, expr, models))
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    with (out/'groups.jsonl').open('x') as f:
        for gid, ids in enumerate(group_list):
            r = ROWS[ids[0]]
            f.write(json.dumps(dict(group=gid, rate=f'0x{r["leaves"][0]:016x}',
                nper=r['raw']['args'][1], environment=r['raw']['reading'], rows=ids,
                trace_indices=[ROWS[i]['trace_index'] for i in ids]))+'\n')
    total = collections.Counter()
    completed = 0
    with (out/'intersections.jsonl').open('x') as f, ProcessPoolExecutor(a.workers) as executor:
        for start_group in range(0, len(group_list), 2*a.workers):
            # Bound submitted/result groups; one slow group cannot accumulate
            # the results of the rest of the entire evidence bank in memory.
            batch = list(enumerate(group_list))[start_group:start_group+2*a.workers]
            for gid, records, counts in executor.map(process, batch, chunksize=1):
                for r in records:
                    f.write(json.dumps(r,separators=(',',':'))+'\n')
                total.update(counts)
                completed += 1
            f.flush()
            checkpoint = dict(groups_done=completed, total_groups=len(group_list),
                              counts=total, seconds=time.monotonic()-start)
            (out/'checkpoint.json').write_text(json.dumps(checkpoint, indent=2)+'\n')
            print(json.dumps(checkpoint), flush=True)
    result = dict(complete=True, group_selection=a.groups, groups=len(group_list),
        profiles=len(PROFILES), counts=total, seconds=time.monotonic()-start,
        arithmetic='exact rational preimages of every store/operation; endpoint parity included',
        scope='FV=0 finite discount; one common arbitrary em per rate/nper/environment and combine profile. Empty intersection rules out that shared helper/profile, not all programs.',
        unsupported_policy='unsupported constraints are counted and supply no conclusion; an empty intersection remains a contradiction of the supported subset',
        input_sha256=hashlib.sha256(Path(a.input).read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    with (out/'summary.json').open('x') as f:
        json.dump(result,f,indent=2)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
