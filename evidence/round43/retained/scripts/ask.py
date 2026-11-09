"""Ask the frozen R42 questions through the question service; keep every raw reply and every wait.

usage: python3 ask.py <questions-dir> <out-dir> [normal|mapping|regression|all]
Questions and both predictions were frozen before this script ran (sha256 in freeze.json).
"""
import json, os, sys, time, collections
sys.path.insert(0, os.path.expanduser('~/projects/DnaCalc/DnaCalc/tools/excel-service'))
from client import Client, number, number_array, boolean

QD, OUT = sys.argv[1], sys.argv[2]
WHAT = sys.argv[3] if len(sys.argv) > 3 else 'all'
QFILES = sys.argv[4].split(',') if len(sys.argv) > 4 else ['LINEST.questions.json', 'LOGEST.questions.json', 'LINEST.sign-questions.json', 'LINEST.zero-questions.json', 'LINEST.pivot-questions.json']
os.makedirs(OUT, exist_ok=True)
excel = Client(timeout=600)
waits = open(os.path.join(OUT, 'waits.tsv'), 'a')

def timed(function, rows, **kw):
    t = time.time()
    try:
        reply = excel.ask(function, rows, priority='bulk', **kw)
        err = ''
    except Exception as e:
        reply = None; err = '%s: %s' % (type(e).__name__, e)
    dt = time.time() - t
    if dt > 3 or err:
        waits.write('%s\t%s\t%d\t%.3f\t%s\n' % (time.strftime('%H:%M:%S', time.gmtime()), function, len(rows), dt, err))
        waits.flush()
    return reply, err

if WHAT in ('normal', 'all'):
    q = json.load(open(os.path.join(QD, 'NORM.S.INV.lower.json')))
    rows = [[number(p)] for p in q['points']]
    replies = []
    for k in range(0, len(rows), 1000):
        r, err = timed('NORM.S.INV', rows[k:k + 1000])
        replies.append(dict(start=k, reply=r, error=err))
    json.dump(dict(questions_sha='see freeze.json', replies=replies), open(os.path.join(OUT, 'NORM.S.INV.lower.replies.json'), 'w'))

if WHAT in ('mapping', 'all'):
    q = json.load(open(os.path.join(QD, 'NORM.S.INV.mapping.json')))
    replies = []
    for w in q['windows']:
        s0 = int(w['start'], 16)
        r, err = timed('NORM.S.INV', [[number('%016x' % (s0 + j))] for j in range(w['count'])])
        replies.append(dict(window=w, reply=r, error=err))
    json.dump(dict(replies=replies), open(os.path.join(OUT, 'NORM.S.INV.mapping.replies.json'), 'w'))

if WHAT in ('regression', 'all'):
    for qf in QFILES:
        q = json.load(open(os.path.join(QD, qf))); fn = q['function']
        groups = collections.defaultdict(list)
        for w in q['witnesses']:
            groups[tuple(w['result_index'])].append(w)
        out = []
        for idx, ws in sorted(groups.items()):
            rows = []
            for w in ws:
                y, x, const, stats = w['args']
                rows.append([number_array([[r[0]] for r in y]), number_array(x), boolean(const), boolean(stats)])
            for k in range(0, len(rows), 100):
                r, err = timed(fn, rows[k:k + 100], result_index=idx)
                out.append(dict(result_index=list(idx), ids=[w['id'] for w in ws[k:k + 100]], reply=r, error=err))
        json.dump(dict(questions=qf, replies=out), open(os.path.join(OUT, qf.replace('.json', '.replies.json')), 'w'))
print('done')
