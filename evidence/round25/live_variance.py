"""Power-of-two guard questions through the repository's canonical Excel worker.

The finite constructions vary magnitude while preserving conditioning. Existing
answers select distinct operand shapes, never a numerical cutoff or coefficient.
"""
import hashlib
import json
import math
from pathlib import Path
import struct
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'checker'))
sys.path.insert(0, str(ROOT.parent / 'DnaCalc/tools/excel-service'))
from client import Client, ServiceError, number
from variance import evaluate, arithmetic
from exact import binary64
from normal import encoding
from fractions import Fraction as Q


def bits(x):
    return struct.pack('>d', x).hex()


def build(constraints):
    shapes = {}
    for row in constraints['centered_required']:
        if row['function'] != 'VAR' or any(isinstance(x, list) for x in row['args']):
            continue
        n = len(row['args'])
        if n in (2, 3, 4, 8):
            shapes.setdefault(n, row)
    result = []
    for function in ('VAR.S', 'VAR.P'):
        rows = []
        for n, source in sorted(shapes.items()):
            xs = [struct.unpack('>d', bytes.fromhex(x.removeprefix('0x')))[0]
                  for x in source['args']]
            for exponent in range(-12, 25, 2):
                values = [math.ldexp(x, exponent) for x in xs]
                exact = [Q(x) for x in values]
                predictions = {}
                for name, algorithm, mode in [('raw_pc64', 7 if function == 'VAR.S' else 8, 'stored_pc64'),
                                               ('centered_pc64', 0, 'stored_pc64'),
                                               ('centered_strict', 0, 'strict')]:
                    profile = dict(algorithm=algorithm, linear_order=0, square_order=0, models=[mode]*5)
                    predictions[name] = encoding(evaluate(profile, exact, function))
                if len(set(predictions.values())) < 2:
                    continue
                rows.append(dict(seed=source['id'], exponent=exponent, n=n,
                                 args=[bits(x) for x in values], predictions=predictions))
        result.append(dict(function=function, rows=rows))
    return result


def main():
    constraints = Path(sys.argv[1]).resolve()
    output = Path(sys.argv[2]).resolve()
    output.mkdir(exist_ok=False)
    sources = [Path(__file__), constraints, ROOT/'checker/variance.py', ROOT/'checker/exact.py',
               ROOT/'checker/normal.py', ROOT.parent/'DnaCalc/tools/excel-service/client.py']
    binding = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    questions = build(json.loads(constraints.read_text()))
    (output/'questions.json').write_text(json.dumps(questions, indent=2)+'\n')
    for packet in questions:
        function = packet['function']
        try:
            reply = Client(timeout=900).ask(function, [[number(x) for x in row['args']] for row in packet['rows']])
        except ServiceError as error:
            reply = error.response
        except TimeoutError as error:
            reply = dict(ok=False, error=str(error))
        (output/(function+'.json')).write_text(json.dumps(reply, indent=2)+'\n')
        print(function, 'calls', len(packet['rows']), 'ok', reply.get('ok'), flush=True)
    assert binding == {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (output/'binding.json').write_text(json.dumps(binding, indent=2)+'\n')


if __name__ == '__main__':
    main()
