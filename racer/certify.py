"""Check exact RN53 rounding intervals against every enumerated final value.

Inputs are racer's research outputs, never additional Excel answers. Fraction
arithmetic decodes every bit of a resident Ext80 value, including underflows.
"""
from fractions import Fraction
from pathlib import Path
import collections
import hashlib
import json
import math
import struct
import sys
from intervals import B64


def number(bits):
    return struct.unpack(">d", int(bits, 16).to_bytes(8, "big"))[0]


def power2(exponent):
    return Fraction(1 << exponent) if exponent >= 0 else Fraction(1, 1 << -exponent)


def exact_root(text):
    kind, bits = text.split(":")
    if kind == "binary64":
        value = number(bits)
        return Fraction(value) if math.isfinite(value) else None
    assert kind == "ext80-le"
    raw = bytes.fromhex(bits)
    assert len(raw) == 10
    significand = int.from_bytes(raw[:8], "little")
    se = int.from_bytes(raw[8:], "little")
    exponent = se & 0x7fff
    if exponent == 0x7fff:
        return None
    sign = -1 if se & 0x8000 else 1
    return sign * significand * power2((exponent or 1) - 16383 - 63)


def interval(expected):
    assert expected.startswith("0x")
    y = number(expected[2:])
    assert math.isfinite(y) and y != 0 and abs(y) >= sys.float_info.min
    exact = Fraction(y)
    prev = math.nextafter(y, -math.inf)
    nxt = math.nextafter(y, math.inf)
    prev = Fraction(prev) if math.isfinite(prev) else -power2(1024)
    nxt = Fraction(nxt) if math.isfinite(nxt) else power2(1024)
    return (prev + exact) / 2, (exact + nxt) / 2, int(expected, 16) % 2 == 0


def pair(q):
    return [str(q.numerator), str(q.denominator)]


def publication(text):
    q = exact_root(text)
    if q is None:
        return 'error:Num'
    rounded = B64.rn(q)
    if rounded is None:
        return 'error:Num'
    if rounded and abs(rounded) < Fraction(sys.float_info.min):
        return 'fault:NotImplemented'
    bits = struct.unpack('>Q', struct.pack('>d', float(rounded)))[0] if rounded else 0
    return f'0x{bits:016x}'


def certify(directory):
    root = Path(directory)
    counts = collections.Counter()
    records = []
    for path in sorted(root.glob("unreachable-roots-*.jsonl")):
        for line in path.open():
            record = json.loads(line)
            expected = record['expected']
            target = number(expected[2:]) if expected.startswith('0x') else None
            normal = target is not None and math.isfinite(target) and abs(target) >= sys.float_info.min
            if not normal:
                outputs = sorted(set(publication(text) for text in record['raw_roots']))
                assert expected not in outputs, (record['order'], expected, outputs)
                counts['certified_unreachable_answers'] += 1
                counts['exact_roots_checked'] += len(record['raw_roots'])
                counts['typed_or_zero_publication_certificates'] += 1
                records.append(dict(order=record['order'],trace_index=record['trace_index'],
                    expected=expected,exact_publications=outputs,
                    reason='Exact rational RN53 followed by application publication kind/zero policy',
                    roots_sha256=hashlib.sha256(json.dumps(record['raw_roots'],separators=(',',':')).encode()).hexdigest()))
                continue
            low, high, closed = interval(expected)
            below = []
            above = []
            nonfinite = 0
            for text in record["raw_roots"]:
                q = exact_root(text)
                if q is None:
                    nonfinite += 1
                    continue
                assert not (low <= q <= high if closed else low < q < high), (record["order"], text)
                if q < low or (q == low and not closed):
                    below.append(q)
                else:
                    assert q > high or (q == high and not closed)
                    above.append(q)
            counts["certified_unreachable_answers"] += 1
            counts["exact_roots_checked"] += len(record["raw_roots"])
            counts["nonfinite_roots"] += nonfinite
            records.append({
                "order": record["order"], "trace_index": record["trace_index"],
                "expected": record["expected"], "rn53_interval": {
                    "lower": pair(low), "upper": pair(high), "closed": closed,
                },
                "roots_below": len(below), "roots_above": len(above),
                "nonfinite_roots": nonfinite,
                "largest_root_below": pair(max(below)) if below else None,
                "smallest_root_above": pair(min(above)) if above else None,
                "roots_sha256": hashlib.sha256(json.dumps(record["raw_roots"], separators=(",", ":")).encode()).hexdigest(),
            })
    summary = json.loads((root / "summary.json").read_text())
    assert len(records) == summary["whole_family_unreachable"]
    with (root / "rounding-interval-certificates.jsonl").open("w") as f:
        for r in sorted(records, key=lambda r: r["order"]):
            f.write(json.dumps(r, separators=(",", ":")) + "\n")
    result = {"counts": dict(counts), "arithmetic": "exact rational; no approximate inverse or join window",
              "claim": f"Every enumerated final tagged value is outside the recorded target's RN53 preimage. The enumeration covers precisely the frozen {summary['programs']}-program family, not all possible programs.",
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (root / "rounding-interval-summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    certify(sys.argv[1])
