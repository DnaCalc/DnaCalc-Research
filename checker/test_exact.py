"""Deliberate independent arithmetic suite; these are not Excel readings."""
from fractions import Fraction as R
import itertools
import math
import random
import struct
import unittest
import copy
import contextlib
import io
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

from exact import (BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P53,
                   Lattice, Span, binary64, inverse, operation, two)
from check import validate_row, quotient, product_preimage, write_result


class ExactArithmetic(unittest.TestCase):
    def test_exhaustive_small_lattices(self):
        checks = 0
        for fmt in (Lattice(2, -2, 3), Lattice(3, -3, 4), Lattice(4, -2, 2)):
            table = [(k, fmt.value(k)) for k in range(-fmt.last-1, fmt.last+2)]
            for (k, x), (n, y) in zip(table, table[1:]):
                self.assertEqual(n, k+1)
                self.assertLess(x, y)
                midpoint = (x+y)/2
                for t in (x, midpoint-(y-x)/16, midpoint, midpoint+(y-x)/16, y):
                    # Brute finite table with the two ideal overflow neighbors.
                    distance = min(abs(t-v) for _, v in table)
                    nearest = [(i, v) for i, v in table if abs(t-v) == distance]
                    index, result = min(nearest, key=lambda iv: iv[0] % 2)
                    if abs(index) > fmt.last:
                        result = None
                    self.assertEqual(fmt.round(t), result)
                    if abs(k) <= fmt.last and x <= t < y:
                        self.assertEqual(fmt.floor_coordinate(t), k)
                    checks += 1
            # Independently enumerate lattice members and their midpoint cells.
            probes = [v for _, v in table]
            probes += [(a[1]+b[1])/2 for a, b in zip(table, table[1:])]
            for a, b in zip(probes[::3], probes[2::3]):
                lo, hi = sorted((a, b))
                if lo == hi:
                    continue
                for lc, hc in itertools.product((False, True), repeat=2):
                    span = Span(lo, hi, lc, hc)
                    members = [(k, v) for k, v in table if abs(k) <= fmt.last and span.contains(v)]
                    bounds = fmt.coordinates(span)
                    self.assertEqual(bounds, (members[0][0], members[-1][0]) if members else None)
                    pre = fmt.pullback(span)
                    if members:
                        first, last = members[0][0], members[-1][0]
                        self.assertEqual(pre, Span((fmt.value(first-1)+fmt.value(first))/2,
                            (fmt.value(last)+fmt.value(last+1))/2, first % 2 == 0, last % 2 == 0))
                    else:
                        self.assertIsNone(pre)
                    checks += 1
        print('exhaustive small-lattice checks:', checks)

    def test_binary64_against_raw_encoding_and_neighbors(self):
        rng = random.Random(0x524e395f43484543)
        values = [0, 1, 2, (1 << 52)-1, 1 << 52, 0x3ff0000000000000,
                  0x7fefffffffffffff] + [rng.randrange(0x7ff0000000000000) for _ in range(2048)]
        for bits in values:
            for signed in (bits, bits | (1 << 63)):
                x = struct.unpack('>d', signed.to_bytes(8, 'big'))[0]
                q = binary64(signed)
                self.assertEqual(q, R.from_float(x))
                self.assertEqual(B.round(q), q)
                k = B.floor_coordinate(q)
                self.assertEqual(B.value(k), q)
                for direction, delta in ((math.inf, 1), (-math.inf, -1)):
                    neighbor = math.nextafter(x, direction)
                    if math.isfinite(neighbor):
                        self.assertEqual(B.value(k+delta), R.from_float(neighbor))

    def test_deliberate_double_rounding_and_mixed_precision(self):
        one = Span(R(1), R(1))
        self.assertEqual(E.pullback(B.pullback(one)),
            Span(1-two(-54)-two(-65), 1+two(-53)+two(-64)))
        b = two(-53)+two(-65)
        self.assertEqual(B.round(1+b), 1+two(-52))
        self.assertEqual(B.round(E.round(1+b)), 1)
        resident = E.round(1+two(-60))
        self.assertEqual(operation(resident, R(-1), 'add', 'strict'), 0)
        self.assertEqual(operation(resident, R(-1), 'add', 'resident_pc64'), two(-60))
        self.assertEqual(operation(resident, R(-1), 'add', 'resident_pc53'), two(-60))
        wide = P53.round(two(1023)+two(1023))
        self.assertEqual(wide, two(1024))
        self.assertEqual(P53.round(wide/2), two(1023))
        self.assertIsNone(B.round(two(1024)))
        self.assertEqual(E.round(resident)-1, two(-60))
        self.assertEqual(B.round(resident)-1, 0)

    def test_open_endpoints_negative_inversion_underflow_overflow(self):
        odd = 1+two(-52)
        span = B.pullback(Span(odd, odd))
        self.assertFalse(span.left)
        self.assertFalse(span.right)
        self.assertEqual(B.round(span.low), R(1))
        self.assertEqual(B.round(span.high), 1+two(-51))
        reflected = inverse(Span(R(1), R(1)), R(-3), 'mul', True, 'resident_pc64')
        self.assertEqual(reflected, E.pullback(Span(R(1), R(1))).linear(R(-1, 3)))
        self.assertIsNone(B.coordinates(Span(R(1), 1+two(-52), False, False)))
        self.assertEqual(B.round(two(-1075)), 0)
        self.assertEqual(B.round(3*two(-1075)), two(-1073))
        self.assertEqual(B.round(-two(-1075)), 0)  # rational zero; sign is out of scope
        maximum = B.value(B.last)
        threshold = (maximum+two(1024))/2
        self.assertIsNone(B.round(threshold))
        self.assertEqual(B.round(threshold-two(970)), maximum)
        self.assertEqual(B.coordinates(Span(-two(1100), two(1100))), (-B.last, B.last))
        with self.assertRaises(ValueError):
            binary64(0x7ff0000000000000)
        with self.assertRaises(ValueError):
            Span(R(-1), R(1)).reciprocal(R(1))

    def test_original_argument_binding_and_defaults(self):
        args = ['0x3fc0000000000000', '0x4000000000000000', '0x3ff0000000000000',
                '0x8000000000000000', '0x4000000000000000']
        row = dict(raw=dict(function='PMT', result_index=None, args=args),
                   leaves=[int(args[0], 16), int(args[2], 16), int(args[3], 16),
                           0x3ff0000000000000, 0xbfd0000000000000], neg_log=0xbfd0000000000000)
        validate_row(row)
        for index in range(4):
            altered = copy.deepcopy(row)
            altered['leaves'][index] ^= 1
            with self.assertRaises(ValueError):
                validate_row(altered)
        for arity in (3, 4):
            altered = copy.deepcopy(row)
            altered['raw']['args'] = args[:arity]
            altered['leaves'][3] = 0
            if arity == 3:
                altered['leaves'][2] = 0
            validate_row(altered)
        altered = copy.deepcopy(row)
        altered['raw']['result_index'] = [1, 1]
        with self.assertRaises(ValueError):
            validate_row(altered)

    def test_quotient_trace_expected_then_actual_columns(self):
        # Constructed arithmetic control: 1.5+oneULP cannot equal DR(q*1.5)
        # for binary64 q. This is not a recorded Excel observation.
        rate_bits, expected_bits = '3ff8000000000000', '3ff8000000000001'
        rate, expected = binary64(rate_bits), binary64(expected_bits)
        span = product_preimage(expected).linear(1/rate)
        self.assertIsNone(B.coordinates(span))
        cert = dict(id='synthetic-parser-gap', rate_bits=rate_bits, expected_bits=expected_bits,
            old_actual_bits=rate_bits, old_quotient_bits='3ff0000000000000',
            lower=span.json()['lower'], upper=span.json()['upper'], closed=span.left,
            neighbor_magnitudes=['3ff0000000000000', '3ff0000000000001'])
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            (base/'raw.tsv').write_text('Q\tsynthetic-parser-gap\t'+rate_bits+
                '\t3ff0000000000000\t'+expected_bits+'\t'+rate_bits+'\n')
            (base/'certs.jsonl').write_text(json.dumps(cert)+'\n')
            with contextlib.redirect_stdout(io.StringIO()):
                quotient(base/'raw.tsv', base/'certs.jsonl', base/'result.json')
            result = json.loads((base/'result.json').read_text())
            self.assertEqual(result['counts']['unreachable_by_binary64_quotient'], 1)
            self.assertEqual(result['counts']['original_quotient_product_replays'], 1)

    def test_source_drift_refuses_a_success_receipt(self):
        with tempfile.TemporaryDirectory() as temp:
            output = Path(temp)/'receipt.json'
            with patch('check.LOADED_SOURCE_SHA256', {'exact.py': 'changed'}):
                with self.assertRaisesRegex(ValueError, 'source changed'):
                    write_result(output, {})
            self.assertFalse(output.exists())


class NormalReferenceTests(unittest.TestCase):
    def test_complement_lattice_discriminator(self):
        from normal import complement_sets
        p,w=two(-54),two(-66)
        exactly=complement_sets(Span(p,p))
        self.assertEqual(exactly[:3],[['0x3ff0000000000000']]*3)
        self.assertEqual(exactly[3:],[['0x3fefffffffffffff','0x3ff0000000000000']]*2)
        above=complement_sets(Span(p,p+w,False,True))
        self.assertEqual(above,[['0x3fefffffffffffff'],['0x3ff0000000000000'],['0x3ff0000000000000'],['0x3fefffffffffffff'],['0x3fefffffffffffff']])

    def test_root_midpoints_and_extended_root_store(self):
        from normal import sqrt_round, encoding
        middle = R(1)+two(-53)
        self.assertEqual(sqrt_round(middle*middle, 'resident_pc64'), middle)
        self.assertEqual(sqrt_round(middle*middle, 'stored_pc64'), R(1))
        for model in ('strict','resident_pc64','resident_pc53'):
            root = sqrt_round(R(2), model)
            fmt = B if model == 'strict' else E if model.endswith('pc64') else P53
            k = fmt.floor_coordinate(root)
            lower, upper = (fmt.value(k-1)+root)/2, (root+fmt.value(k+1))/2
            self.assertLess(lower*lower,R(2))
            self.assertGreater(upper*upper,R(2))
        self.assertEqual(encoding(two(-1023)), '0x0000000000000000')
        self.assertEqual(encoding(two(-1022)), '0x0010000000000000')

    def test_concrete_transcendental_bindings_reject_altered_input_or_result(self):
        from normal import Replay
        nodes = [dict(op='input',bits='0000000000000000',value='binary64:0000000000000000'),
                 dict(op='exp',args=[0],provider=3,concrete_leaf='3:binary64:0000000000000000',value='binary64:3ff0000000000000')]
        helper = {'3:binary64:0000000000000000':dict(provider=3,input='binary64:0000000000000000',output='binary64:3ff0000000000000')}
        m = Replay(nodes,{},helper);m.node('input',R(0),bits='0000000000000000');m.exp(0,3)
        for key,value in [('input','binary64:0000000000000001'),('output','binary64:3ff0000000000001'),('provider',2)]:
            altered = copy.deepcopy(helper);altered[next(iter(altered))][key] = value
            m = Replay(nodes,{},altered);m.node('input',R(0),bits='0000000000000000')
            with self.assertRaises(ValueError):m.exp(0,3)

    def test_directed_mathematical_reference_uses_two_precisions(self):
        from reference import normal_bounds
        self.assertEqual(normal_bounds('0000000000000000')['rounded_binary64'],'0x3fe0000000000000')
        for bits in ('3ff0000000000000','bff0000000000000','bfe537a6f4de9bd4','c042c00000000000'):
            a,b = normal_bounds(bits,192),normal_bounds(bits,384)
            self.assertEqual(a['rounded_binary64'],b['rounded_binary64'])
            self.assertLessEqual(R(*map(int,a['lower'])),R(*map(int,b['lower'])))
            self.assertGreaterEqual(R(*map(int,a['upper'])),R(*map(int,b['upper'])))


if __name__ == '__main__':
    unittest.main()
