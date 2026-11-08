import unittest
from fractions import Fraction
from racer.question_designer import (Boundary, vector_path, choose_split,
                                     eliminate, Census, two_sum, product_residue)


class Questions(unittest.TestCase):
    def test_boundary_retains_opposite_endpoints(self):
        edge = Boundary(0, 1 << 60, False, True)
        count = 0
        while (mid := edge.question()) is not None:
            edge.observe(mid, mid >= 12345678901234567)
            count += 1
        self.assertEqual((edge.low, edge.high), (12345678901234566, 12345678901234567))
        self.assertEqual(count, 60)

    def test_boundary_rejects_wrong_or_third_answers(self):
        with self.assertRaises(ValueError): Boundary(1, 2, 0, 0)
        edge = Boundary(0, 8, 'a', 'b')
        with self.assertRaises(ValueError): edge.observe(3, 'a')
        with self.assertRaises(ValueError): edge.observe(4, 'c')
        self.assertEqual((edge.low, edge.high), (0, 8))

    def test_vector_endpoints_and_cancellation(self):
        self.assertEqual(vector_path([1., -2.], [3., 4.], 0, 4), [1., -2.])
        self.assertEqual(vector_path([1., -2.], [3., 4.], 4, 4), [3., 4.])
        self.assertEqual(vector_path([1., -2.], [3., 4.], 2, 4), [2., 1.])

    def test_balanced_question_and_elimination(self):
        predictions = dict(a=[0, 0, 0], b=[0, 1, 0], c=[0, 2, 1], d=[1, 2, 1])
        self.assertEqual(choose_split(predictions)['index'], 1)
        self.assertEqual(eliminate(predictions, list(predictions), 2, 1), ['c', 'd'])
        self.assertIsNone(choose_split({'a':[0], 'b':[0]}))
        self.assertEqual(eliminate(predictions, ['a'], 0, 9), [])

    def test_census_preserves_types_and_first_question(self):
        census = Census()
        for q, a in enumerate([True, 1, True, '1']): census.observe(q, a)
        report = census.report()
        self.assertEqual(report['distinct'], 3)
        self.assertEqual(report['discovery_positions'], [1, 2, 4])
        self.assertEqual(report['first_questions']['true'], 0)

    def test_instrumentation_is_exact(self):
        a, b = 1., 2. ** -54
        total, error = two_sum(a, b)
        self.assertEqual(Fraction(total) + Fraction(error), Fraction(a) + Fraction(b))
        a, b = 1. + 2. ** -52, 1. - 2. ** -52
        product, error = product_residue(a, b)
        self.assertEqual(Fraction(product) + error, Fraction(a) * Fraction(b))


if __name__ == '__main__': unittest.main()
