import unittest
from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from racer.question_designer import (Boundary, vector_path, choose_split,
                                     eliminate, Census, two_sum, product_residue, DecisionPath)


class Questions(unittest.TestCase):
    def test_path_preserves_third_decisions_and_all_seen_changes(self):
        path = DecisionPath([(0, 'a'), (64, 'b')])
        indices = path.questions(7)
        self.assertEqual(indices, [8,16,24,32,40,48,56])
        path.observe((i, 'a' if i < 17 else 'c' if i < 43 else 'b') for i in indices)
        while indices := path.questions(64):
            path.observe((i, 'a' if i < 17 else 'c' if i < 43 else 'b') for i in indices)
        self.assertEqual(path.changes(), [(16,17),(42,43)])
        self.assertTrue(path.report()['adjacent'])
        with self.assertRaises(ValueError): path.observe([(0,'b')])
        before=dict(path.observations)
        with self.assertRaises(ValueError):path.observe([(65,'c'),(0,'b')])
        self.assertEqual(path.observations,before)
        with self.assertRaises(ValueError):DecisionPath([(0,'a'),(0,'b'),(1,'c')])

    def test_equal_path_endpoints_make_no_island_claim(self):
        path = DecisionPath([(0, 'a'), (100, 'a')])
        self.assertEqual(path.questions(), [])
        self.assertIn('islands', path.report()['limits'])

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

    def test_cli_boundary_emits_the_next_vector(self):
        with tempfile.TemporaryDirectory()as directory:
            source=Path(directory)/'state.json';target=Path(directory)/'question.json'
            source.write_text(json.dumps(dict(low=0,high=8,low_answer=False,high_answer=True,
                observations=[dict(index=4,answer=True)],left=[0.],right=[8.],last=8)))
            subprocess.run([sys.executable,'-m','racer.question_designer','boundary',str(source),str(target)],
                           cwd=Path(__file__).resolve().parents[1],check=True,capture_output=True,timeout=30)
            result=json.loads(target.read_text())
            self.assertEqual((result['low'],result['high'],result['next_index']),(0,4,2))
            self.assertEqual(result['question'],[2.])

    def test_cli_census_is_typed_and_never_overwrites_evidence(self):
        with tempfile.TemporaryDirectory()as directory:
            source=Path(directory)/'answers.jsonl';target=Path(directory)/'census.json'
            source.write_text(''.join(json.dumps(dict(question=i,answer=a))+'\n'for i,a in enumerate([True,1,True])))
            command=[sys.executable,'-m','racer.question_designer','census',str(source),str(target)]
            subprocess.run(command,cwd=Path(__file__).resolve().parents[1],check=True,capture_output=True,timeout=30)
            original=target.read_bytes();result=json.loads(original)
            self.assertEqual(result['distinct'],2)
            self.assertEqual(result['discovery_positions'],[1,2])
            again=subprocess.run(command,cwd=Path(__file__).resolve().parents[1],check=False,capture_output=True,timeout=30)
            self.assertNotEqual(again.returncode,0)
            self.assertEqual(target.read_bytes(),original)


if __name__ == '__main__': unittest.main()
