"""Reusable questions for discrete numerical decisions (original MIT research).

Candidates supply predictions before answers arrive. Answers only eliminate
candidates; neither this module nor its census fits a formula to observations.
An adjacent bracket identifies a transition on the chosen path, not monotonicity
or a unique boundary elsewhere. Existing plateau_probe handles constant loads.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import json
import math


def identity(value):
    """Stable typed identity; callers should use bit strings for numeric answers."""
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


@dataclass
class Boundary:
    """Bisect an integer-indexed path whose endpoints decide differently.

    The caller evaluates path(index), asks its oracle and passes the discrete
    answer to observe. A path can interpolate a vector, step an input's binary64
    bits, or vary one structural parameter. Nonmonotone paths can contain other
    transitions: only this retained opposite-answer bracket is promised.
    """
    low: int
    high: int
    low_answer: object
    high_answer: object

    def __post_init__(self):
        if type(self.low) is not int or type(self.high) is not int or self.low >= self.high:
            raise ValueError('ordered integer endpoints required')
        if identity(self.low_answer) == identity(self.high_answer):
            raise ValueError('endpoints must decide differently')

    def question(self):
        return None if self.high - self.low == 1 else (self.low + self.high) // 2

    def observe(self, index, answer):
        if index != self.question() or index is None:
            raise ValueError('answer is not for the pending midpoint')
        key = identity(answer)
        if key == identity(self.low_answer):
            self.low = index
        elif key == identity(self.high_answer):
            self.high = index
        else:
            raise ValueError('third decision on a two-decision path')


def vector_path(left, right, index, last=1 << 52):
    """Fixed binary64 interpolation, with exact original endpoints retained."""
    if len(left) != len(right) or last < 1 or not 0 <= index <= last:
        raise ValueError('invalid vector path')
    if any(not math.isfinite(v) for v in [*left, *right]):
        raise ValueError('finite coordinates required')
    if index == 0:
        return list(left)
    if index == last:
        return list(right)
    t = index / last
    return [(1 - t) * a + t * b for a, b in zip(left, right)]


class DecisionPath:
    """Refine every observed change, retaining intermediate third decisions.

    Equal endpoints do not establish a constant interval. A coarse census must
    precede refinement; unobserved islands remain possible. Indices refer to
    the declared path, which need not be injective in binary64 coordinates.
    """
    def __init__(self, observations):
        self.observations = {}
        self.observe(observations)
        if len(self.observations) < 2 or any(type(i) is not int for i in self.observations):
            raise ValueError('at least two integer-indexed observations required')

    def changes(self):
        ordered = sorted(self.observations)
        return [(a, b) for a, b in zip(ordered, ordered[1:])
                if identity(self.observations[a]) != identity(self.observations[b])]

    def questions(self, points=64):
        if type(points) is not int or points < 1:
            raise ValueError('positive integer point count required')
        return sorted({a + (b-a)*k//(min(points, b-a-1)+1)
                       for a, b in self.changes() if b-a > 1
                       for k in range(1, min(points, b-a-1)+1)})

    def observe(self, observations):
        pending = dict(self.observations)
        for index, answer in observations:
            if type(index) is not int:
                raise ValueError('integer index required')
            if index in pending and identity(pending[index]) != identity(answer):
                raise ValueError('conflicting repeated observation')
            pending[index] = answer
        self.observations = pending

    def report(self):
        return dict(observations=sorted(self.observations.items()), changes=self.changes(),
                    adjacent=all(b-a == 1 for a,b in self.changes()),
                    limits='Only observed transitions on this path; equal endpoints do not exclude islands.')


def choose_split(predictions, candidates=None):
    """Minimise the largest surviving answer bucket, then sum of squared sizes.

    predictions maps candidate names to equal-length frozen answer lists. The
    stable input order breaks ties. Equivalent candidates remain represented;
    no separating question returns None, never an identification claim.
    """
    names = list(predictions) if candidates is None else list(candidates)
    if not names:
        raise ValueError('no candidates')
    lengths = {len(predictions[name]) for name in names}
    if len(lengths) != 1:
        raise ValueError('different question lists')
    best = None
    for i in range(next(iter(lengths))):
        counts = Counter(identity(predictions[name][i]) for name in names)
        if len(counts) < 2:
            continue
        score = (max(counts.values()), sum(n * n for n in counts.values()), i)
        if best is None or score < best[0]:
            best = (score, i, dict(counts))
    return None if best is None else dict(index=best[1], buckets=best[2], score=best[0][:2])


def eliminate(predictions, candidates, index, answer):
    key = identity(answer)
    return [name for name in candidates if identity(predictions[name][index]) == key]


class Census:
    """Map distinct outputs and their first questions, without a pass numerator."""
    def __init__(self):
        self.counts = Counter()
        self.first = {}
        self.observations = 0
        self.discoveries = []

    def observe(self, question, answer):
        key = identity(answer)
        self.observations += 1
        if key not in self.counts:
            self.first[key] = question
            self.discoveries.append(self.observations)
        self.counts[key] += 1

    def report(self):
        return dict(observations=self.observations, distinct=len(self.counts),
                    counts=dict(self.counts), first_questions=self.first,
                    discovery_positions=self.discoveries,
                    status='finite_observation_map')


def two_sum(a, b):
    """Error-free sum when intermediate values stay finite and do not underflow."""
    total = a + b
    part = total - a
    return total, (a - (total - part)) + (b - part)


def product_residue(a, b):
    """Exact rational residue for midpoint question design, including underflow."""
    product = a * b
    if not all(math.isfinite(v) for v in (a, b, product)):
        raise ValueError('finite product required')
    return product, Fraction(a) * Fraction(b) - Fraction(product)


def main():
    """File interface: frozen splitting, streaming census, or one path question.

    split INPUT OUTPUT: INPUT has predictions and optional questions.
    census INPUT.jsonl OUTPUT: every record has question and answer.
    boundary INPUT OUTPUT: INPUT has Boundary's four fields, optional sequential
    observations [{index,answer}], and optional left/right/last vector path.
    OUTPUT is exclusively created and binds the exact input/source bytes.
    """
    import argparse
    import hashlib
    from pathlib import Path
    parser=argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('mode',choices=['split','census','boundary'])
    parser.add_argument('input',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();raw=args.input.read_bytes()
    if args.mode=='census':
        census=Census()
        for line in raw.decode().splitlines():
            record=json.loads(line);census.observe(record['question'],record['answer'])
        result=census.report()
    else:
        document=json.loads(raw)
        if args.mode=='split':
            result=choose_split(document['predictions'],document.get('candidates'))
            if result is not None and 'questions'in document:result['question']=document['questions'][result['index']]
            result=dict(selection=result,status='undetermined_within_question_pool'if result is None else'separating_question')
        else:
            boundary=Boundary(document['low'],document['high'],document['low_answer'],document['high_answer'])
            for record in document.get('observations',[]):boundary.observe(record['index'],record['answer'])
            index=boundary.question()
            result=dict(low=boundary.low,high=boundary.high,low_answer=boundary.low_answer,
                        high_answer=boundary.high_answer,next_index=index,
                        status='adjacent_indices'if index is None else'next_question')
            if index is not None and 'left'in document:
                result['question']=vector_path(document['left'],document['right'],index,document.get('last',1<<52))
    result['input_sha256']=hashlib.sha256(raw).hexdigest()
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x')as stream:json.dump(result,stream,indent=2,allow_nan=False)


if __name__=='__main__':main()
