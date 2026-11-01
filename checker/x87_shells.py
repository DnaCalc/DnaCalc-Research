"""Exact operand-lattice reachability of the original-rate final x87 product.

This is a necessary structural diagnostic, never a producer or fitted quotient.
No lattice point is selected. Sources: IEEE754 and original Excel observations.
"""
import argparse
from collections import Counter
from pathlib import Path
import hashlib
import json

from exact import (BINARY64 as B, EXTENDED64 as E, EXTENDED53 as P,
                   Lattice, Span, binary64)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def certify(trace, gaps, output):
    sources=[Path(__file__),Path(__file__).with_name('exact.py')]
    frozen={str(p):digest(p)for p in [*sources,Path(trace),Path(gaps)]}
    ids=[json.loads(line)['id']for line in Path(gaps).read_text().splitlines()]
    if len(ids)!=24 or len(set(ids))!=24:
        raise ValueError('expected the original 24 distinct certified gap IDs')
    rows=[json.loads(line)for line in Path(trace).read_text().splitlines()]
    proofs=[];counts=Counter();seen=set()
    for row in rows:
        if row['id']not in ids:continue
        if row['id']in seen:raise ValueError('duplicate original witness')
        seen.add(row['id'])
        if row['function']!='PMT' or row['result_index']is not None:
            raise ValueError('changed original identity')
        r=abs(binary64(row['args'][0]));y=abs(binary64(row['expected']))
        if not r or not y>=B.value(B.half)or row['actual']==row['expected']:
            raise ValueError('normal nonzero original mismatch required')
        publication=B.pullback(Span(y,y))
        shells={}
        for precision,fmt in ((53,P),(64,E)):
            interval=fmt.pullback(publication).linear(1/r)
            entries={}
            for operand_precision in range(53,65):
                lattice=Lattice(operand_precision,-16382,16383)
                domain=lattice.coordinates(interval)
                count=domain[1]-domain[0]+1 if domain else 0
                entries[str(operand_precision)]=count
            # Prove a failed binary64 shell with its exact adjacent bracket.
            if entries['53']:
                raise ValueError('an alleged original binary64 quotient gap became reachable')
            low=B.floor_coordinate(interval.low)
            below,above=B.value(low),B.value(low+1)
            if not below<=interval.low<=interval.high<=above:
                raise ValueError('incorrect adjacent quotient gap bracket')
            if interval.contains(below)or interval.contains(above):
                raise ValueError('binary64 neighbor reaches target')
            if not entries['64']:
                raise ValueError('no extended quotient: inspect before claiming reachability')
            minimal=next(int(p)for p,c in entries.items()if c)
            counts[f'pc{precision}_min_operand_precision_{minimal}']+=1
            counts[f'pc{precision}_binary64_gaps']+=1
            counts[f'pc{precision}_extended_reachable']+=1
            shells[f'pc{precision}']=dict(quotient_magnitude_interval=interval.json(),
                lattice_point_counts=entries,minimal_significand_precision=minimal,
                binary64_adjacent_bracket=[[str(q.numerator),str(q.denominator)]for q in (below,above)])
        proofs.append(dict(witness=row,shells=shells))
    if seen!=set(ids):raise ValueError('missing original certified gap witness')
    if frozen!={p:digest(p)for p in frozen}:raise ValueError('sources/evidence changed')
    receipt=dict(complete=True,rows=len(proofs),counts=dict(counts),proofs=proofs,
        input_sha256=frozen,
        claim='Every original gap excludes every finite binary64 quotient under either final x87 precision setting. Every gap has at least one 80-bit quotient under each setting; only lattice counts/preimages are retained.',
        limits='Existential reachability is not a matched Excel function, chosen intermediate, executable producer or evidence of internal Excel precision. No fit, selector or product change.')
    Path(output).write_text(json.dumps(receipt,indent=2)+'\n')
    return receipt


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('trace','gaps','output'):p.add_argument(name)
    a=p.parse_args();r=certify(a.trace,a.gaps,a.output)
    print(json.dumps({k:v for k,v in r.items()if k!='proofs'}))
