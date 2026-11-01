"""All four-factor product/quotient trees and fixed timing identities, MIT."""
from itertools import permutations
import json


def expressions():
    def trees(order):
        if len(order)==1:
            powers=tuple(int(i==order[0])for i in range(4))
            yield order[0],powers
            return
        for split in range(1,len(order)):
            for a,ap in trees(order[:split]):
                for b,bp in trees(order[split:]):
                    yield ('mul',a,b),tuple(x+y for x,y in zip(ap,bp))
                    yield ('div',a,b),tuple(x-y for x,y in zip(ap,bp))
    found={json.dumps(tree,separators=(',',':'))for order in permutations(range(4))
           for tree,powers in trees(order)if powers==(1,1,-1,-1)}
    result=[json.loads(encoded)for encoded in sorted(found)]
    if len(result)!=60:raise ValueError('incomplete four-factor grammar')
    return result


def templates(Builder):
    result=[]
    for mode in ('growth','discount'):
        for numerator in range(2):
            for association in range(68):
                g=Builder(7);rate,pv,fv,timing,growth,em=0,2,3,4,5,6
                if numerator==0:
                    a,b=(pv,fv)if mode=='growth'else(fv,pv)
                    amount=g.op('add',g.op('mul',a,growth),b)
                else:
                    a=pv if mode=='growth'else fv
                    amount=g.op('add',g.op('add',pv,fv),g.op('mul',a,em))
                rt=g.op('mul',rate,timing,store=False)
                tf=g.op('add',g.one,rt)
                op=g.op
                if association<60:
                    tree=expressions()[association]
                    leaves=[amount,rate,em,tf]
                    def emit(t,root=False):
                        if isinstance(t,int):return leaves[t]
                        a=emit(t[1]);b=emit(t[2])
                        return op(t[0],a,b,store=not root)
                    final=emit(tree,True)
                    kind='four-factor'
                else:
                    tree=None;kind='timing';form=association-60
                    if form==0:
                        q=op('div',amount,em);d=op('div',rt,tf)
                        adjusted=op('sub',q,op('mul',q,d));final=op('mul',adjusted,rate,store=False)
                    elif form==1:
                        q=op('div',amount,em);d=op('div',rt,tf);qr=op('mul',q,rate)
                        final=op('sub',qr,op('mul',qr,d),store=False)
                    elif form==2:
                        q=op('div',amount,em);d=op('div',rt,tf)
                        factor=op('sub',rate,op('mul',rate,d));final=op('mul',q,factor,store=False)
                    elif form==3:
                        q=op('div',amount,em);product=op('mul',rate,rt)
                        factor=op('sub',rate,op('div',product,tf));final=op('mul',q,factor,store=False)
                    elif form==4:
                        q=op('div',amount,em);reciprocal=op('div',g.one,rate)
                        denominator=op('add',reciprocal,timing);final=op('div',q,denominator,store=False)
                    elif form in (5,6):
                        denominator=op('add',em,op('mul',em,rt))
                        if form==5:final=op('mul',op('div',amount,denominator),rate,store=False)
                        else:final=op('mul',amount,op('div',rate,denominator),store=False)
                    else:
                        q=op('div',amount,em);d=op('div',rt,tf)
                        complement=op('sub',g.one,d)
                        final=op('mul',op('mul',q,complement),rate,store=False)
                # Independently choose PC53/PC64 for the final arithmetic result,
                # while its full-width inputs retain the caller's computation.
                g.nodes[final]['phase']=1
                if mode=='growth':final=op('neg',final,store=False)
                result.append(dict(id=len(result),mode=mode,numerator=numerator,
                    association=association,kind=kind,expression=tree,graph=g.graph(final,final)))
    return result
