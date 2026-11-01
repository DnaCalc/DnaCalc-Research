"""Independent finite consumer grammar, derived by leaf-set partitions, MIT."""
from functools import lru_cache
import json


def expressions():
    @lru_cache(None)
    def partitioned(leaves):
        if len(leaves)==1:return (leaves[0],)
        result=[]
        for mask in range(1,(1<<len(leaves))-1):
            left=tuple(x for i,x in enumerate(leaves)if mask>>i&1)
            right=tuple(x for i,x in enumerate(leaves)if not(mask>>i&1))
            for a in partitioned(left):
                for b in partitioned(right):
                    result.extend((('mul',a,b),('div',a,b)))
        return tuple(result)
    def degree(tree,sign=1):
        if isinstance(tree,int):return {tree:sign}
        a=degree(tree[1],sign);b=degree(tree[2],-sign if tree[0]=='div'else sign)
        return {i:a.get(i,0)+b.get(i,0)for i in range(4)}
    result=sorted({json.dumps(t,separators=(',',':'))for t in partitioned((0,1,2,3))
                   if degree(t)=={0:1,1:1,2:-1,3:-1}})
    if len(result)!=60:raise ValueError('independent expression count differs')
    return [json.loads(t)for t in result]


def templates(context,append,finish):
    result=[];exprs=expressions()
    for mode in ('growth','discount'):
        for numerator in range(2):
            for association in range(68):
                ns=context(7);one=7
                def op(kind,a=0,b=0,stored=True):return append(ns,kind,a,b,stored)
                balance,other=(2,3)if mode=='growth'else(3,2)
                if numerator==0:
                    product=op('mul',balance,5);amount=op('add',product,other)
                else:
                    total=op('add',2,3);product=op('mul',balance,6);amount=op('add',total,product)
                rt=op('mul',0,4,stored=False);tf=op('add',one,rt)
                if association<60:
                    tree=exprs[association]
                    terminals=(amount,0,6,tf)
                    def interpret(t,outer=False):
                        if isinstance(t,int):return terminals[t]
                        left=interpret(t[1]);right=interpret(t[2]);return op(t[0],left,right,stored=not outer)
                    final=interpret(tree,True);kind='four-factor'
                else:
                    tree=None;kind='timing';case=association-60
                    if case in (0,1,2):
                        q=op('div',amount,6);ratio=op('div',rt,tf)
                        if case==0:
                            term=op('mul',q,ratio);adjusted=op('sub',q,term);final=op('mul',adjusted,0,stored=False)
                        elif case==1:
                            product=op('mul',q,0);term=op('mul',product,ratio);final=op('sub',product,term,stored=False)
                        else:
                            term=op('mul',0,ratio);factor=op('sub',0,term);final=op('mul',q,factor,stored=False)
                    elif case==3:
                        q=op('div',amount,6);product=op('mul',0,rt);term=op('div',product,tf)
                        factor=op('sub',0,term);final=op('mul',q,factor,stored=False)
                    elif case==4:
                        q=op('div',amount,6);reciprocal=op('div',one,0);denominator=op('add',reciprocal,4)
                        final=op('div',q,denominator,stored=False)
                    elif case in (5,6):
                        term=op('mul',6,rt);denominator=op('add',6,term)
                        if case==5:
                            quotient=op('div',amount,denominator);final=op('mul',quotient,0,stored=False)
                        else:
                            rate=op('div',0,denominator);final=op('mul',amount,rate,stored=False)
                    else:
                        q=op('div',amount,6);ratio=op('div',rt,tf);factor=op('sub',one,ratio)
                        adjusted=op('mul',q,factor);final=op('mul',adjusted,0,stored=False)
                ns[final]['phase']=1
                if mode=='growth':final=op('neg',final,stored=False)
                result.append(dict(id=len(result),mode=mode,numerator=numerator,association=association,
                    kind=kind,expression=tree,graph=finish(ns,final,final)))
    return result
