from nl import *
import itertools
d=load(); t=d["t"]; x=d["x"]; p=d["p"]
AS=[[2.5090809287301226727e3,3.3430575583588128105e4,6.7265770927008700853e4,4.5921953931549871457e4,1.3731693765509461125e4,1.9715909503065514427e3,1.3314166789178437745e2,3.3871328727963666080],
 [5.2264952788528545610e3,2.8729085735721942674e4,3.9307895800092710610e4,2.1213794301586595867e4,5.3941960214247511077e3,6.8718700749205790830e2,4.2313330701600911252e1,1.0]]
def horner(c,r,mode):
    v=np.full_like(r,c[0])
    for a in c[1:]:
        if mode=="x87": v=add(mul(v,r),a)
        elif mode=="ext": v=np.asarray(v,dtype=LD)*np.asarray(r,dtype=LD)+LD(a)
        else: v=v*r+a
    return v
for name,m in [("B",(np.abs(t)>0.7072)&(np.abs(t)<=0.85)),("A",(np.abs(t)<0.7070)&(np.abs(t)>0.5)),("C",(t>0.8501))]:
    q=p[m]-0.5; xm=x[m]
    for rmode in ["d","x87"]:
        r=(0.180625-q*q) if rmode=="d" else add(0.180625,-mul(q,q))
        for mode in ["d","x87","ext"]:
            N=horner(AS[0],r,mode); D=horner(AS[1],r,mode)
            if mode=="ext":
                c1=st(LD(q)*N/D); c2=st(LD(q)*st(N)/st(D)); c3=st(LD(q)*LD(st(N/D)))
            else:
                c1=div(mul(q,N),D); c2=mul(q,div(N,D)); c3=div(q*N,D) if mode=="d" else c1
            print(name,m.sum(),"r",rmode,"h",mode," (qN)/D %.4f  q(N/D) %.4f  alt %.4f"%((c1==xm).mean(),(c2==xm).mean(),(c3==xm).mean()))
