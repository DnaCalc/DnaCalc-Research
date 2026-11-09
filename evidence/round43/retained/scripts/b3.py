import sys
sys.path.insert(0,"/home/ubuntu/codex-runs/r42/pylib")
import mpmath as mp
from nl import *
mp.mp.dps=50
d=load(); t=d["t"]; x=d["x"]; p=d["p"]
Pc=["2.5090809287301226727e3","3.3430575583588128105e4","6.7265770927008700853e4","4.5921953931549871457e4","1.3731693765509461125e4","1.9715909503065514427e3","1.3314166789178437745e2","3.3871328727963666080"]
Qc=["5.2264952788528545610e3","2.8729085735721942674e4","3.9307895800092710610e4","2.1213794301586595867e4","5.3941960214247511077e3","6.8718700749205790830e2","4.2313330701600911252e1","1.0"]
Cc=["7.74545014278341407640e-4","2.27238449892691845833e-2","2.41780725177450611770e-1","1.27045825245236838258","3.64784832476320460504","5.76949722146069140550","4.63033784615654529590","1.42343711074968357734"]
Dc=["1.05075007164441684324e-9","5.47593808499534494600e-4","1.51986665636164571966e-2","1.48103976427480074590e-1","6.89767334985100004550e-1","1.67638483018380384940","2.05319162663775882187","1.0"]
s2=mp.sqrt(2)
A=[float(mp.mpf(c)/(2*s2)) for c in Pc]; B=[float(mp.mpf(c)) for c in Qc]
C=[float(mp.mpf(c)) for c in Cc]; Dd=[float(mp.mpf(c)*s2) for c in Dc]
def horner(c,r,mode):
    v=np.full_like(r,c[0])
    for a in c[1:]: v=add(mul(v,r),a) if mode=="x87" else v*r+a
    return v
for name,m in [("B",(np.abs(t)>0.7072)&(np.abs(t)<=0.85)),("A",(np.abs(t)<0.7070)&(np.abs(t)>0.3))]:
    tt=t[m]; xm=x[m]
    for mode in ["d","x87"]:
        r=0.180625-0.25*tt*tt if mode=="d" else add(0.180625,-mul(0.25,mul(tt,tt)))
        z1=horner(A,r,mode); z2=horner(B,r,mode)
        for form in ["(t*z1)/z2","t*(z1/z2)"]:
            y=(tt*z1)/z2 if form[0]=="(" else tt*(z1/z2)
            if mode=="x87": y=div(mul(tt,z1),z2) if form[0]=="(" else mul(tt,div(z1,z2))
            print(name,mode,form,"x=y*s2 %.4f"%(mul(y,S2)==xm).mean(), " x=d*s2 %.4f"%((y*S2)==xm).mean())
m=t>0.8501; tt=t[m]; xm=x[m]
for mode in ["d","x87"]:
    ln2=float(mp.log(2))
    w=np.sqrt(ln2-np.log(1.0-tt)) if mode=="d" else np.sqrt(add(ln2,-np.log(1.0-tt)))
    r=w-1.6
    z1=horner(C,r,mode); z2=horner(Dd,r,mode)
    y=z1/z2 if mode=="d" else div(z1,z2)
    print("C",mode,"x=y*s2 %.4f"%(mul(y,S2)==xm).mean(), "plain %.4f"%((y*S2)==xm).mean())
    # also: AS241 tail at 1-p directly
    w2=np.sqrt(-np.log(1.0-p[m])); r2=w2-1.6
    cc=[float(mp.mpf(c)) for c in Cc]; dd=[float(mp.mpf(c)) for c in Dc]
    xx=div(horner(cc,r2,mode),horner(dd,r2,mode)) if mode=="x87" else horner(cc,r2,mode)/horner(dd,r2,mode)
    print("C AS241(1-p)",mode,"%.4f"%(xx==xm).mean())
