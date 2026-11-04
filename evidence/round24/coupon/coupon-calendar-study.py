import json,struct,pathlib,datetime,calendar,collections,itertools,hashlib
root=pathlib.Path(__file__).parent
rows=[json.loads(l)for l in (root/'coupdays-bank.jsonl').read_text().splitlines()]
def dec(x):return struct.unpack('>d',bytes.fromhex(x[2:]))[0] if x.startswith('0x') else None
def dt(x):return datetime.datetime(1899,12,30)+datetime.timedelta(days=int(x))
def eom(d):return d.day==calendar.monthrange(d.year,d.month)[1]
def shift(d,mo,end=False):
 k=d.year*12+d.month-1+mo;y,m=divmod(k,12);m+=1
 return d.replace(year=y,month=m,day=calendar.monthrange(y,m)[1] if end else min(d.day,calendar.monthrange(y,m)[1]))
programs=list(itertools.product(range(6),range(7)))
counts=[[0]*5 for _ in programs];kills=[None]*len(programs);excluded=collections.Counter()
for r in rows:
 a=[dec(v)for v in r['args']]
 if len(a)<4 or any(v is None for v in a)or a[3]!=1 or a[2]not in(1,2,4)or a[0]>=a[1]or a[0]<61:
  excluded['outside_actual_actual_finite_calendar_scope']+=1;continue
 try: settle,maturity=map(dt,a[:2])
 except (ValueError,OverflowError):excluded['outside_python_calendar']+=1;continue
 months=12//int(a[2]);want=dec(r['expected']);old=dec(r['actual'])
 for i,(walk,measure) in enumerate(programs):
  after=maturity;k=1
  while True:
   reference=maturity if walk<3 else after
   end=eom(maturity) if walk%3==0 else False if walk%3==1 else eom(reference)
   before=shift(reference,-k*months if walk<3 else -months,end)
   if before<=settle:break
   k+=1;after=before
  right,left=(after,before) if measure==0 else (shift(before,months),before) if measure==1 else (after,shift(after,-months)) if measure==2 else (shift(before,months,eom(before)),before) if measure==3 else (after,shift(after,-months,eom(after))) if measure==4 else (shift(before,months,eom(maturity)),before) if measure==5 else (after,shift(after,-months,eom(maturity)))
  got=float((right-left).days);matched=got==want;was=old==want
  for j,v in enumerate([1,matched,was,matched and not was,was and not matched]):counts[i][j]+=int(v)
  if not matched and kills[i]is None:kills[i]={'id':r['id'],'args':a,'got':got,'expected':want,'baseline':old}
result={'complete':True,'programs':programs,'counts_all_exact_baseline_gains_losses':counts,'first_contradictions':kills,'excluded':dict(excluded),'zero_loss_improvements':[i for i,c in enumerate(counts)if c[3]>0 and c[4]==0],'limits':'42 coherent calendar shift/walk programs. Research-only dates>=61, basis1 numeric inputs. No answer selector or fitted threshold. Original date/month shift definitions; no product edit.'}
(root/'coupon-calendar-summary.json').write_text(json.dumps(result,indent=2)+'\n')
print([(i,programs[i],c)for i,c in enumerate(counts)if c[3]>0 and c[4]==0])
