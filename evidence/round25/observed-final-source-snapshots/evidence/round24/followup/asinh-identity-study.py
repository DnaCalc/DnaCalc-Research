"""Fixed ASINH identity rounding, using existing independent arithmetic and LN."""
import pathlib,json,hashlib,struct,sys,subprocess,math
sys.path.insert(0,str(pathlib.Path.cwd()/"checker"))
from exact import binary64,operation,BINARY64 as B
from normal import sqrt_round,encoding
raw=pathlib.Path(sys.argv[1]);exe=pathlib.Path(sys.argv[2]);rows=[json.loads(l)for l in (raw/"asinh-bank.jsonl").open()]
models=["strict","stored_pc64","resident_pc64","stored_pc53","resident_pc53"]
questions=[];links=[];excluded={};guards=[]
for i,row in enumerate(rows):
 a=row["args"]
 if len(a)!=1 or not isinstance(a[0],str)or not a[0].startswith("0x"):
  excluded[i]="shape";continue
 try:x=binary64(a[0])
 except ValueError:excluded[i]="nonfinite";continue
 f=struct.unpack(">d",bytes.fromhex(a[0][2:]))[0]
 if 0<abs(f)<sys.float_info.min:excluded[i]="existing input refusal";continue
 if not math.isfinite(f*f):guards.append(i);assert row["actual"]=="error:Num";continue
 a=abs(x)
 for m,model in enumerate(models):
  square=operation(a,a,"mul",model);radicand=operation(square,1,"add",model);root=sqrt_round(radicand,model);inner=B.round(operation(a,root,"add",model));bits=encoding(inner)
  identity=f"asinh-{i}-{m}";questions.append(dict(id=identity,args=[bits]));links.append((i,m,identity))
request=raw/"asinh-ln-requests-final.json";request.write_text(json.dumps(dict(function="LN",witnesses=questions)))
out=raw/"asinh-ln-predictions-final.jsonl";subprocess.run([str(exe),"invoke-inputs",str(request),str(out)],check=True)
predictions=[json.loads(l)for l in out.open()];assert len(predictions)==len(links)
counts=[[0]*5 for _ in models];first_losses=[None]*len(models);changes=[]
for position,((i,m,identity),prediction) in enumerate(zip(links,predictions)):
 assert prediction["id"]==identity
 assert prediction["args"]==questions[position]["args"]
 row=rows[i];got=prediction["calc_func_prediction"]
 if got.startswith("0x")and int(row["args"][0],16)>>63 and int(got,16):got=f"0x{int(got,16)^(1<<63):016x}"
 if m==0:assert got==row["actual"],(row["id"],got,row["actual"])
 was=row["actual"]==row["expected"];now=got==row["expected"]
 for k,b in enumerate((True,now,was,now and not was,was and not now)):counts[m][k]+=int(b)
 if was and not now and first_losses[m]is None:first_losses[m]=dict(original=row,prediction=got)
 if got!=row["actual"]:changes.append(dict(row=i,id=row["id"],program=m,baseline=row["actual"],expected=row["expected"],prediction=got))
provenance=raw/"asinh-discovery-provenance.json"
provenance.write_bytes(subprocess.check_output([str(exe),"provenance"]))
result=dict(complete=True,input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [raw/"asinh-bank.jsonl",request,out,provenance]},transport_arguments_equal=True,models=models,counts_all_exact_baseline_gains_losses=counts,first_losses=first_losses,excluded=excluded,unchanged_square_overflow_rows=len(guards),zero_loss_improvements=[m for m,c in enumerate(counts)if c[3]and not c[4]],source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [pathlib.Path(__file__),pathlib.Path("checker/exact.py"),pathlib.Path("checker/normal.py"),exe]},limits="Five fixed precision/store programs for sign(x)*LN(abs(x)+sqrt(x*x+1)); existing square-overflow guard and subnormal admission remain unchanged. Independent rational arithmetic/squared-midpoint roots; worksheetLN is the original conditional numerical leaf. Strict identity must reproduce every admitted registered ASINH. No fitted coefficient, threshold, row selector or new Excel answer.")
(raw/"asinh-identity-summary-final.json").write_text(json.dumps(result,indent=2)+"\n");(raw/"asinh-identity-changes-final.jsonl").write_text("".join(json.dumps(r)+"\n"for r in changes));print(result)
