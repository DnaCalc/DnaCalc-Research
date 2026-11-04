"""R24 counterfactuals through the existing independent variance checker."""
import collections, hashlib, json, pathlib, sys
sys.path.insert(0,str(pathlib.Path.cwd()/"checker"))
from variance import evaluate
from exact import binary64
from normal import encoding
raw=pathlib.Path(sys.argv[1]);scoped=raw/"variance-current-scoped/scoped.jsonl"
counts=collections.defaultdict(lambda:[0,0,0,0,0]);first_losses={};controls=0
for line in scoped.open():
 row=json.loads(line);old=row["raw"];values=list(map(binary64,row["values_bits"]))
 for algorithm in (0,7,8):
  profile=dict(algorithm=algorithm,linear_order=0,square_order=0,models=["strict"]*5)
  root=evaluate(profile,values,old["function"]);got="error:Num" if root is None else encoding(root)
  if algorithm==0:
   assert got==old["actual"],(old["id"],got,old["actual"])
   controls+=1
  was=old["actual"]==old["expected"];now=got==old["expected"];key=old["function"]+":"+str(algorithm)
  for i,b in enumerate((True,now,was,now and not was,was and not now)):counts[key][i]+=int(b)
  if was and not now and key not in first_losses:first_losses[key]=dict(id=old["id"],formula=old["formula"],expected=old["expected"],baseline=old["actual"],prediction=got)
result=dict(complete=True,registered_controls=controls,programs={"0":"original centered two-pass","7":"(sum(x*x)-sum(x)^2/n)/d","8":"(n*sum(x*x)-sum(x)^2)/(n*d)"},counts_all_exact_baseline_gains_losses=dict(counts),first_losses=first_losses,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()for p in [pathlib.Path(__file__),pathlib.Path("checker/variance.py"),pathlib.Path("checker/exact.py"),pathlib.Path("checker/normal.py")]},input_sha256=hashlib.sha256(scoped.read_bytes()).hexdigest(),limits="Existing independent exact rational variance grammar, three fixed strict formulas only; no fitted selector or new numerical program. Original numeric scope/exclusions from variance-prepare. Short original literals refute a15-digit-limit explanation but do not prove Excel operand bits; new operand echoes are still questions. Every changing uncentered formula loses existing matches; no production adoption.")
(raw/"variance-current-counterfactuals.json").write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result))
