"""Local fixed GEOMEAN composition hypotheses; no new Excel inputs/answers."""
from pathlib import Path
from collections import Counter
import hashlib,json,struct,subprocess,sys
sys.path.insert(0,str(Path('/home/ubuntu/projects/DnaCalc/DnaCalc-Research/racer')))
from geomean import coordinates
source,executable,output=map(Path,sys.argv[1:]);executable=executable.resolve();output.mkdir()
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
frozen={str(p):sha(p)for p in(source,executable,Path(__file__))}
def bits(x):return '0x'+struct.pack('>d',x).hex()
def mul(xs,order):
 if order==2:
  if len(xs)==1:return xs[0]
  k=len(xs)//2;return mul(xs[:k],2)*mul(xs[k:],2)
 result=1.0
 for x in xs if order==0 else xs[::-1]:result*=x
 return result
requests=[];bindings=[];scope=Counter()
for index,line in enumerate(source.open()):
 raw=json.loads(line);inputs,reason=coordinates(raw)
 if reason:scope[reason]+=1;continue
 scope['scoped']+=1
 xs=[struct.unpack('>d',v.to_bytes(8,'big'))[0]for v in inputs]
 for order in range(3):
  qid=f'r18-local-power-{index}-{order}'
  args=[bits(mul(xs,order)),bits(1.0/len(xs))]
  requests.append(dict(id=qid,args=args))
  bindings.append(dict(id=qid,trace_index=index,order=order,source_id=raw['id'],expected=raw['expected'],
   frozen_actual=raw['actual'],args=args,not_excel_evidence=True))
request=output/'requests.json';request.write_text(json.dumps(dict(function='POWER',witnesses=requests,selection='Local composition hypotheses only; never submit as Excel evidence.'))+'\n')
subprocess.run([str(executable),'invoke-inputs',str(request),str(output/'predictions.jsonl')],check=True,timeout=300)
actuals=[json.loads(line)for line in(output/'predictions.jsonl').open()];assert len(actuals)==len(bindings)
counts=[Counter()for _ in range(3)];failures=[[]for _ in range(3)]
with(output/'application-route-inputs.jsonl').open('x')as out:
 for binding,actual in zip(bindings,actuals):
  assert actual['id']==binding['id']and actual['args']==binding['args']and actual['not_excel_evidence']is True
  got=actual['calc_func_prediction'];now=got==binding['expected'];old=binding['frozen_actual']==binding['expected'];c=counts[binding['order']]
  c['observations']+=1;c['matched']+=now;c['failed']+=not now;c['gains_vs_frozen']+=now and not old;c['losses_vs_frozen']+=old and not now
  if not now and len(failures[binding['order']])<4:failures[binding['order']].append(dict(binding,calculated_power=got))
  out.write(json.dumps(dict(id=binding['id'],function='POWER',args=binding['args'],result_index=None,numeric_cell_inputs=True,actual=got,not_excel_evidence=True))+'\n')
subprocess.run([str(executable),'function-routes',str(output/'application-route-inputs.jsonl'),str(output/'application-routes.jsonl')],check=True,timeout=300)
assert frozen=={p:sha(Path(p))for p in frozen}
result=dict(complete=True,profiles=3,scope=scope,counts=counts,counterexamples=failures,input_sha256=frozen,
 generated_request_sha256=sha(request),prediction_sha256=sha(output/'predictions.jsonl'),
 limits='Three fixed binary64 product orders (forward,reverse,balanced), then the existing registered POWER(product,RN53(1/n)) call. Typed inputs and actual POWER workbook routes are audited. This is a conditional shared-helper composition hypothesis, not an Excel internal call trace or a product gain. No row-dependent selector, fitted parameter or fresh validation; original observations are discovery. Frozen actual comparisons are scoped historical application observations, not the current full bank.')
with(output/'summary.json').open('x')as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(dict(scope=scope,counts=counts)))
