"""Seasonal state identification after the nonseasonal reconstruction.

Original tier-1 question research. Future phase selections expose the final
trend and zero-mean seasonal states; inverse transitions test initialisation
without choosing a polynomial degree or adding fitted pieces.
"""
import argparse,itertools,random
from .ets import *

def design():
 rng=random.Random(397801);cases=[]
 for period in [2,3,4,5,6,8]:
  for cycles in [2,3,4]:
   n=period*cycles
   if not 8<=n<=32:continue
   pattern=[rng.randrange(-16,17)/4 for _ in range(period)]
   for kind in ['periodic','trend','impulse','noise']:
    y=[pattern[j%period]+(j/8 if kind=='trend'else 0.)for j in range(n)]
    if kind=='impulse':y[n//2]+=1/16
    if kind=='noise':y=[v+rng.randrange(-8,9)/32 for v in y]
    cases.append(dict(y=y,season=period,kind=kind,cycles=cycles))
 write(OUT/'step7-seasonal-design.json',dict(cases=cases,producer_source=Path(__file__).read_text()))
 c,s=service();write(OUT/'step7-seasonal-stat.json',ask_cases(c,s,cases,'FORECAST.ETS.STAT',range(1,9)))
 # A common horizon list includes two observations of every phase.
 write(OUT/'step7-seasonal-forecast.json',ask_cases(c,s,cases,'FORECAST.ETS',range(10)))
 write(OUT/'step7-seasonal-confidence.json',ask_cases(c,s,cases,'FORECAST.ETS.CONFINT',[(1.,.95)]))
 print('seasonal cases',len(cases),flush=True)

def initial(y,p,trend_kind,season_kind):
 n=len(y);cycles=n//p
 if trend_kind=='raw':t=initial_trend(y)
 elif trend_kind=='zero':t=0.
 elif trend_kind=='firsttwo':t=sequential_total(y[j+p]-y[j]for j in range(p))/(p*p)
 elif trend_kind=='end':t=(sequential_total(y[-p:])-sequential_total(y[:p]))/(p*(n-p))
 elif trend_kind=='pairs':t=sequential_total(y[j+p]-y[j]for j in range(n-p))/(p*(n-p))
 elif trend_kind=='cycles':
  means=[sequential_total(y[j*p:(j+1)*p])/p for j in range(cycles)];t=initial_trend(means)/p
 elif trend_kind=='within':
  num=0.;den=0.
  for q in range(p):
   xs=list(range(q,n,p));mx=sequential_total(xs)/len(xs);my=sequential_total(y[j]for j in xs)/len(xs)
   for j in xs:num+=(j-mx)*(y[j]-my);den+=(j-mx)*(j-mx)
  t=num/den
 else:raise ValueError(trend_kind)
 if season_kind=='first':se=list(y[:p])
 elif season_kind=='first_detrend':se=[y[j]-j*t for j in range(p)]
 elif season_kind=='all_detrend':se=[sequential_total(y[j]-j*t for j in range(q,n,p))/len(range(q,n,p))for q in range(p)]
 elif season_kind=='all':se=[sequential_total(y[j]for j in range(q,n,p))/len(range(q,n,p))for q in range(p)]
 else:raise ValueError(season_kind)
 center=sequential_total(se)/p;se=[v-center for v in se];l=y[0]-se[0]
 return l,t,se

def fit(y,p,a,b,g,init):
 l,t,se=init;se=se.copy();errors=[0.];predictions=[y[0]]
 for j,v in enumerate(y[1:],1):
  q=j%p;prediction=(l+t)+se[q];e=v-prediction;l=l+(t+a*e);t=t+b*e;se[q]=se[q]+g*e;errors.append(e);predictions.append(prediction)
 center=sequential_total(se)/p;l+=center;se=[v-center for v in se]
 return l,t,se,errors,predictions

def analyze():
 d=read(OUT/'step7-seasonal-stat.json');st=grouped(d);fo=grouped(read(OUT/'step7-seasonal-forecast.json'));programs=[];inverse=[]
 for i,z in enumerate(d['cases']):
  y=z['y'];p=z['season'];a,b,g=[value(st[i][str(k)]['bits'])for k in [1,2,3]];f=[value(fo[i][str(k)]['bits'])for k in range(10)]
  t=(f[p]-f[0])/p;l=sequential_total(f[h]-h*t for h in range(p))/p;se=[0.]*p
  for h in range(p):se[(len(y)-1+h)%p]=f[h]-h*t-l
  trace=[l,t,se.copy()]
  if a+g!=1:
   for j in range(len(y)-1,0,-1):
    q=j%p;e=(y[j]-(l+se[q]))/(1-a-g);t=t-b*e;se[q]=se[q]-g*e;l=l-t-a*e
   inverse.append(dict(case=i,parameters=[a,b,g],initial_level=l,initial_trend=t,initial_seasons=se,initial_prediction=l+se[0],observed_y0=y[0],final=trace))
 for tk,sk in itertools.product(['zero','raw','firsttwo','end','pairs','cycles','within'],['first','first_detrend','all','all_detrend']):
  exact=[];worst=0.;first=[]
  for i,z in enumerate(d['cases']):
   a,b,g=[value(st[i][str(k)]['bits'])for k in [1,2,3]];l,t,se,e,pr=fit(z['y'],z['season'],a,b,g,initial(z['y'],z['season'],tk,sk))
   predicted=[(l+h*t)+se[(len(z['y'])-1+h)%z['season']]for h in range(10)];actual=[value(fo[i][str(h)]['bits'])for h in range(10)]
   err=max(abs(x-y)/max(1.,abs(y))for x,y in zip(predicted,actual));worst=max(worst,err)
   if all(bits(x)==bits(y)for x,y in zip(predicted,actual)):exact.append(i)
   elif len(first)<3:first.append(dict(case=i,kind=z['kind'],relative_error=err,predicted=[bits(v)for v in predicted],actual=[bits(v)for v in actual]))
  programs.append(dict(trend=tk,seasons=sk,all_exact=len(exact)==len(d['cases']),exact_cases=exact,worst_relative_error=worst,first_refutations=first))
 write(OUT/'step7-initial-programs.json',dict(programs=programs,inverse_states=inverse,producer_source=Path(__file__).read_text(),limits='Inverse states are approximate cancellation-sensitive diagnostics; they are never counted as bit-exact answers. STAT parameters condition these candidate recurrences; postfit parameter resets can invalidate the inversion.'))
 print('all-exact programs',[(p['trend'],p['seasons'])for p in programs if p['all_exact']]);print('smallest full-set residual',min(p['worst_relative_error']for p in programs))

def basis():
 cases=[]
 for n,p in [(8,2),(12,2),(12,3),(12,4),(16,4)]:
  pattern=[float(2*j-p+1)for j in range(p)];base=[pattern[j%p]for j in range(n)]
  for at in [-1,*range(n)]:
   y=base.copy()
   if at>=0:y[at]+=1/16
   cases.append(dict(y=y,season=p,kind='basis',at=at))
 write(OUT/'step7-basis-design.json',dict(cases=cases,producer_source=Path(__file__).read_text()))
 c,s=service();write(OUT/'step7-basis-stat.json',ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]))
 write(OUT/'step7-basis-forecast.json',ask_cases(c,s,cases,'FORECAST.ETS',range(10)))
 print('basis cases',len(cases),flush=True)

def basis_analysis(prefix='step7-basis'):
 from fractions import Fraction
 d=read(OUT/(prefix+'-stat.json'));st=grouped(d);fo=grouped(read(OUT/(prefix+'-forecast.json')));rows=[];baselines={}
 for i,z in enumerate(d['cases']):
  y=z['y'];p=z['season'];n=len(y);a,b,g=[value(st[i][str(k)]['bits'])for k in [1,2,3]];f=[value(fo[i][str(k)]['bits'])for k in range(10)]
  t=(f[p]-f[0])/p;l=sequential_total(f[h]-h*t for h in range(p))/p;se=[0.]*p
  for h in range(p):se[(n-1+h)%p]=f[h]-h*t-l
  denom=1-a-g
  for j in range(n-1,0,-1):
   q=j%p;e=(y[j]-(l+se[q]))/denom;t-=b*e;se[q]-=g*e;l=l-t-a*e
  center=sequential_total(se)/p;l+=center;se=[v-center for v in se]
  key=(n,p)
  if z['at']==-1:baselines[key]=(l,t,se)
  base=baselines[key];ts=(t-base[1])*16;ss=[(v-u)*16 for v,u in zip(se,base[2])]
  row=dict(case=i,n=n,period=p,at=z['at'],parameters=[a,b,g],initial=[l,t,se],initial_prediction_error=l+se[0]-y[0],inverse_amplification=abs(denom)**-(n-1),trend_coefficient=str(Fraction(ts).limit_denominator(65536)),season_coefficients=[str(Fraction(v).limit_denominator(65536))for v in ss],raw_trend_coefficient=ts,raw_season_coefficients=ss)
  rows.append(row)
 write(OUT/(prefix+'-kernels.json'),dict(rows=rows,producer_source=Path(__file__).read_text(),limits='Diagnostic inverse kernels, not exact proofs. Conditioning and rounding can invalidate recovered rational coefficients; every proposed shape requires fresh forward predictions. The positive-trend design separates the postfit reset that can invalidate the earlier zero-trend inversions.'))
 for n,p in baselines:
  print('n/period',n,p)
  for z in rows:
   if(z['n'],z['period'])==(n,p):print(z['at'],z['parameters'],z['trend_coefficient'],z['season_coefficients'],'init-error',z['initial_prediction_error'],'amp',z['inverse_amplification'])

def zero_residual_proof():
 from .ets_proof import freeze,finish
 rng=random.Random(397932);cases=[];forecasts=[];stats=[];confidence=[]
 for i in range(128):
  p=rng.choice([2,4,8]);cycles=rng.randrange(max(2,(8+p-1)//p),32//p+1);n=p*cycles;ex=rng.randrange(-300,301)
  pattern=[rng.randrange(-128,129)/8 for _ in range(p)]
  if 96<=i<112:pattern=[pattern[0]]*p
  if i>=112:pattern=pattern[:p//2]*2
  drift=rng.randrange(-16,17)/16 if i<80 else 0.;offset=rng.randrange(-512,513)/4
  y=[math.ldexp(offset+pattern[j%p]+drift*j,ex)for j in range(n)]
  init=initial(y,p,'firsttwo','first_detrend')
  if init[0]+init[2][0]!=y[0]:raise ValueError('initial observation not exact')
  l,t,se,e,pr=fit(y,p,.751,.001,.001,init)
  if any(v!=0. for v in e):raise ValueError('constructed innovations not zero')
  cases.append(dict(y=y,season=p,cycles=cycles,drift=math.ldexp(drift,ex),kind='exact_seasonal_affine'))
  forecasts.extend(bits((l+h*t)+se[(n-1+int(h))%p])for h in [0.,1.,2.,5.])
  stats.extend(bits(v)for v in [.751,.001,.001,0.,0.,0.,0.,1.])
  confidence.extend(bits(0.)for _ in range(4))
 name='step7-zero-residual'
 for kind,pred in [('forecast',forecasts),('stat',stats),('confidence',confidence)]:
  freeze(name+'-'+kind,dict(cases=cases,predictions=pred,extra_source=Path(__file__).read_text(),method='Fixed complete-cycle dyadic region: first-two-season trend, detrended first-season states, one innovations recurrence. No measured parameters are inputs. Every constructed innovation is zero before asking Excel; includes trends, constants and shorter-period aliases.'))
 finish(name+'-forecast',cases,forecasts,'FORECAST.ETS',[0.,1.,2.,5.])
 finish(name+'-stat',cases,stats,'FORECAST.ETS.STAT',range(1,9))
 finish(name+'-confidence',cases,confidence,'FORECAST.ETS.CONFINT',[(h,.95)for h in [1.,2.,5.,10.]])

def detection_map():
 rng=random.Random(398117);cases=[]
 for n,trend in itertools.product([8,12,16,24,32],[0.,1/8,1.,16.,128.,65536.]):
  cases.append(dict(y=[(-1. if j%2 else 1.)+trend*j for j in range(n)],kind='alternating',design=[n,trend]))
 for n,p in itertools.product([16,24,32],[3.,4.,4.5,5.,6.,6.5,8.]):
  cases.append(dict(y=[math.sin(2*math.pi*j/p)for j in range(n)],kind='sinusoid',design=[n,p]))
 noise=[rng.randrange(-128,129)/16 for _ in range(32)]
 for ex in range(-12,13):
  amplitude=2.**ex;cases.append(dict(y=[v+(-amplitude if j%2 else amplitude)for j,v in enumerate(noise)],kind='noise_amplitude',design=[ex]))
 for p in range(3,11):
  pattern=[rng.randrange(-16,17)/4 for _ in range(p)]
  for n in [2*p-1,2*p,3*p]:cases.append(dict(y=[pattern[j%p]for j in range(n)],kind='length_boundary',design=[n,p]))
 write(OUT/'step8-detection-design.json',dict(cases=cases,producer_source=Path(__file__).read_text()))
 c,s=service();d=ask_cases(c,s,cases,'FORECAST.ETS.SEASONALITY',[0]);write(OUT/'step8-detection-answers.json',d)
 for z,a in zip(cases,d['answers']):print(z['kind'],z['design'],value(a['bits'])if a['type']=='number'else a)

def zero_edges():
 from .ets_proof import freeze,finish
 cases=[dict(y=[v]*n,season=p,kind='zero_or_exponent_edge')for n,p in [(8,2),(8,4),(16,8),(24,4),(24,8),(32,2),(32,8)]for v in [0.,2.**-400,2.**400]]
 name='step7-zero-edges'
 specs=[('stat',[bits(v)for _ in cases for v in [.751,.001,.001,0.,0.,0.,0.,1.]],'FORECAST.ETS.STAT',list(range(1,9))),
        ('forecast',[bits(z['y'][0])for z in cases for _ in range(4)],'FORECAST.ETS',[0.,1.,2.,5.]),
        ('confidence',[bits(0.)for _ in cases for _ in range(3)],'FORECAST.ETS.CONFINT',[(1.,.95),(1.,math.nextafter(1.,0.)),(10.,.95)])]
 for kind,pred,function,variants in specs:freeze(name+'-'+kind,dict(cases=cases,predictions=pred,extra_source=Path(__file__).read_text(),method='Explicit complete-cycle all-zero and constant histories at both retained magnitude boundaries; includes the rounded upper confidence endpoint.'))
 for kind,pred,function,variants in specs:finish(name+'-'+kind,cases,pred,function,variants)

def automatic_two_proof():
 from .ets_proof import freeze,finish
 rng=random.Random(398263);cases=[];forecasts=[];stats=[];confidence=[]
 for i in range(96):
  n=2*rng.randrange(4,17);ex=rng.randrange(-300,301);amplitude=rng.randrange(1,65)/8;offset=rng.randrange(-4096,4097)/4
  trend=rng.choice([0.,1/8,1.,16.,128.,65536.])*rng.choice([-1.,1.])
  y=[math.ldexp(offset+(-amplitude if j%2 else amplitude)+trend*j,ex)for j in range(n)]
  l,t,se,e,pr=fit(y,2,.751,.001,.001,initial(y,2,'firsttwo','first_detrend'))
  if any(v!=0. for v in e):raise ValueError('constructed automatic fit is not exact')
  cases.append(dict(y=y,season=1,kind='exact_alternating_with_trend'))
  forecasts.extend(bits((l+h*t)+se[(n-1+int(h))%2])for h in [0.,1.,2.,5.])
  stats.extend(bits(v)for v in [.751,.001,.001,0.,0.,0.,0.,1.]);confidence.extend(bits(0.)for _ in range(4))
 name='step8-automatic-two'
 specs=[('forecast',forecasts,'FORECAST.ETS',[0.,1.,2.,5.]),('stat',stats,'FORECAST.ETS.STAT',list(range(1,9))),('confidence',confidence,'FORECAST.ETS.CONFINT',[(h,.95)for h in [1.,2.,5.,10.]]),('period',[bits(2.)]*len(cases),'FORECAST.ETS.SEASONALITY',[0])]
 for kind,pred,function,variants in specs:freeze(name+'-'+kind,dict(cases=cases,predictions=pred,extra_source=Path(__file__).read_text(),method='Unconditional exact alternating affine histories: period two and the same zero-innovation fit under automatic seasonality; fresh offsets, signs, exponents and trend/signal ratios.'))
 for kind,pred,function,variants in specs:finish(name+'-'+kind,cases,pred,function,variants)

def detection_analysis():
 d=read(OUT/'step8-detection-answers.json');actual=[value(a['bits'])if a['type']=='number'else None for a in d['answers']];programs=[]
 for preprocess,norm,search in itertools.product(['mean','linear','difference'],['biased','unbiased','cosine','pearson'],['first','global']):
  predictions=[];bad=[]
  for i,z in enumerate(d['cases']):
   y=z['y'];n=len(y);mean=sequential_total(y)/n;t=initial_trend(y)
   x=[v-mean for v in y]if preprocess=='mean'else[v-(mean+(j-(n-1)/2)*t)for j,v in enumerate(y)]if preprocess=='linear'else[v-u for u,v in zip(y,y[1:])]
   if preprocess=='difference':
    mx=sequential_total(x)/len(x);x=[v-mx for v in x]
   c=[1.];energy=sequential_total(v*v for v in x)
   for lag in range(1,len(x)):
    left=x[:-lag];right=x[lag:];num=sequential_total(a*b for a,b in zip(left,right))
    if norm=='biased':den=energy
    elif norm=='unbiased':den=energy*(len(x)-lag)/len(x)
    elif norm=='cosine':den=math.sqrt(sequential_total(v*v for v in left)*sequential_total(v*v for v in right))
    else:
     ml=sequential_total(left)/len(left);mr=sequential_total(right)/len(right)
     num=sequential_total((a-ml)*(b-mr)for a,b in zip(left,right));den=math.sqrt(sequential_total((v-ml)**2 for v in left)*sequential_total((v-mr)**2 for v in right))
    c.append(num/den if den else 0.)
   peaks=[k for k in range(2,min(n//3,len(c)-2)+1)if c[k]>0 and c[k]>c[k-1]and c[k]>=c[k+1]]
   pred=(peaks[0]if search=='first'else max(peaks,key=lambda k:c[k]))if peaks else 0
   predictions.append(pred)
   if pred!=actual[i]:bad.append(dict(case=i,kind=z['kind'],design=z['design'],predicted=pred,actual=actual[i],correlations=c))
  programs.append(dict(preprocessing=preprocess,normalization=norm,search=search,predictions=predictions,refutations=bad,all_exact=not bad))
 write(OUT/'step8-correlation-programs.json',dict(programs=programs,producer_source=Path(__file__).read_text(),limits='Fixed simple correlation candidates, no fitted thresholds. Every refutation is retained; a candidate must survive fresh questions before implementation.'))
 for z in programs:print(z['preprocessing'],z['normalization'],z['search'],'refutations',len(z['refutations']),[(r['case'],r['predicted'],r['actual'])for r in z['refutations'][:4]])

def automatic_lattice_edges():
 from .ets_proof import freeze,finish
 cases=[];forecasts=[];stats=[];confidence=[]
 for n in [8,12,16,24,32]:
  bound=2.**50/n**3;dc=2.**math.floor(math.log2(bound/4));slope=2.**math.floor(math.log2(bound/(4*n)))
  for offset,trend in [(dc,0.),(-dc,0.),(0.,slope),(dc/2,-slope)]:
   y=[offset+(1. if j%2==0 else -1.)+trend*j for j in range(n)]
   l,t,se,e,pr=fit(y,2,.751,.001,.001,initial(y,2,'firsttwo','first_detrend'))
   if any(v!=0. for v in e):raise ValueError('edge construction not exact')
   cases.append(dict(y=y,season=1,kind='automatic_lattice_edge',offset=offset,trend=trend))
   forecasts.extend(bits((l+h*t)+se[(n-1+int(h))%2])for h in [0.,1.,2.,5.])
   stats.extend(bits(v)for v in [.751,.001,.001,0.,0.,0.,0.,1.]);confidence.extend(bits(0.)for _ in range(4))
 name='step8-automatic-lattice-edges'
 specs=[('period',[bits(2.)]*len(cases),'FORECAST.ETS.SEASONALITY',[0]),('forecast',forecasts,'FORECAST.ETS',[0.,1.,2.,5.]),('stat',stats,'FORECAST.ETS.STAT',list(range(1,9))),('confidence',confidence,'FORECAST.ETS.CONFINT',[(h,.95)for h in [1.,2.,5.,10.]])]
 for kind,pred,function,variants in specs:freeze(name+'-'+kind,dict(cases=cases,predictions=pred,extra_source=Path(__file__).read_text(),method='DC and trend values near the existing n-cubed integer-lattice limit, amplitude one, before any detector-specific bound. All histories satisfy the explicit fit region.'))
 for kind,pred,function,variants in specs:finish(name+'-'+kind,cases,pred,function,variants)
 neighbors=[]
 for z in cases[::4]:
  for direction in [-math.inf,math.inf]:
   y=[math.nextafter(v,direction)if j%2==0 else v for j,v in enumerate(z['y'])]
   neighbors.append(dict(y=y,season=1,kind='single_bit_neighbor_outside_exact_lattice'))
 c,s=service();write(OUT/'step8-lattice-neighbor-answers.json',ask_cases(c,s,neighbors,'FORECAST.ETS.SEASONALITY',[0]))

def nonzero_controls(mixed=False):
 from .ets_proof import freeze,finish
 from .ets_confidence import variance
 rng=random.Random(397159 if mixed else 397147);cases=[]
 for i in range(64):
  ex=rng.randrange(-250,251);offset=rng.randrange(-512,513)/8;amplitude=rng.randrange(1,65)/8;trend=rng.choice([64.,128.,256.]);size=rng.choice([-1.,1.])*rng.randrange(1,17)/16
  second=rng.randrange(-8,9)/16 if mixed else 0.
  y=[offset+(-amplitude if j%2==0 else amplitude)+trend*j for j in range(12)]
  for j,change in [(7,size+7*second),(9,-2*size-3*second),(11,size)]:y[j]+=change
  cases.append(dict(y=[math.ldexp(v,ex)for v in y],season=2,initial_level=math.ldexp(offset,ex),initial_trend=math.ldexp(trend,ex),initial_seasons=[math.ldexp(-amplitude,ex),math.ldexp(amplitude,ex)]))
 name='step7-nonzero-mixed'if mixed else'step7-nonzero-controls'
 c,s=service();parameters=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);write(OUT/('step7-nonzero-mixed-parameters.json'if mixed else'step7-nonzero-parameters.json'),parameters);st=grouped(parameters)
 forecasts=[];stats=[];confidence=[];z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits'])
 for cs,ss in zip(cases,st):
  a,b,g=[value(ss[str(k)]['bits'])for k in [1,2,3]];l,t,se,e,pr=fit(cs['y'],2,a,b,g,(cs['initial_level'],cs['initial_trend'],cs['initial_seasons']))
  forecasts.extend(bits((l+h*t)+se[(11+int(h))%2])for h in [0.,1.,2.,5.]);stats.extend(bits(v)for v in statistics_from_errors(cs['y'],e,pr))
  confidence.extend(bits(z*math.sqrt(sequential_total(v*v for v in e)/12*variance(a,b,h,39,g)))for h in [1.,2.,5.,10.])
 extra={k:Path(__file__).with_name(k).read_text()for k in ['ets_seasonal.py','ets_confidence.py']}
 for kind,pred in [('forecast',forecasts),('stat',stats),('confidence',confidence)]:freeze(name+'-'+kind,dict(cases=cases,predictions=pred,extra_sources=extra,parameters=parameters,method='Conditional on STAT1-3 and a hypothesised initial-state nullspace: late odd-phase perturbations [1,-2,1] leave initial seasons and the regression slope unchanged in the recovered n12/p2 basis map. Large positive trend separates postfit reset. This tests initialization/recursion/normalization together; no nonzero seasonal product extension follows without a forward proof.'))
 finish(name+'-forecast',cases,forecasts,'FORECAST.ETS',[0.,1.,2.,5.])
 finish(name+'-stat',cases,stats,'FORECAST.ETS.STAT',[4,5,6,7])
 finish(name+'-confidence',cases,confidence,'FORECAST.ETS.CONFINT',[(h,.95)for h in [1.,2.,5.,10.]])

def gamma_map():
 cases=[]
 # A large linear trend separates the already measured postfit beta reset.
 # Reuse complete discovery designs rather than adding points to fit losses.
 for filename in ['step7-seasonal-stat.json','step7-basis-stat.json']:
  for original in read(OUT/filename)['cases']:
   z=dict(original);z['y']=[v+128.*j for j,v in enumerate(z['y'])];z['origin']=filename;cases.append(z)
 # One-dimensional projective map of the known n12/p2 initial-state nullspace.
 # A=1, B=r, C=-(3+7r)/11 keeps the recovered regression slope unchanged.
 for numerator in range(-48,49):
  ratio=numerator/4;changes=[1.,ratio,-(3.+7.*ratio)/11.]
  y=[(-2. if j%2==0 else 2.)+128.*j for j in range(12)]
  for j,delta in zip([7,9,11],changes):y[j]+=delta
  cases.append(dict(y=y,season=2,kind='gamma_projective_map',ratio=ratio,changes=changes))
 write(OUT/'step7-gamma-map-design.json',dict(cases=cases,producer_source=Path(__file__).read_text(),limits='Discovery map; ratio cases use rounded C and do not claim exact initial-state cancellation. Parameters are observations, not a fitted search program.'))
 c,s=service();d=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);write(OUT/'step7-gamma-map-answers.json',d)
 st=grouped(d);points=sorted(set(tuple(value(v[str(k)]['bits'])for k in [1,2,3])for v in st))
 write(OUT/'step7-gamma-map-points.json',dict(points=points,cases=len(cases),limits='Observed reported triples; resetting or rounded parameters need separate separation before identifying raw grid values.'))
 print('gamma map cases',len(cases),'observed triples',len(points),points,flush=True)

def gamma_fresh(structured=False):
 from .ets_proof import freeze,finish
 rng=random.Random(397673 if structured else 397671);cases=[];pred=[];objectives=[]
 if structured:
  # A simple handwritten family suggested by all mapped axes: inner quarter
  # steps at gamma=.001, the complementary outer edge, a near-zero alpha
  # edge, and the high-alpha beta axis. Complete every Cartesian pair before
  # asking fresh questions; do not append the individual failed selections.
  grid=[[a,b,.001]for a in [.751,.501,.251,.126]for b in [.001,a-.001]]
  grid += [[.002,.001,g]for g in [.125,.25,.5,.75]]
  grid += [[a,b,.999-a]for a in [.1,.5,.75,.9]for b in [.001,a-.001]]
  grid += [[.998,b,.001]for b in [.099,.249,.499,.749,.899]]
 else:grid=read(OUT/'step7-gamma-map-points.json')['points']
 # Measured zero-objective starting point is first; all other points keep a
 # fixed declaration order. No raw triple is added after this proof freezes.
 start=[.751,.001,.001];grid=[start]+[p for p in grid if p!=start]
 for i in range(100):
  ex=rng.randrange(-200,201);offset=rng.randrange(-64,65)/8;amplitude=rng.randrange(1,33)/8;trend=128.
  v=rng.randrange(-16,17)/16
  u=rng.randrange(-32,33)/16 if i<64 else -7*v+rng.choice([-1.,1.])*2.**-rng.randrange(1,13)
  y=[offset+(-amplitude if j%2==0 else amplitude)+trend*j for j in range(12)]
  for j,delta in [(7,u+7*v),(9,-2*u-3*v),(11,u)]:y[j]+=delta
  y=[math.ldexp(x,ex)for x in y];init=(math.ldexp(offset,ex),math.ldexp(trend,ex),[math.ldexp(-amplitude,ex),math.ldexp(amplitude,ex)])
  def score(par):
   a,b,g=par;l,t,se=init[0],init[1],init[2].copy();errors=[0.]
   for j,actual in enumerate(y[1:],1):
    q=j%2;e=actual-((l+t)+se[q]);errors.append(e)
    if j>=8:l+=t
    else:l=l+(t+a*e);t+=b*e;se[q]+=g*e
   return sequential_total(e*e for e in errors)
  scores=[score(p)for p in grid];choice=min(range(len(grid)),key=lambda j:scores[j]);pred.extend(bits(x)for x in grid[choice]);objectives.append(scores)
  cases.append(dict(y=y,season=2,kind='fresh_gamma_nullspace',initial=init,u=u,v=v))
 name='step7-gamma-structured-fresh'if structured else'step7-gamma-fresh'
 freeze(name,dict(cases=cases,predictions=pred,grid=grid,objectives=objectives,extra_source=Path(__file__).read_text(),method='Fixed 25-point Cartesian axis/edge family'if structured else'Fixed 21 observed triples',limits='Hybrid validation SSE, known n12/p2 initial-state nullspace. Fresh dyadic mixtures, exponents and near-pole directions. This proves only this bounded family if all parameters match; it does not identify general seasonal initialization or unobserved grid points. Full versus tail SSE is undetermined here because the prefix error term is parameter-independent. Grid values that no tested input separates from neighbors are undetermined, not recovered constants.'))
 finish(name,cases,pred,'FORECAST.ETS.STAT',[1,2,3])

def gamma_plateau():
 from .plateau_probe import probe,Builder,Native,build_library
 import numpy as np
 rng=random.Random(397691);g=Builder(16)
 op=lambda name,a,b=0:g.op(name,a,b,store=False)
 constants={k:g.constant(bits(v))for k,v in [('alpha',.002),('beta',.001),('gamma',.125)]}
 l,t,se=12,13,[14,15]
 for j in range(1,12):
  q=j%2;p=op('add',op('add',l,t),se[q]);e=op('sub',j,p)
  l=op('add',l,op('add',t,op('mul',constants['alpha'],e)))
  t=op('add',t,op('mul',constants['beta'],e));se[q]=op('add',se[q],op('mul',constants['gamma'],e))
 two=g.constant(bits(2.));center=op('div',op('add',se[0],se[1]),two)
 l=op('add',l,center);se=[op('sub',s,center)for s in se];out=op('add',op('add',l,t),se[0])
 shape=dict(graph=g.graph(out,out),constants=constants,mask=0,cw=0x123f,function='FORECAST.ETS',argument_inputs=list(range(12)))
 cases=[];inputs=[]
 for i in range(8192):
  ex=rng.randrange(-200,201);offset=rng.randrange(-128,129)/8;amplitude=rng.randrange(1,33)/8;trend=128.;u=rng.randrange(-1024,1025)/1024
  y=[offset+(-amplitude if j%2==0 else amplitude)+trend*j for j in range(12)]
  for j,delta in [(7,u),(9,-2*u),(11,u)]:y[j]+=delta
  y=[math.ldexp(v,ex)for v in y];init=[math.ldexp(offset,ex),math.ldexp(trend,ex),math.ldexp(-amplitude,ex),math.ldexp(amplitude,ex)]
  cases.append(dict(y=y,season=2,kind='plateau_nullspace',initial=init));inputs.append([int(bits(v),16)for v in y+init])
 library=OUT/'seasonal-plateau-native.so';build_library(library);native=Native(library);inputs=np.array(inputs,dtype=np.uint64)
 # Verify the diagnostic instruction tree against the Python binary64 tree.
 baseline=native.grid(shape['graph'],inputs,shape['cw'],shape['mask'])
 for i,z in enumerate(cases):
  l,t,se,e,pr=fit(z['y'],2,.002,.001,.125,(z['initial'][0],z['initial'][1],z['initial'][2:]))
  if int(bits((l+t)+se[0]),16)!=int(baseline[i,0]):raise ValueError('plateau instruction tree differs')
 constants_result=probe(shape,inputs,native,keep=4);selected={}
 index={tuple(int(v)for v in row):i for i,row in enumerate(inputs)}
 for z in constants_result:
  for neighbor in z['neighbors']:
   for example in neighbor['examples']:
    key=tuple(int(v,16)for v in example['inputs']);selected[index[key]]=cases[index[key]]
 chosen=list(selected.values());pred=[f'{int(baseline[i,0]):016x}'for i in selected]
 write(OUT/'step7-gamma-plateau-frozen.json',dict(shape=shape,constants=constants_result,cases=chosen,predictions=pred,search_cases=len(cases),producer_source=Path(__file__).read_text(),limits='Conditional recurrence constant separation on the measured initial-state nullspace. No separator within this finite search means undetermined; gamma search boundaries and general initialization are not inferred.'))
 if chosen:
  c,s=service();parameters=ask_cases(c,s,chosen,'FORECAST.ETS.STAT',[1,2,3]);answers=ask_cases(c,s,chosen,'FORECAST.ETS',[1.]);write(OUT/'step7-gamma-plateau-parameters.json',parameters);write(OUT/'step7-gamma-plateau-answers.json',answers)
  parameter_failures=[i for i,z in enumerate(grouped(parameters))if [z[str(k)]['bits']for k in [1,2,3]]!=[bits(v)for v in [.002,.001,.125]]]
  failures=[dict(case=i,predicted=p,actual=a)for i,(p,a)in enumerate(zip(pred,answers['answers']))if a.get('bits')!=p]
  write(OUT/'step7-gamma-plateau-proof.json',dict(cases=len(chosen),parameter_failures=parameter_failures,forecast_failures=failures,constants=constants_result));print('plateau selected',len(chosen),'parameter failures',len(parameter_failures),'forecast failures',len(failures),flush=True)

def initialization_basis():
 rng=random.Random(397797);cases=[]
 for n,p in [(8,2),(12,2),(16,2),(24,2),(32,2),(16,4),(24,4),(32,4),(32,8)]:
  pattern=[float(2*q-p+1)for q in range(p)];base=[128.*j+pattern[j%p]+rng.randrange(-8,9)/16 for j in range(n)]
  for at in [-1,*range(n)]:
   y=base.copy()
   if at>=0:y[at]+=1/16
   cases.append(dict(y=y,season=p,kind='positive_trend_initialization_basis',at=at))
 write(OUT/'step7-positive-basis-design.json',dict(cases=cases,producer_source=Path(__file__).read_text(),limits='Discovery derivatives only. Large linear trend separates postfit reset; random residual background avoids singular parameter boundaries where possible. Ill-conditioned inversions remain diagnostics.'))
 c,s=service();write(OUT/'step7-positive-basis-stat.json',ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]));write(OUT/'step7-positive-basis-forecast.json',ask_cases(c,s,cases,'FORECAST.ETS',range(10)))
 print('positive basis cases',len(cases),flush=True)

def moving_initial(y,p,ma_order='direct',adjustment='centered',level_order='centered'):
 n=len(y);width=2*p;count=min(3*p,n-width);offset=p//2-1
 if count<p:raise ValueError('short initialization is a separate stage')
 window=sequential_total(y[:width]);totals=[0.]*p;counts=[0]*p
 for j in range(count):
  if ma_order=='direct':window=sequential_total(y[j:j+width])
  elif j:window=(window-y[j-1])+y[j+width-1]if ma_order=='subtract_add'else window+(y[j+width-1]-y[j-1])
  delta=y[j+offset]-window/width;q=(j+offset)%p;totals[q]+=delta;counts[q]+=1
 raw=[v/c for v,c in zip(totals,counts)];center=sequential_total(raw)/p;se=[v-center for v in raw]
 adjusted=[v-se[j%p]for j,v in enumerate(y)]if adjustment=='centered'else[(v-raw[j%p])+center for j,v in enumerate(y)]
 t=initial_trend(adjusted);l=y[0]-se[0]if level_order=='centered'else(y[0]-raw[0])+center
 return l,t,se

def initialization_fresh(revised=False,partial=False):
 from .ets_proof import freeze,finish
 from .ets_confidence import variance
 rng=random.Random(397817 if partial else 397813 if revised else 397811);cases=[]
 for i in range(128):
  n=(9 if i<64 else [8,10,11,12][(i-64)%4])if partial else rng.randrange(8,33);ex=rng.randrange(-200,201);offset=3*rng.randrange(-256,257)/16;amplitude=3*rng.randrange(1,65)/16;trend=192.
  # Multiples of three make the m=3 phase division exact before the OLS
  # stage. This is an algebraic input region, fixed before reading answers.
  y=[offset+(-amplitude if j%2==0 else amplitude)+trend*j+3*rng.randrange(-16,17)/16 for j in range(n)]
  cases.append(dict(y=[math.ldexp(v,ex)for v in y],season=2,kind='fresh_period_two_moving_average_initialization'))
 name='step7-initialization-partial-fresh'if partial else'step7-initialization-window-fresh'if revised else'step7-initialization-fresh';c,s=service();parameters=ask_cases(c,s,cases,'FORECAST.ETS.STAT',[1,2,3]);write(OUT/(name+'-parameters.json'),parameters);st=grouped(parameters);forecasts=[];stats=[];confidence=[];z=value(read(OUT/'step6-normal-quantiles.json')['reply']['rows'][0]['bits'])
 for cs,ss in zip(cases,st):
  a,b,g=[value(ss[str(k)]['bits'])for k in [1,2,3]];l,t,se,e,pr=fit(cs['y'],2,a,b,g,moving_initial(cs['y'],2))
  forecasts.extend(bits((l+h*t)+se[(len(cs['y'])-1+int(h))%2])for h in [0.,1.,2.,5.]);stats.extend(bits(v)for v in statistics_from_errors(cs['y'],e,pr));confidence.extend(bits(z*math.sqrt(sequential_total(v*v for v in e)/len(cs['y'])*variance(a,b,h,39,g)))for h in [1.,2.,5.,10.])
 for kind,pred in [('forecast',forecasts),('stat',stats),('confidence',confidence)]:freeze(name+'-'+kind,dict(cases=cases,predictions=pred,parameters=parameters,extra_sources={k:Path(__file__).with_name(k).read_text()for k in ['ets_seasonal.py','ets_confidence.py']},method='Conditional on measured parameters, the fixed 2p-point forward moving average, min(3p,n-2p) observations with separate phase counts, observation offset p//2-1, centered phase means and OLS of adjusted whole history. Fresh period-two noisy histories and exponents, including odd lengths. No optimizer or reset decision is supplied by this proof.'))
 finish(name+'-forecast',cases,forecasts,'FORECAST.ETS',[0.,1.,2.,5.]);finish(name+'-stat',cases,stats,'FORECAST.ETS.STAT',[4,5,6,7]);finish(name+'-confidence',cases,confidence,'FORECAST.ETS.CONFINT',[(h,.95)for h in [1.,2.,5.,10.]])

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('command',choices=['design','analyze','basis','basis-analysis','zero','detection','zero-edges','automatic-two','detection-analysis','automatic-edges','nonzero','nonzero-mixed','gamma-map','gamma-fresh','gamma-structured','gamma-plateau','initialization-basis','initialization-analysis','initialization-fresh','initialization-window-fresh','initialization-partial-fresh']);a=p.parse_args();{'design':design,'analyze':analyze,'basis':basis,'basis-analysis':basis_analysis,'zero':zero_residual_proof,'detection':detection_map,'zero-edges':zero_edges,'automatic-two':automatic_two_proof,'detection-analysis':detection_analysis,'automatic-edges':automatic_lattice_edges,'nonzero':nonzero_controls,'nonzero-mixed':lambda:nonzero_controls(True),'gamma-map':gamma_map,'gamma-fresh':gamma_fresh,'gamma-structured':lambda:gamma_fresh(True),'gamma-plateau':gamma_plateau,'initialization-basis':initialization_basis,'initialization-analysis':lambda:basis_analysis('step7-positive-basis'),'initialization-fresh':initialization_fresh,'initialization-window-fresh':lambda:initialization_fresh(True),'initialization-partial-fresh':lambda:initialization_fresh(True,True)}[a.command]()
