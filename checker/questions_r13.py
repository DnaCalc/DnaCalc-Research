"""Independent preservation, input-derivation and prediction audit for R13."""
import argparse
from collections import Counter
import csv
from pathlib import Path
import hashlib
import json

from exact import BINARY64 as B, EXTENDED64 as E, binary64
from x87 import Exact, Instructions, bits64, decode, encode, require


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def records(path):return [json.loads(line)for line in Path(path).read_text().splitlines()]


def standalone_graph(function,args,reciprocal):
    """Independently construct the advertised worksheet-boundary hypothesis."""
    require(type(reciprocal)is bool,'reciprocal flag is not logical')
    require(len(args)==(2 if function=='POWER'else 1),'wrong standalone arity')
    nodes=[dict(op='input',a=i)for i in range(len(args))]
    one=len(nodes);nodes.append(dict(op='const',bits='3ff0000000000000'))
    def emit(kind,left=0,right=0,spill=True):
        n=dict(op=kind,a=left,b=right)
        if spill:n['site']=sum('site'in n for n in nodes)
        nodes.append(n);return len(nodes)-1
    if function=='LN':
        ln2=emit('fldln2',spill=False);result=emit('fyl2x',ln2,0)
    else:
        if function=='EXP':
            log2e=emit('fldl2e',spill=False);argument=emit('mul',0,log2e)
        else:
            require(function=='POWER','unknown standalone function')
            argument=emit('fyl2x',1,0)
        negative=bool(int(args[0 if function=='EXP'else 1],16)>>63)
        if reciprocal and negative:argument=emit('neg',argument,spill=False)
        integer=emit('frndint',argument,spill=False)
        fraction=emit('sub',argument,integer)
        reduced=emit('f2xm1',fraction)
        mantissa=emit('add',one,reduced)
        result=emit('fscale',mantissa,integer)
        if reciprocal and negative:result=emit('div',one,result)
    return dict(nodes=nodes,outputs=[result,result],sites=sum('site'in n for n in nodes))


def check(packet,trace,gaps,library,output):
    packet=Path(packet);root=Path(__file__).resolve().parents[1]
    manifest=read(packet/'manifest.json');batches=read(packet/'batches.json')
    require(manifest['new_excel_answers']==0,'a prepared packet cannot contain new Excel answers')
    files={str(p.relative_to(packet)):sha(p)for p in packet.rglob('*')if p.is_file()and p!=packet/'manifest.json'}
    require(files==manifest['files_sha256'],'packet file list or checksum differs')
    frozen={str(p):sha(p)for p in (Path(__file__),root/'checker/x87.py',root/'checker/exact.py',Path(trace),Path(gaps),Path(library),packet/'manifest.json')}
    counts=Counter();listing=[]
    for batch in batches:
        doc=read(packet/batch['path']);calls=doc['witnesses'];fn=batch['function']
        require(doc['function']==fn and len(calls)==batch['calls'],'batch function/count differs')
        seen=set();ids=set()
        for row in calls:
            require(set(row)<=set(('id','args','result_index','input_transport')),'an answer or unexpected field is in a request')
            key=(json.dumps(row['args'],sort_keys=True),json.dumps(row.get('result_index')),row.get('input_transport','cell_value2_bulk'))
            require(key not in seen and row['id']not in ids,'duplicate request in a batch')
            require(row.get('input_transport','cell_value2_bulk')in ('cell_value2_bulk','formula_literals'),'unknown transport')
            seen.add(key);ids.add(row['id']);counts[batch['origin']]+=1
            listing.append(dict(origin=batch['origin'],function=fn,**row))
    for name,n in (('round11',64),('structural-r11-v2',52),('round12-pmt-groups',64)):
        require(counts[name]==n,'unanswered packet obligation omitted')
        source=root/'packets'/name;copied=packet/'prior'/name
        src={str(p.relative_to(source)):sha(p)for p in source.rglob('*')if p.is_file()}
        dst={str(p.relative_to(copied)):sha(p)for p in copied.rglob('*')if p.is_file()}
        require(src==dst,'immutable original packet changed')
    with (packet/'READING-LIST.tsv').open(newline='')as file:
        tsv=list(csv.DictReader(file,delimiter='\t'))
    require(len(tsv)==len(listing)==manifest['calls'],'reading list does not cover every batch')
    for i,(a,b)in enumerate(zip(tsv,listing),1):
        require(a['ordinal']==str(i)and all(a[k]==b[k]for k in ('origin','function','id')),'reading list order/identity differs')
        require(a['transport']==b.get('input_transport','cell_value2_bulk')and json.loads(a['arguments_json'])==b['args'],'reading list changed arguments/transport')
    source={r['id']:r for r in records(trace)}
    gap_ids=[r['id']for r in records(gaps)]
    byid={r['id']:r for r in listing if r['origin']=='r13'}
    require(len(byid)==counts['r13'],'duplicate R13 ID')
    selection=read(packet/'selection-r13.json')
    pmt=[r for r in selection if r['function']=='PMT']
    require([r['source_id']for r in pmt]==gap_ids and len(pmt)==24,'original gap rechecks incomplete')
    exact=Exact(Instructions(library));checked=Counter();covered=set()
    for s in selection:
        request=byid[s['id']];original=source[s['source_id']];covered.add(s['id'])
        require(request['function']==s['function'],'selection function differs')
        args=request['args'];d=s['derivation'];r=binary64(original['args'][0])
        base=bits64(B.round(E.round(1+r)));n=int(original['args'][1],16)
        if s['function']=='PMT':require(args==original['args'],'changed original gap inputs')
        elif s['function']=='LN':require(args==[f'0x{base:016x}'],'wrong stored base LN input')
        elif s['function']=='POWER':
            require(args[0]==f'0x{base:016x}'and int(args[1],16)in (n,n^(1<<63)),'wrong POWER base/exponent')
        elif s['function']=='EXP':
            require(d['cw']in (0x123f,0x133f)and int(d['periods'],16)in (n,n^(1<<63)),'wrong exponent derivation')
            cw=d['cw'];zero=encode(0)
            rom=exact.instructions.run('fldln2',zero,zero,cw)
            log=exact.instructions.run('fyl2xp1',rom,encode(r),cw)
            t=exact.operation('mul',encode(binary64(d['periods'])),log,cw,False)
            require(t.hex()==d['resident_t'],'native exponent transcript differs')
            require(args==[f'0x{bits64(B.round(decode(t)),bool(t[9]&128)):016x}'],'EXP input does not match fixed native derivation')
        else:raise ValueError('unknown R13 audit function')
        checked['input_derivations']+=1
    require(covered==set(byid),'an R13 request has no selection/derivation')
    seen=set()
    for prediction in records(packet/'predictions/native-elementary.jsonl'):
        request=byid[prediction['id']];graph=prediction['graph'];cw=prediction['cw']
        require(request['function']==prediction['function']and prediction['function']in ('LN','POWER','EXP'),'unexpected standalone prediction')
        require(graph==standalone_graph(prediction['function'],request['args'],prediction['reciprocal']),'standalone function graph differs')
        key=(prediction['id'],prediction['reciprocal'],cw)
        require(key not in seen,'duplicate elementary hypothesis');seen.add(key)
        values=prediction['predictions_by_store_mask']
        require(len(values)==1<<graph['sites'],'missing standalone store mask')
        for mask,proposed in enumerate(values):
            raw=exact.graph(graph,[int(x,16)for x in request['args']],cw,mask)[:10]
            require(proposed==f'0x{bits64(B.round(decode(raw)),bool(raw[9]&128)):016x}','native standalone prediction differs')
            checked['standalone_predictions']+=1
    expected={(r['id'],reciprocal,cw)for r in byid.values()if r['function']!='PMT'
              for reciprocal in ((False,)if r['function']=='LN'else(False,True))for cw in (0x123f,0x133f)}
    require(seen==expected,'incomplete standalone hypothesis coverage')
    require(counts==manifest['origin_calls'],'manifest origin counts differ')
    require(frozen=={p:sha(p)for p in frozen},'sources/evidence changed')
    result=dict(complete=True,calls=len(listing),origin_calls=dict(counts),checked=dict(checked),
        native_microcode_replays=dict(exact.instructions.counts),inputs_sha256=frozen,
        packet_file_sha256=files,new_excel_answers=0,
        limits='Prepared questions and hypothesis predictions only. Original180 obligations preserved exactly, fixed derivations checked. No answer, semantic pass, surviving PMT model or CPU portability inferred.')
    Path(output).write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('packet','trace','gaps','library','output'):parser.add_argument(key)
    a=parser.parse_args();r=check(a.packet,a.trace,a.gaps,a.library,a.output)
    print(json.dumps({k:v for k,v in r.items()if k not in ('packet_file_sha256','inputs_sha256')}))
