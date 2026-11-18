#!/usr/bin/env python3
"""Read-only chr22 native Hail de novo versus original Mendelian-screen comparison.
Uses boto3 and the local chunks_chr22_v1.tsv; writes TSV summaries locally only.
"""
import csv
import gzip
import io
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import urlparse
import boto3

BUCKET='ucsf-spark'
ROOT='01-sfari/analysis/trio-wgs/hail-mendelian/results/'
OLD=ROOT+'chr22_nvme_v1/'
NEW=ROOT+'chr22_hail_native_v1/'
OUT=Path('chr22_native_comparison')
OUT.mkdir(exist_ok=True)
s3=boto3.client('s3', region_name='us-east-1')

def get_bytes(key):
    return s3.get_object(Bucket=BUCKET,Key=key)['Body'].read()

def read_json(key):
    return json.loads(get_bytes(key))

def rows_from_bgz(key):
    raw=get_bytes(key)
    with gzip.open(io.BytesIO(raw), 'rt') as fh:
        return list(csv.DictReader(fh,delimiter='\t'))

def require_keys(row, fields, context):
    missing=set(fields)-set(row)
    if missing: raise ValueError(f'{context}: missing fields: {sorted(missing)}')

def key(r):
    return (r['contig'],int(r['pos']),r['ref'],r['alt'],r['child_id'])

def locus(k): return k[:4]

def snv(k): return len(k[2])==len(k[3])==1

def as_float(value):
    try: return float(value)
    except (TypeError,ValueError): return None

def write_tsv(path,rows,fields):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

with open('chunks_chr22_v1.tsv',newline='') as f:
    manifest=list(csv.DictReader(f,delimiter='\t'))
if len(manifest)!=17: raise ValueError(f'Expected 17 manifest chunks, found {len(manifest)}')
if len({r['chunk_id'] for r in manifest})!=17: raise ValueError('Duplicate chunk IDs')

old_events={}; new_events={}; overview=[]; confidence_total=Counter();
for c in manifest:
    chunk=c['chunk_id']; print('Loading',chunk,flush=True)
    old_meta=read_json(OLD+chunk+'.completed.json')
    new_meta=read_json(NEW+chunk+'.hail_denovo.completed.json')
    if old_meta.get('status')!='complete' or new_meta.get('status')!='complete':
        raise ValueError('Incomplete run manifest for '+chunk)
    if new_meta.get('chunk_id')!=chunk or new_meta.get('matrix')!=c['matrix']:
        raise ValueError('Native matrix provenance mismatch: '+chunk)
    old_rows=rows_from_bgz(OLD+chunk+'.candidate_trios.tsv.bgz')
    new_rows=rows_from_bgz(NEW+chunk+'.hail_denovo.tsv.bgz')
    old_need=['contig','pos','ref','alt','child_id','child_AB','PREQC_AC']
    new_need=['contig','pos','ref','alt','child_id','p_de_novo','confidence','prior','child_GT','father_GT','mother_GT']
    for r in old_rows: require_keys(r,old_need,'original/'+chunk)
    for r in new_rows: require_keys(r,new_need,'native/'+chunk)
    old_chunk={key(r):r for r in old_rows}; new_chunk={key(r):r for r in new_rows}
    if len(old_chunk)!=len(old_rows): raise ValueError('Duplicate original keys: '+chunk)
    if len(new_chunk)!=len(new_rows): raise ValueError('Duplicate native keys: '+chunk)
    if old_meta['candidate_trio_events']!=len(old_rows): raise ValueError('Original event count mismatch: '+chunk)
    if old_meta['candidate_sites']!=len({locus(k) for k in old_chunk}): raise ValueError('Original sites mismatch: '+chunk)
    if new_meta['events']!=len(new_rows) or new_meta['sites']!=len({locus(k) for k in new_chunk}):
        raise ValueError('Native event/site count mismatch: '+chunk)
    conf=Counter(r['confidence'] for r in new_rows)
    if dict(conf)!={str(k):int(v) for k,v in new_meta['confidence'].items()}:
        raise ValueError('Native confidence count mismatch: '+chunk)
    confidence_total.update(conf)
    overlap=set(old_chunk)&set(new_chunk)
    overview.append(dict(chunk_id=chunk,old_events=len(old_rows),new_events=len(new_rows),shared_events=len(overlap),old_only=len(old_rows)-len(overlap),native_only=len(new_rows)-len(overlap),native_high=conf.get('HIGH',0),native_medium=conf.get('MEDIUM',0),native_low=conf.get('LOW',0)))
    old_events.update(old_chunk);new_events.update(new_chunk)

old_keys=set(old_events);native_keys=set(new_events); shared=old_keys&native_keys
old_only=old_keys-native_keys;native_only=native_keys-old_keys
assert len(old_events)==72261, f'Original events expected 72261 got {len(old_events)}'
assert len({locus(k) for k in old_keys})==60403

print('\n=== CHROMOSOME 22 OVERALL ===')
print('Original candidate trio events:',f'{len(old_keys):,}')
print('Native candidate trio events:',f'{len(native_keys):,}')
print('Shared exact variant-child events:',f'{len(shared):,}')
print('Original-only events:',f'{len(old_only):,}')
print('Native-only events:',f'{len(native_only):,}')
print('Original unique sites:',f'{len({locus(k) for k in old_keys}):,}')
print('Native unique sites:',f'{len({locus(k) for k in native_keys}):,}')
print('Shared unique sites:',f'{len({locus(k) for k in old_keys}&{locus(k) for k in native_keys}):,}')
print('Native confidence counts:',dict(confidence_total))
print('Native event classification SNV:',sum(snv(k) for k in native_keys),'non-SNV:',sum(not snv(k) for k in native_keys))

print('\n=== OVERLAP BY NATIVE CONFIDENCE ===')
for cls in sorted(confidence_total):
    cls_keys={k for k in native_keys if new_events[k]['confidence']==cls}
    overlap_cls=cls_keys&old_keys
    print(f'{cls:8}: total={len(cls_keys):,} shared={len(overlap_cls):,} native-only={len(cls_keys)-len(overlap_cls):,}')

print('\n=== OUTLIER CHILDREN ===')
outlier_rows=[]
for child in ['SDSM-V8W','SDSM-2FNA']:
    a={k for k in old_keys if k[4]==child};b={k for k in native_keys if k[4]==child}
    confidence=Counter(new_events[k]['confidence'] for k in b)
    row=dict(child_id=child,old_events=len(a),native_events=len(b),shared=len(a&b),old_only=len(a-b),native_only=len(b-a),high=confidence.get('HIGH',0),medium=confidence.get('MEDIUM',0),low=confidence.get('LOW',0))
    outlier_rows.append(row)
    print(row)

# Examine native confidence against low AB for exactly matching original events
print('\n=== MATCHED EVENTS: ALLELE BALANCE BY NATIVE CONFIDENCE ===')
for cls in sorted(confidence_total):
    ab=[as_float(old_events[k]['child_AB']) for k in shared if new_events[k]['confidence']==cls]
    ab=[x for x in ab if x is not None]
    if ab: print(cls,'n=',len(ab),'median_AB=',round(statistics.median(ab),3),'AB<0.30=',f'{sum(x<.3 for x in ab)/len(ab):.1%}')

per_child_old=Counter(k[4] for k in old_keys)
per_child_new=Counter(k[4] for k in native_keys)
children=sorted(set(per_child_old)|set(per_child_new))
child_rows=[]
for child in children:
    a={k for k in old_keys if k[4]==child};b={k for k in native_keys if k[4]==child}
    conf=Counter(new_events[k]['confidence'] for k in b)
    child_rows.append(dict(child_id=child,original=len(a),native=len(b),shared=len(a&b),native_high=conf.get('HIGH',0),native_medium=conf.get('MEDIUM',0),native_low=conf.get('LOW',0)))
write_tsv(OUT/'per_chunk.tsv',overview,list(overview[0]))
write_tsv(OUT/'per_child.tsv',child_rows,list(child_rows[0]))
write_tsv(OUT/'outlier_children.tsv',outlier_rows,list(outlier_rows[0]))
# Event-level classifications without redundant genotype data; no patient IDs beyond existing cohort IDs
comparison=[]
for k in sorted(old_keys|native_keys):
    a=old_events.get(k);b=new_events.get(k)
    comparison.append(dict(contig=k[0],pos=k[1],ref=k[2],alt=k[3],child_id=k[4],in_original=bool(a),in_native=bool(b),confidence=b['confidence'] if b else '',p_de_novo=b['p_de_novo'] if b else '',prior=b['prior'] if b else '',original_child_AB=a['child_AB'] if a else '',original_PREQC_AC=a['PREQC_AC'] if a else ''))
write_tsv(OUT/'event_overlap.tsv',comparison,list(comparison[0]))
print('\nSaved local comparison files to',OUT.resolve())
print('NOTE: Hail confidence categories and posterior probabilities do not independently validate calls.')
