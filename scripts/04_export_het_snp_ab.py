#!/usr/bin/env python3
"""Export all heterozygous biallelic SNP allele balances for selected samples from a prepared Hail MT.
One immutable BGZF TSV plus JSON provenance per chunk. No de novo or AC filtering.
"""
import argparse
import json
import os
import shutil
import uuid
from datetime import datetime, timezone

import boto3
import hail as hl
import hailtop.fs as hfs


def main():
    p = argparse.ArgumentParser()
    for name in ('chunk-id','matrix','prep-manifest','samples','out-prefix','tmp-base'):
        p.add_argument('--'+name, required=True)
    p.add_argument('--pass-only', action='store_true', help='Keep only VCF PASS/unfiltered variants')
    a = p.parse_args()
    selected = list(dict.fromkeys(s.strip() for s in a.samples.split(',') if s.strip()))
    if len(selected) < 3:
        raise ValueError('Provide at least 3 sample IDs')
    output = a.out_prefix + '.het_snp_ab.tsv.bgz'
    marker = a.out_prefix + '.het_snp_ab.completed.json'
    if hfs.exists(output) or hfs.exists(marker):
        raise FileExistsError(f'Output already exists for {a.chunk_id}; do not overwrite: {output}')
    with hfs.open(a.prep_manifest) as f:
        prep = json.load(f)
    if prep.get('status') != 'complete' or prep.get('matrix') != a.matrix:
        raise ValueError('Prepared matrix provenance check failed')
    if a.chunk_id not in a.matrix:
        raise ValueError('Chunk ID does not match MT path')

    run_id = uuid.uuid4().hex
    local_tmp = '/tmp/ab-' + run_id
    shared_tmp = a.tmp_base.rstrip('/') + '/ab-' + run_id
    os.makedirs(local_tmp, exist_ok=False)
    hl.init(local='local[8]', tmp_dir=shared_tmp, local_tmpdir=local_tmp,
            app_name='hail-het-ab-' + a.chunk_id,
            spark_conf={'spark.hadoop.fs.s3.impl': 'org.apache.hadoop.fs.s3a.S3AFileSystem',
                        'spark.hadoop.fs.s3a.impl': 'org.apache.hadoop.fs.s3a.S3AFileSystem',
                        'spark.hadoop.fs.s3a.connection.maximum': '128'})
    try:
        mt = hl.read_matrix_table(a.matrix)
        required = {'GT','AD','DP','GQ'}
        missing = required - set(mt.entry.dtype.fields)
        if missing:
            raise ValueError(f'Missing required entry fields: {sorted(missing)}')
        if list(mt.col_key.dtype.fields) != ['s']:
            raise ValueError('Expected column key s')
        mt = mt.filter_cols(hl.literal(set(selected)).contains(mt.s))
        found = set(mt.aggregate_cols(hl.agg.collect(mt.s)))
        if found != set(selected):
            raise ValueError(f'Samples missing from MT: {sorted(set(selected)-found)}')
        mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]))
        if a.pass_only:
            mt = mt.filter_rows(hl.is_missing(mt.filters) | (hl.len(mt.filters) == 0))
        mt = mt.filter_entries(
            hl.is_defined(mt.GT) & mt.GT.is_het() &
            hl.is_defined(mt.AD) & hl.or_else(hl.len(mt.AD) == 2, False) &
            hl.is_defined(mt.DP) & (mt.DP >= 20) &
            hl.is_defined(mt.GQ) & (mt.GQ >= 30) &
            hl.or_else(hl.sum(mt.AD) > 0, False)
        )
        mt = mt.annotate_entries(AB=hl.float64(mt.AD[1]) / hl.sum(mt.AD))
        ht = mt.entries()
        ht = ht.filter(hl.is_defined(ht.AB))
        ht = ht.select(sample_id=ht.s, ref=ht.alleles[0], alt=ht.alleles[1],
                       DP=ht.DP, GQ=ht.GQ, AD_ref=ht.AD[0], AD_alt=ht.AD[1], AB=ht.AB,
                       contig=ht.locus.contig, pos=ht.locus.position)
        ht = ht.key_by()
        ht = ht.select('sample_id','contig','pos','ref','alt','DP','GQ','AD_ref','AD_alt','AB')
        # Keep the computation small, only selected columns. Export is distributed, not to_pandas().
        ht.export(output)
        bucket_key = marker[5:].split('/', 1)
        if len(bucket_key) != 2:
            raise ValueError('Expected S3 output prefix')
        boto3.client('s3',region_name='us-east-1').put_object(
            Bucket=bucket_key[0], Key=bucket_key[1],
            Body=json.dumps({'status':'complete','chunk_id':a.chunk_id, 'matrix':a.matrix,
                'samples':selected,'pass_only':a.pass_only,'min_dp':20,'min_gq':30,
                'output':output,'finished_utc':datetime.now(timezone.utc).isoformat()},indent=2).encode(),
            ContentType='application/json')
        print(f'Completed {a.chunk_id}: {output}',flush=True)
    finally:
        hl.stop()
        shutil.rmtree(local_tmp, ignore_errors=True)
        try:
            hfs.rmtree(shared_tmp)
        except Exception as err:
            print(f'WARNING: scratch cleanup failed: {err}',flush=True)

if __name__ == '__main__':
    main()
