#!/usr/bin/env python3
"""Extract cohort-QC'ed, LD-pruned genotype markers from one existing Hail MT.

Intended for Hail 0.2.138 and AWS Batch. Does not modify source MatrixTables.
"""
import argparse
import gzip
import json
import os
import shutil
import uuid
from datetime import datetime, timezone

import boto3
import hail as hl
import hailtop.fs as hfs


def s3parts(uri):
    if not uri.startswith('s3://'):
        raise ValueError(f'Expected s3:// URI: {uri}')
    return uri[5:].split('/', 1)


def main():
    p = argparse.ArgumentParser()
    for field in ('chunk-id', 'chromosome', 'matrix', 'prep-manifest', 'samples', 'out-prefix', 'tmp-base'):
        p.add_argument('--' + field, required=True)
    p.add_argument('--window-modulus', type=int, default=25000)
    p.add_argument('--window-width', type=int, default=1000)
    p.add_argument('--ld-r2', type=float, default=0.2)
    p.add_argument('--min-maf', type=float, default=0.10)
    p.add_argument('--min-call-rate', type=float, default=0.98)
    a = p.parse_args()
    samples = list(dict.fromkeys(v.strip() for v in a.samples.split(',') if v.strip()))
    if len(samples) != 6:
        raise ValueError('Expected exactly six sample IDs')
    if a.window_modulus < 1 or not 0 < a.window_width <= a.window_modulus:
        raise ValueError('Invalid positional thinning window')
    if not 0 < a.ld_r2 < 1:
        raise ValueError('Invalid LD r2')
    if not 0 < a.min_maf < 0.5:
        raise ValueError('Invalid MAF')
    output = a.out_prefix + '.kinship_markers.tsv.bgz'
    marker = a.out_prefix + '.kinship_markers.completed.json'
    output_exists, marker_exists = hfs.exists(output), hfs.exists(marker)
    if marker_exists and output_exists:
        with hfs.open(marker) as f:
            prior = json.load(f)
        if (prior.get('status') == 'complete' and
                prior.get('chunk_id') == a.chunk_id and
                prior.get('chromosome') == a.chromosome and
                prior.get('matrix') == a.matrix and
                prior.get('samples') == samples and
                prior.get('output') == output and
                prior.get('filters') == {
                    'min_dp': 10, 'min_gq': 20,
                    'min_call_rate': a.min_call_rate,
                    'min_maf': a.min_maf,
                    'position_modulus': a.window_modulus,
                    'position_window_width': a.window_width,
                    'ld_r2': a.ld_r2,
                    'ld_pruning_window_bp': 500000,
                }):
            print(f'REUSE COMPLETED: {a.chunk_id}', flush=True)
            return
        raise ValueError('Existing completed output has different provenance or filters: ' + a.out_prefix)
    if marker_exists or output_exists:
        raise FileExistsError('Partial immutable output exists; inspect before retry: ' + a.out_prefix)
    with hfs.open(a.prep_manifest) as f:
        prep = json.load(f)
    if prep.get('status') != 'complete' or prep.get('matrix') != a.matrix:
        raise ValueError('Matrix provenance does not match completed preparation')
    if a.chunk_id not in a.matrix:
        raise ValueError('Matrix URI does not match expected chunk')
    run_id = uuid.uuid4().hex
    local_tmp = '/tmp/kin-' + run_id
    shared_tmp = a.tmp_base.rstrip('/') + '/kin-' + run_id
    os.makedirs(local_tmp, exist_ok=False)
    hl.init(local='local[8]', tmp_dir=shared_tmp, local_tmpdir=local_tmp,
            app_name='kinship-markers-' + a.chunk_id,
            spark_conf={'spark.hadoop.fs.s3.impl': 'org.apache.hadoop.fs.s3a.S3AFileSystem',
                        'spark.hadoop.fs.s3a.impl': 'org.apache.hadoop.fs.s3a.S3AFileSystem',
                        'spark.hadoop.fs.s3a.connection.maximum': '128'})
    try:
        mt = hl.read_matrix_table(a.matrix)
        if 'GT' not in mt.entry.dtype.fields or 's' not in mt.col_key.dtype.fields:
            raise ValueError('Matrix lacks GT or sample column s')
        # Efficient deterministic pre-thinning, prior to cohort-wide variant QC.
        mt = mt.filter_rows(
            (mt.locus.contig == a.chromosome) &
            (hl.len(mt.alleles) == 2) &
            hl.is_snp(mt.alleles[0], mt.alleles[1]) &
            (mt.locus.position % a.window_modulus < a.window_width)
        )
        # Entry filters applied to all 16,320 samples, NOT only the selected six.
        mt = mt.filter_entries(
            hl.is_defined(mt.GT) & mt.GT.is_diploid() &
            hl.is_defined(mt.DP) & (mt.DP >= 10) &
            hl.is_defined(mt.GQ) & (mt.GQ >= 20)
        )
        # An empty chunk is valid: Hail LD pruning cannot build a zero-row BlockMatrix.
        # Export the same header-only BGZF-compatible gzip table used by the summariser.
        # The final manifest is written only after the output is present.
        n_after_thinning = mt.count_rows()
        n_after_qc = 0
        n_after_prune = 0
        if n_after_thinning:
            mt = hl.variant_qc(mt)
            mt = mt.filter_rows(
                (mt.variant_qc.call_rate >= a.min_call_rate) &
                (mt.variant_qc.AF[1] >= a.min_maf) &
                (mt.variant_qc.AF[1] <= 1 - a.min_maf)
            )
            n_after_qc = mt.count_rows()
        if n_after_qc:
            # LD pruning uses all cohort samples. Pruning is independent per chunk.
            pruned = hl.ld_prune(mt.GT, r2=a.ld_r2, bp_window_size=500000)
            mt = mt.filter_rows(hl.is_defined(pruned[mt.row_key]))
            n_after_prune = mt.count_rows()

        mt = mt.filter_cols(hl.literal(set(samples)).contains(mt.s))
        found = set(mt.aggregate_cols(hl.agg.collect(mt.s)))
        if found != set(samples):
            raise ValueError('Missing sample IDs: ' + repr(sorted(set(samples) - found)))

        if n_after_prune:
            ht = mt.entries()
            ht = ht.filter(hl.is_defined(ht.GT))
            ht = ht.select(
                chromosome=ht.locus.contig,
                pos=ht.locus.position,
                ref=ht.alleles[0],
                alt=ht.alleles[1],
                sample_id=ht.s,
                dosage=ht.GT.n_alt_alleles(),
            )
            ht = ht.key_by()
            ht = ht.select('chromosome', 'pos', 'ref', 'alt', 'sample_id', 'dosage')
            ht.export(output)
        else:
            # No valid SNPs in this chunk. This is an expected data condition.
            # Python gzip creates a valid stream readable by the summariser.
            out_bucket, out_key = s3parts(output)
            header = b'chromosome\tpos\tref\talt\tsample_id\tdosage\n'
            boto3.client('s3', region_name='us-east-1').put_object(
                Bucket=out_bucket, Key=out_key,
                Body=gzip.compress(header), ContentType='application/gzip',
                IfNoneMatch='*'
            )
            print(f'EMPTY MARKER CHUNK: {a.chunk_id} '
                  f'(thinned={n_after_thinning}, qc={n_after_qc})', flush=True)
        bucket, key = s3parts(marker)
        boto3.client('s3', region_name='us-east-1').put_object(
            Bucket=bucket, Key=key, ContentType='application/json',
            Body=json.dumps({
                'status': 'complete', 'chunk_id': a.chunk_id,
                'chromosome': a.chromosome, 'matrix': a.matrix,
                'samples': samples, 'output': output,
                'markers_after_thinning': n_after_thinning,
                'markers_after_qc': n_after_qc,
                'markers_after_ld_pruning': n_after_prune,
                'filters': {'min_dp': 10, 'min_gq': 20,
                            'min_call_rate': a.min_call_rate,
                            'min_maf': a.min_maf,
                            'position_modulus': a.window_modulus,
                            'position_window_width': a.window_width,
                            'ld_r2': a.ld_r2,
                            'ld_pruning_window_bp': 500000},
                'finished_utc': datetime.now(timezone.utc).isoformat(),
            }, indent=2).encode())
        print(f'COMPLETE: {a.chunk_id} => {output}', flush=True)
    finally:
        hl.stop()
        shutil.rmtree(local_tmp, ignore_errors=True)
        try:
            hfs.rmtree(shared_tmp)
        except Exception as e:
            print('WARNING: scratch cleanup: ' + str(e), flush=True)


if __name__ == '__main__':
    main()
