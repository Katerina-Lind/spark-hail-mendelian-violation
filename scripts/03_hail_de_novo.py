#!/usr/bin/env python3
"""Run Hail's native de_novo caller on an existing immutable chr22 MatrixTable."""
import argparse
import csv
import json
import os
import shutil
import uuid
from datetime import datetime, timezone

import boto3
import hail as hl
import hailtop.fs as hfs


def pedigree_from_tsv(path):
    required = {'family_id', 'sample_id', 'father_sample_id', 'mother_sample_id'}
    trios = []
    seen = set()
    with hfs.open(path, 'r') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError('Unexpected pedigree TSV columns')
        for record in reader:
            child = record['sample_id'].strip()
            dad = record['father_sample_id'].strip()
            mom = record['mother_sample_id'].strip()
            fam = record['family_id'].strip()
            if not all((child, dad, mom, fam)) or child in (dad, mom) or dad == mom or child in seen:
                raise ValueError('Invalid or duplicate pedigree child: ' + child)
            seen.add(child)
            trios.append(hl.Trio(s=child, pat_id=dad, mat_id=mom,
                                 fam_id=fam, is_female=None))
    if len(trios) != 7027:
        raise ValueError(f'Expected 7027 complete trios, observed {len(trios)}')
    return hl.Pedigree(trios)


def main():
    parser = argparse.ArgumentParser()
    for flag in ('matrix', 'prep-manifest', 'pedigree', 'lcr-bed', 'chunk-id', 'out-prefix', 'tmp-base'):
        parser.add_argument('--' + flag, required=True)
    a = parser.parse_args()
    target = a.out_prefix + '.hail_denovo.tsv.bgz'
    marker = a.out_prefix + '.hail_denovo.completed.json'
    if hfs.exists(target) or hfs.exists(marker):
        raise FileExistsError('Immutable output exists; select a NEW output prefix')
    with hfs.open(a.prep_manifest, 'r') as f:
        prep = json.load(f)
    if prep.get('status') != 'complete' or prep.get('matrix') != a.matrix or prep.get('schema_version') != 'glnexus_split_autosomal_snv_indel_v1':
        raise ValueError('MatrixTable provenance does not match preparation manifest')
    if a.chunk_id not in a.matrix:
        raise ValueError('Chunk ID does not match matrix path')
    run_id = uuid.uuid4().hex
    local_tmp = f'/tmp/dnv-native-{run_id}'
    os.makedirs(local_tmp, exist_ok=False)
    shared_tmp = a.tmp_base.rstrip('/') + '/native-' + run_id
    hl.init(local='local[8]', tmp_dir=shared_tmp, local_tmpdir=local_tmp,
            app_name='hail-native-denovo-' + a.chunk_id,
            spark_conf={'spark.hadoop.fs.s3.impl': 'org.apache.hadoop.fs.s3a.S3AFileSystem',
                        'spark.hadoop.fs.s3a.impl': 'org.apache.hadoop.fs.s3a.S3AFileSystem',
                        'spark.hadoop.fs.s3a.connection.maximum': '128'})
    try:
        mt = hl.read_matrix_table(a.matrix)
        needed = {'GT', 'AD', 'DP', 'GQ', 'PL'}
        if not needed.issubset(set(mt.entry.dtype.fields)):
            raise ValueError('Missing fields required by Hail de_novo: ' + repr(needed - set(mt.entry.dtype.fields)))
        if list(mt.row_key.dtype.fields) != ['locus', 'alleles']:
            raise ValueError('Unexpected row key; native de_novo requires locus, alleles')
        if list(mt.col_key.dtype.fields) != ['s']:
            raise ValueError('Unexpected column key')
        ped = pedigree_from_tsv(a.pedigree)
        family_lookup = hl.literal({t.s: t.fam_id for t in ped.trios})
        # Apply the SAME site-level eligibility as the original screen; leave
        # entry-level quality decisions to native hl.de_novo.
        mt = hl.variant_qc(mt)
        mt = mt.filter_rows(mt.variant_qc.AC[1] < 20)
        mt = mt.drop('variant_qc')
        lcr = hl.import_bed(a.lcr_bed, reference_genome='GRCh38')
        mt = mt.filter_rows(~hl.is_defined(lcr[mt.locus]))
        mt = mt.filter_rows(hl.is_defined(mt.qual) & (mt.qual >= 150))
        mt = mt.filter_rows(hl.is_missing(mt.filters) | ~mt.filters.contains('MONOALLELIC'))
        # Zero external AF prior: native Hail will still use cohort allele
        # frequency and its global prior. Not equivalent to using gnomAD.
        ht = hl.de_novo(mt, ped, pop_frequency_prior=hl.float64(0.0))
        ht = ht.select(
            child_id=ht.id,
            father_id=ht.father.s,
            mother_id=ht.mother.s,
            family_id=family_lookup[ht.id],
            p_de_novo=ht.p_de_novo,
            confidence=ht.confidence,
            prior=ht.prior,
            child_GT=hl.str(ht.proband_entry.GT),
            father_GT=hl.str(ht.father_entry.GT),
            mother_GT=hl.str(ht.mother_entry.GT),
            child_DP=ht.proband_entry.DP,
            father_DP=ht.father_entry.DP,
            mother_DP=ht.mother_entry.DP,
            child_GQ=ht.proband_entry.GQ,
            father_GQ=ht.father_entry.GQ,
            mother_GQ=ht.mother_entry.GQ,
            child_AD=ht.proband_entry.AD,
            father_AD=ht.father_entry.AD,
            mother_AD=ht.mother_entry.AD,
        )
        # Avoid counting/recalculating the caller more than once.
        ht = ht.checkpoint(shared_tmp + '/native-calls.ht', overwrite=False)
        by_confidence = ht.aggregate(hl.agg.counter(ht.confidence))
        events = sum(by_confidence.values())
        sites = ht.aggregate(hl.agg.count_distinct(hl.tuple([ht.locus, ht.alleles])))
        out = ht.annotate(contig=ht.locus.contig, pos=ht.locus.position,
                          ref=ht.alleles[0], alt=ht.alleles[1])
        out = out.key_by().drop('locus', 'alleles')
        # Retain scalar fields for a robust TSV export.
        out = out.select('contig','pos','ref','alt','child_id','family_id',
                         'father_id','mother_id','p_de_novo','confidence','prior',
                         'child_GT','father_GT','mother_GT',
                         'child_DP','father_DP','mother_DP',
                         'child_GQ','father_GQ','mother_GQ',
                         child_AD=hl.str(out.child_AD),
                         father_AD=hl.str(out.father_AD),
                         mother_AD=hl.str(out.mother_AD))
        out.export(target)
        p = target.split('/', 3)
        if p[0] != 's3:':
            raise ValueError('Output target must be an S3 URI')
        boto3.client('s3').put_object(
            Bucket=p[2], Key=marker.split('/', 3)[3],
            Body=json.dumps({'status':'complete', 'chunk_id':a.chunk_id,
                'matrix':a.matrix, 'prep_manifest':a.prep_manifest,
                'out_tsv':target, 'events':events, 'sites':sites,
                'confidence':by_confidence, 'caller':'hl.de_novo',
                'external_af_prior':'constant_zero',
                'finished_utc':datetime.now(timezone.utc).isoformat()},
                indent=2).encode('utf-8'), ContentType='application/json')
        print(f'{a.chunk_id}: {events} calls, {sites} sites, {by_confidence}', flush=True)
    finally:
        hl.stop()
        shutil.rmtree(local_tmp, ignore_errors=True)
        try:
            hfs.rmtree(shared_tmp)
        except Exception as error:
            print(f"WARNING: could not remove run scratch {shared_tmp}: {error}")

if __name__ == '__main__':
    main()
