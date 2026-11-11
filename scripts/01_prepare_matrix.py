#!/usr/bin/env python3
"""Prepare a durable, versioned Hail MatrixTable from one GLnexus pVCF chunk.

Biological transformations replicate the student's import/autosome/split/SNV-indel
logic. Do NOT change an existing version in place. Supply an immutable input ID
(e.g. S3 VersionId or SHA256 recorded separately), and a distinct output path
for each distinct input/preprocessing version.
"""
import argparse
from datetime import datetime, timezone
import json
import os
import re
import shutil
import sys
import uuid

import hail as hl
import hailtop.fs as hfs

# Kept verbatim from student script; assumes GT is first FORMAT field.
PARTIAL_GT_PATTERN = r'(?<=\t)(?:[0-9]+[/|]\.|\.[/|][0-9]+)(?=:)'
PREP_SCHEMA_VERSION = 'glnexus_split_autosomal_snv_indel_v1'


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--vcf', required=True)
    p.add_argument('--matrix', required=True, help='Permanent, version-specific .mt path')
    p.add_argument('--manifest', required=True, help='Permanent provenance JSON path')
    p.add_argument('--input-version', required=True, help='Immutable S3 VersionId/checksum or a frozen-input manifest ID')
    p.add_argument('--local-cores', type=int, default=None)
    p.add_argument('--tmp-dir', help='Shared Hail temporary base (prefer S3 for distributed Spark)')
    p.add_argument('--local-tmpdir', default=os.environ.get('TMPDIR', '/tmp'))
    return p.parse_args()


def validate_existing(a):
    mt_exists = hfs.exists(a.matrix)
    manifest_exists = hfs.exists(a.manifest)
    if not mt_exists and not manifest_exists:
        return False
    if mt_exists != manifest_exists:
        raise RuntimeError('Incomplete preparation: MatrixTable and manifest do not both exist. Repair manually, never silently reuse.')
    with hfs.open(a.manifest, 'r') as stream:
        record = json.load(stream)
    expected = {'input_vcf': a.vcf, 'input_version': a.input_version,
                'matrix': a.matrix, 'schema_version': PREP_SCHEMA_VERSION,
                'reference_genome': 'GRCh38', 'status': 'complete'}
    for key, val in expected.items():
        if record.get(key) != val:
            raise RuntimeError(f'Existing preparation has mismatched {key}: {record.get(key)!r} != {val!r}')
    # Hail will check the underlying metadata on read; no full count needed.
    _ = hl.read_matrix_table(a.matrix)
    print(f'REUSE existing prepared MatrixTable: {a.matrix}', flush=True)
    return True


def main():
    a = arguments()
    if a.matrix.rstrip('/') == a.manifest.rstrip('/'):
        raise ValueError('MatrixTable and manifest must have different paths')
    if not a.matrix.endswith('.mt'):
        raise ValueError('--matrix must end in .mt')
    if a.local_cores is not None and a.local_cores < 1:
        raise ValueError('--local-cores must be >= 1')
    local_dir = os.path.join(a.local_tmpdir, 'hail-prep-' + uuid.uuid4().hex)
    os.makedirs(local_dir, exist_ok=False)
    shared_dir = (a.tmp_dir.rstrip('/') + '/' + os.path.basename(local_dir)) if a.tmp_dir else local_dir
    initialised = False
    try:
        kw = {'tmp_dir': shared_dir, 'local_tmpdir': local_dir, 'app_name': 'prepare-glnexus-matrix'}
        if a.local_cores:
            kw['local'] = f'local[{a.local_cores}]'
        hl.init(**kw, spark_conf={
    "spark.hadoop.fs.s3.impl":
        "org.apache.hadoop.fs.s3a.S3AFileSystem",
    "spark.hadoop.fs.s3a.connection.maximum": "128",
    "spark.hadoop.fs.s3a.impl":
        "org.apache.hadoop.fs.s3a.S3AFileSystem",
})
        initialised = True
        if validate_existing(a):
            return
        # Refuse even a partial previous output; in multi-worker runs each chunk
        # must have one writer and a deterministic, unique versioned path.
        if hfs.exists(a.matrix) or hfs.exists(a.manifest):
            raise FileExistsError('Partial preparation already exists')
        mt = hl.import_vcf(
            a.vcf, reference_genome='GRCh38', force_bgz=True,
            array_elements_required=False,
            find_replace=(PARTIAL_GT_PATTERN, './.'), block_size=128)
        required = {'GT', 'DP', 'AD', 'GQ', 'PL'}
        absent = required - set(mt.entry.dtype.fields)
        if absent:
            raise ValueError(f'Missing expected FORMAT fields: {sorted(absent)}')
        if 'AQ' not in mt.info.dtype.fields:
            raise ValueError('Missing required GLnexus INFO/AQ')
        # Count input records only once in the preparation stage.
        input_records = mt.count_rows()
        mt = mt.filter_rows(hl.literal({f'chr{i}' for i in range(1, 23)}).contains(mt.locus.contig))
        mt = mt.annotate_entries(_source_gq=mt.GQ)
        mt = hl.split_multi_hts(mt, keep_star=False)
        mt = mt.annotate_entries(GQ=mt._source_gq).drop('_source_gq')
        mt = mt.filter_rows(hl.is_snp(mt.alleles[0], mt.alleles[1]) |
                            hl.is_indel(mt.alleles[0], mt.alleles[1]))
        # Persist a raw-QC-independent dataset. No genotype filtering occurs here.
        # Do not set overwrite=True: a failed concurrent write must be visible.
        mt.write(a.matrix, overwrite=False)
        # Validate readability before marking dataset complete.
        saved = hl.read_matrix_table(a.matrix)
        cols = saved.count_cols()
        with hfs.open(a.manifest, 'w') as stream:
            json.dump({
                'status': 'complete', 'schema_version': PREP_SCHEMA_VERSION,
                'input_vcf': a.vcf, 'input_version': a.input_version,
                'matrix': a.matrix, 'reference_genome': 'GRCh38',
                'import_records': input_records, 'input_samples': cols,
                'hail_version': hl.__version__,
                'finished_utc': datetime.now(timezone.utc).isoformat(),
            }, stream, indent=2, sort_keys=True)
        print(f'CREATED MatrixTable: {a.matrix} ({input_records} input records, {cols} samples)', flush=True)
    finally:
        if initialised:
            hl.stop()
        shutil.rmtree(local_dir, ignore_errors=True)
        if shared_dir != local_dir:
            try:
                hfs.rmtree(shared_dir)
            except Exception as exc:
                print(f'WARNING: could not remove shared scratch: {exc}', file=sys.stderr)


if __name__ == '__main__':
    main()
