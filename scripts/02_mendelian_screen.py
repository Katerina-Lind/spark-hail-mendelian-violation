#!/usr/bin/env python3
"""Screen a PREPARED GLnexus MatrixTable for autosomal Mendelian code-2 candidates.

Inputs:
  --matrix      Preprocessed, persistent Hail MatrixTable (.mt)
  --pedigree    tab-delimited pedigrees.tsv containing sample_id, father_sample_id,
                mother_sample_id and family_id
  --lcr-bed     GRCh38 low-complexity-region BED (required unless --skip-lcr)
  --prep-manifest  provenance JSON created by 01_prepare_matrix.py
  --out-prefix  unique output prefix per chunk (local or Hail-supported URI)

Outputs:
  OUT.candidate_trios.tsv.bgz   One row per allele / qualifying trio.
  OUT.candidate_sites.vcf.bgz   Optional biallelic candidate VCF (omit with --no-vcf).
  OUT.completed.json            Publication-last success marker with counts and provenance.

The eligible genotype pattern is father=0/0, mother=0/0, child=0/1,
corresponding to Hail Mendelian error code 2 on autosomes. This is a
candidate screen, not proof of a de novo mutation.

The old GATK-specific SOR, QD, MQ, ReadPosRankSum and VQSR-PASS filters
are intentionally NOT applied because they cannot be recreated from the
GLnexus pVCF. QUAL >=150 is retained literally but the QUAL score in GLnexus
is not necessarily calibrated like the old GATK callset's QUAL.
"""

import argparse
import csv
from datetime import datetime, timezone
import json
import os
import re
import shutil
import sys
from urllib.parse import unquote, urlparse
import uuid

import hail as hl
import hailtop.fs as hfs

# This substitution is applied to VCF text BEFORE parsing with Hail.
# Only the GT portion of a sample field (after a TAB, before a colon) can
# match. It replaces e.g. ./1, 1/., .|2, 2|. with ./., avoiding the
# undefined semantics of GLnexus half-calls in Hail's diploid call type.
# The original input VCF is never modified. This particular pVCF has GT
# as its first FORMAT tag and other FORMAT tags after GT.
PARTIAL_GT_PATTERN = r'(?<=\t)(?:[0-9]+[/|]\.|\.[/|][0-9]+)(?=:)'


def cli():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix', required=True, help='Persistent .mt created by 01_prepare_matrix.py')
    ap.add_argument('--prep-manifest', required=True, help='Matching .mt provenance JSON')
    ap.add_argument('--pedigree', required=True, help='pedigrees.tsv, local or cloud URI')
    ap.add_argument('--out-prefix', required=True, help='Unique output prefix for this chunk (local or cloud URI)')
    ap.add_argument('--chunk-id', help='Optional filename-safe ID appended to --out-prefix; useful in Nextflow')
    lcr = ap.add_mutually_exclusive_group(required=True)
    lcr.add_argument('--lcr-bed', help='Same GRCh38 LCR BED used by the 2021 pipeline')
    lcr.add_argument('--skip-lcr', action='store_true', help='Explicitly omit original LCR exclusion')
    ap.add_argument('--local-cores', type=int, default=None,
                    help='Spark local cores, e.g. 4; otherwise use HAIL_LOCAL_CORES or Hail default')
    ap.add_argument('--tmp-dir', help='Base shared Hail temp directory (local or Hail-supported cloud URI)')
    ap.add_argument('--local-tmpdir', help='Base local scratch directory, e.g. /scratch (default: TMPDIR or /tmp)')
    ap.add_argument('--checkpoint-dir', help='Base Hail checkpoint directory; defaults to run-specific shared tmp')
    ap.add_argument('--no-checkpoints', action='store_true',
                    help='Disable intermediate final-site and candidate-table checkpoints (not advised at scale)')
    ap.add_argument('--keep-tmp', action='store_true', help='Retain run-specific temp/checkpoint data for debugging')
    ap.add_argument('--min-qual', type=float, default=150.0)
    ap.add_argument('--max-ac-exclusive', type=int, default=20)
    ap.add_argument('--no-vcf', action='store_true', help='Export primary candidate TSV only; skip large VCF')
    ap.add_argument('--overwrite', action='store_true', help='Replace existing outputs for this unique chunk prefix')
    return ap.parse_args()


def _local_path(uri):
    """Return local filesystem path or None for remote Hail URI."""
    parsed = urlparse(uri)
    if parsed.scheme == 'file':
        if parsed.netloc not in ('', 'localhost'):
            raise ValueError(f'Nonlocal file:// authority not supported: {uri}')
        return unquote(parsed.path)
    if not parsed.scheme:
        return uri
    return None


def _append_path(base, name):
    """Join a run-unique child path without changing a cloud URI scheme."""
    return base.rstrip('/') + '/' + name


def _publish_file(stage_path, final_path):
    """Local rename is atomic; cloud copy is NOT atomic (marker is published last)."""
    stage_local = _local_path(stage_path)
    final_local = _local_path(final_path)
    if stage_local is not None and final_local is not None:
        os.replace(stage_local, final_local)
    else:
        hfs.copy(stage_path, final_path)
        hfs.remove(stage_path)


def _cleanup_stage(stage_paths):
    for path in stage_paths:
        try:
            hfs.remove(path)
        except Exception as exc:
            print(f'WARNING: could not clean staging artifact {path}: {exc}', file=sys.stderr)


def _publish_results(staged, manifest_stage, manifest_final, overwrite, obsolete_outputs=()):
    """Publish data first, then a success marker. Never mark partial outputs as complete.

    For S3, copies are not transactional. Consumers must verify completed.json.
    Existing outputs with --overwrite may be replaced; use distinct prefixes
    for concurrently running chunks/retries (Nextflow work dirs are ideal).
    """
    targets = [final for _, final in staged] + list(obsolete_outputs) + [manifest_final]
    existing = [path for path in targets if hfs.exists(path)]
    if existing and not overwrite:
        raise FileExistsError(f'Existing output(s): {existing}; use unique --out-prefix or --overwrite')
    try:
        # Invalidate any old success marker before modifying output files.
        if hfs.exists(manifest_final):
            hfs.remove(manifest_final)
        for obsolete in obsolete_outputs:
            hfs.remove(obsolete)
        for stage, final in staged:
            _publish_file(stage, final)
        # This marker is the last file, after every export and metadata check.
        _publish_file(manifest_stage, manifest_final)
    except Exception:
        # Prevent partly published data from being mistaken for successful work.
        # On object storage, a failed copy may have left an incomplete object.
        for path in targets:
            try:
                hfs.remove(path)
            except Exception as exc:
                print(f'WARNING: unable to clean incomplete output {path}: {exc}', file=sys.stderr)
        raise


def _heap_and_spark_summary():
    sc = hl.spark_context()
    runtime = sc._jvm.java.lang.Runtime.getRuntime()
    heap_gb = round(runtime.maxMemory() / (1024 ** 3), 2)
    print(f'Spark master: {sc.master}; maximum Java heap: {heap_gb:.2f} GiB', flush=True)
    print(f'Spark configured driver memory: {sc.getConf().get("spark.driver.memory", "unset")}', flush=True)
    if heap_gb < 4:
        print('WARNING: JVM heap is <4 GiB. Set PYSPARK_SUBMIT_ARGS before starting Python.',
              file=sys.stderr, flush=True)
    return heap_gb


def _prepare_run(args):
    """Select unique outputs and isolated per-run Hail scratch/checkpoint paths."""
    chunk_id = args.chunk_id
    if chunk_id and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', chunk_id):
        raise ValueError('--chunk-id must contain only letters, digits, dot, underscore or dash')
    if args.local_cores is None and os.environ.get('HAIL_LOCAL_CORES'):
        args.local_cores = int(os.environ['HAIL_LOCAL_CORES'])
    if args.local_cores is not None and args.local_cores < 1:
        raise ValueError('--local-cores must be a positive integer')

    args.final_prefix = args.out_prefix + ('.' + chunk_id if chunk_id else '')
    args.run_id = uuid.uuid4().hex
    args.stage_prefix = args.final_prefix + '.__stage_' + args.run_id
    args.final_tsv = args.final_prefix + '.candidate_trios.tsv.bgz'
    args.final_vcf = args.final_prefix + '.candidate_sites.vcf.bgz'
    args.final_manifest = args.final_prefix + '.completed.json'
    args.stage_tsv = args.stage_prefix + '.candidate_trios.tsv.bgz'
    args.stage_vcf = args.stage_prefix + '.candidate_sites.vcf.bgz'
    args.stage_manifest = args.stage_prefix + '.completed.json'

    root = args.local_tmpdir or os.environ.get('TMPDIR') or '/tmp'
    if _local_path(root) is None:
        raise ValueError('--local-tmpdir must be a local directory, not s3:// or gs://')
    args.local_run_tmp = os.path.join(os.path.abspath(_local_path(root)), 'dnv-' + args.run_id)
    os.makedirs(args.local_run_tmp, exist_ok=False)
    args.shared_run_tmp = _append_path(args.tmp_dir, 'dnv-' + args.run_id) if args.tmp_dir else args.local_run_tmp
    if args.shared_run_tmp != args.local_run_tmp:
        hfs.mkdir(args.shared_run_tmp)
    args.checkpoint_run_dir = (_append_path(args.checkpoint_dir, 'dnv-' + args.run_id)
                               if args.checkpoint_dir else _append_path(args.shared_run_tmp, 'checkpoints'))
    if args.checkpoint_run_dir != args.shared_run_tmp:
        hfs.mkdir(args.checkpoint_run_dir)
    out_dir = os.path.dirname(args.final_prefix) or '.'
    if not hfs.exists(out_dir):
        hfs.mkdir(out_dir)

    # Check before ingesting a large VCF, not after expensive Hail calculations.
    # Even with --no-vcf, reject/reconcile a stale VCF left by a previous run.
    target_files = [args.final_tsv, args.final_manifest, args.final_vcf]
    existing = [p for p in target_files if hfs.exists(p)]
    if existing and not args.overwrite:
        raise FileExistsError(f'Output already exists: {existing}; use unique --out-prefix or --overwrite')
    return target_files


def read_pedigree(path):
    required = {'family_id', 'sample_id', 'father_sample_id', 'mother_sample_id'}
    output = []
    seen = set()
    with hfs.open(path, 'r') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError('Pedigree must contain columns: ' + ', '.join(sorted(required)))
        for i, row in enumerate(reader, start=2):
            kid = row['sample_id'].strip()
            dad = row['father_sample_id'].strip()
            mom = row['mother_sample_id'].strip()
            fam = row['family_id'].strip()
            if not all((kid, dad, mom, fam)):
                raise ValueError(f'Incomplete trio in pedigree line {i}')
            if kid in (dad, mom) or dad == mom:
                raise ValueError(f'Invalid relationship for {kid} in pedigree line {i}')
            if kid in seen:
                raise ValueError(f'Duplicate child sample_id: {kid}')
            seen.add(kid)
            output.append(hl.Trio(s=kid, pat_id=dad, mat_id=mom,
                                  fam_id=fam, is_female=None))
    if not output:
        raise ValueError('No trios in pedigree')
    return output


def valid_ab(entry):
    """ALT / (REF+ALT); missing if depths are not a valid positive pair."""
    ad = entry.AD
    good = (hl.is_defined(ad) & (hl.len(ad) == 2)
            & hl.is_defined(ad[0]) & hl.is_defined(ad[1])
            & (ad[0] >= 0) & (ad[1] >= 0)
            & (hl.sum(ad) > 0))
    return hl.or_missing(good, hl.float64(ad[1]) / hl.float64(hl.sum(ad)))


def good_depth(entry):
    return hl.is_defined(entry.DP) & (entry.DP >= 10) & (entry.DP <= 200)


def main():
    args = cli()
    hail_started = False
    stage_paths = []
    try:
        _prepare_run(args)
        stage_paths = [args.stage_tsv, args.stage_vcf, args.stage_manifest]
        init_kwargs = {'tmp_dir': args.shared_run_tmp, 'local_tmpdir': args.local_run_tmp,
                       'app_name': f'glnexus-mv-{args.chunk_id or "chunk"}'}
        if args.local_cores is not None:
            init_kwargs['local'] = f'local[{args.local_cores}]'
        hl.init(**init_kwargs, spark_conf={
    "spark.hadoop.fs.s3.impl":
        "org.apache.hadoop.fs.s3a.S3AFileSystem",
    "spark.hadoop.fs.s3a.impl":
        "org.apache.hadoop.fs.s3a.S3AFileSystem",
})
        hail_started = True
        heap_gb = _heap_and_spark_summary()
        print(f'Run ID: {args.run_id}; output prefix: {args.final_prefix}', flush=True)
        print(f'Hail shared temporary dir: {args.shared_run_tmp}', flush=True)
        print(f'Local scratch dir: {args.local_run_tmp}', flush=True)
        print(f'Checkpoint dir: {args.checkpoint_run_dir}', flush=True)
        _process_chunk(args, heap_gb)
    finally:
        _cleanup_stage(stage_paths)
        if hail_started:
            hl.stop()
        if 'args' in locals() and hasattr(args, 'local_run_tmp') and not args.keep_tmp:
            shutil.rmtree(args.local_run_tmp, ignore_errors=True)
            if getattr(args, 'shared_run_tmp', args.local_run_tmp) != args.local_run_tmp:
                try:
                    hfs.rmtree(args.shared_run_tmp)
                except Exception as exc:
                    print(f'WARNING: failed to remove run scratch {args.shared_run_tmp}: {exc}', file=sys.stderr)
            if (args.checkpoint_dir and hasattr(args, 'checkpoint_run_dir')
                    and not args.checkpoint_run_dir.startswith(args.shared_run_tmp + '/')):
                try:
                    hfs.rmtree(args.checkpoint_run_dir)
                except Exception as exc:
                    print(f'WARNING: failed to remove checkpoints {args.checkpoint_run_dir}: {exc}', file=sys.stderr)


def _process_chunk(args, heap_gb):
    """All biological thresholds and trio criteria below are unchanged."""

    trios = read_pedigree(args.pedigree)
    ped = hl.Pedigree(trios)
    print(f'Loaded {len(trios):,} complete trios from pedigree.', flush=True)

    # Reuse durable preprocessing: never parse the original pVCF in this script.
    with hfs.open(args.prep_manifest, 'r') as handle:
        prep = json.load(handle)
    if prep.get('status') != 'complete' or prep.get('matrix') != args.matrix:
        raise ValueError('Preparation manifest does not match requested MatrixTable')
    if prep.get('schema_version') != 'glnexus_split_autosomal_snv_indel_v1':
        raise ValueError('Incompatible preprocessing schema version')
    input_records = prep['import_records']
    mt = hl.read_matrix_table(args.matrix)
    required_entry = {'GT', 'DP', 'AD', 'GQ', 'PL'}
    absent = required_entry - set(mt.entry.dtype.fields)
    if absent or 'AQ' not in mt.info.dtype.fields:
        raise ValueError(f'Prepared MatrixTable is missing FORMAT {absent} or INFO/AQ')

    # Check the exact pedigree sample set before variant QC (original behaviour).
    vcf_samples = set(mt.aggregate_cols(hl.agg.collect(mt.s)))
    ped_samples = {x for trio in trios for x in (trio.s, trio.pat_id, trio.mat_id)}
    missing = sorted(ped_samples - vcf_samples)
    if missing:
        raise ValueError(f'{len(missing)} pedigree members absent from MatrixTable: {missing[:15]}')
    extras = vcf_samples - ped_samples
    if extras:
        mt = mt.filter_cols(hl.literal(ped_samples).contains(mt.s))
    print(f'Input samples: {len(vcf_samples):,}; pedigree members: {len(ped_samples):,}.', flush=True)

    # Recompute original *trio-cohort, pre-QC* AC from fully called GTs;
    # do not trust GLnexus full-cohort AQ/QUAL or old INFO/AC after downcoding.
    mt = hl.variant_qc(mt)
    mt = mt.annotate_rows(_preqc_ac=mt.variant_qc.AC[1])
    mt = mt.filter_rows(mt._preqc_ac < args.max_ac_exclusive)
    mt = mt.drop('variant_qc')

    if not args.skip_lcr:
        lcr = hl.import_bed(args.lcr_bed, reference_genome='GRCh38')
        mt = mt.filter_rows(~hl.is_defined(lcr[mt.locus]))

    # Exactly the old numeric QUAL cutoff; not equivalent calibration to
    # GATK QUAL, but the threshold itself is reproducible in this pVCF.
    mt = mt.filter_rows(hl.is_defined(mt.qual) & (mt.qual >= args.min_qual))

    # GLnexus FILTER='.' is NOT proven GATK PASS; omit the old SNP PASS check.
    # Exclude nonstandard MONOALLELIC regions (not ordinary diploid loci).
    mt = mt.filter_rows(hl.is_missing(mt.filters) |
                        ~mt.filters.contains('MONOALLELIC'))

    # Convert to trio columns: one entry contains child, father, mother.
    # Unlike the old code, all THREE calls must pass for the SAME trio.
    tm = hl.trio_matrix(mt, ped, complete_trios=True)
    kid, dad, mom = tm.proband_entry, tm.father_entry, tm.mother_entry
    kid_ab, dad_ab, mom_ab = valid_ab(kid), valid_ab(dad), valid_ab(mom)
    tm_snv = hl.is_snp(tm.alleles[0], tm.alleles[1])
    tm_indel = hl.is_indel(tm.alleles[0], tm.alleles[1])

    parent_snv = lambda entry, ab: (
        entry.GT.is_hom_ref() & (entry.GQ >= 30) & (ab <= 0.05)
    )
    parent_indel = lambda entry, ab: (
        entry.GT.is_hom_ref() & (entry.GQ >= 30) & (entry.DP >= 16) & (ab <= 0.05)
    )
    child_snv = kid.GT.is_het() & (kid.GQ >= 99) & (kid_ab >= 0.22) & (kid_ab <= 0.78)
    child_indel = kid.GT.is_het() & (kid.GQ >= 99) & (kid_ab >= 0.20) & (kid_ab <= 0.80)

    passes = (
        good_depth(kid) & good_depth(dad) & good_depth(mom)
        & (
            (tm_snv & parent_snv(dad, dad_ab) & parent_snv(mom, mom_ab) & child_snv)
            | (tm_indel & parent_indel(dad, dad_ab) & parent_indel(mom, mom_ab) & child_indel)
        )
    )
    tm = tm.filter_entries(hl.or_else(passes, False))

    # Entry table contains ONLY qualifying child-parent-parent trios.
    # The old script retained other Mendel-error types' entries by mistake;
    # this version keeps only members of genuinely qualifying code-2 trios.
    ht = tm.entries()
    candidate = ht.select(
        family_id=ht.fam_id,
        child_id=ht.id,
        father_id=ht.father.s,
        mother_id=ht.mother.s,
        child_GT=hl.str(ht.proband_entry.GT),
        father_GT=hl.str(ht.father_entry.GT),
        mother_GT=hl.str(ht.mother_entry.GT),
        child_DP=ht.proband_entry.DP,
        father_DP=ht.father_entry.DP,
        mother_DP=ht.mother_entry.DP,
        child_GQ=ht.proband_entry.GQ,
        father_GQ=ht.father_entry.GQ,
        mother_GQ=ht.mother_entry.GQ,
        child_AB=valid_ab(ht.proband_entry),
        father_AB=valid_ab(ht.father_entry),
        mother_AB=valid_ab(ht.mother_entry),
        # Record unfiltered depth components so later stringent parental ALT
        # support filters can run on merged TSVs without reopening wide VCFs.
        child_AD_ref=ht.proband_entry.AD[0],
        child_AD_alt=ht.proband_entry.AD[1],
        father_AD_ref=ht.father_entry.AD[0],
        father_AD_alt=ht.father_entry.AD[1],
        mother_AD_ref=ht.mother_entry.AD[0],
        mother_AD_alt=ht.mother_entry.AD[1],
        # Preserve site QC provenance from BEFORE trio-entry masking.
        PREQC_AC=ht._preqc_ac,
        QUAL=ht.qual,
        GLNEXUS_AQ=ht.info.AQ[ht.a_index - 1],
    )

    # Set of relevant samples per candidate allele. Candidate rarity makes
    # each per-site set small, avoiding a dense 16k x 7k trio cross-product.
    by_site = candidate.key_by('locus', 'alleles')
    by_site = by_site.select(
        child=by_site.child_id,
        father=by_site.father_id,
        mother=by_site.mother_id,
    ).collect_by_key()
    by_site = by_site.select(
        keep_samples=hl.set(hl.flatten(by_site.values.map(
            lambda x: hl.array([x.child, x.father, x.mother])
        ))),
        n_candidate_trios=hl.len(by_site.values),
    )

    mt = mt.semi_join_rows(by_site)
    mt = mt.annotate_rows(
        _eligible_samples=by_site[mt.row_key].keep_samples,
        _n_candidate_trios=by_site[mt.row_key].n_candidate_trios,
    )
    # Mask entries for everyone except members of a PASSING trio. All QC
    # criteria were applied jointly above, so the 3 qualifying calls survive.
    mt = mt.filter_entries(mt._eligible_samples.contains(mt.s))

    # Retain the original row-level mean GQ>=50 and final AC>0 checks.
    # (With >=99 child GQ and >=30 both parents, mean GQ>=50 is generally
    # redundant for qualifying trios, but maintained for transparency.)
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(
        (mt.variant_qc.gq_stats.mean >= 50) & (mt.variant_qc.AC[1] > 0)
    )

    if not args.no_checkpoints:
        # Small (post-Mendelian-screen) MatrixTable checkpoint avoids
        # recomputing the wide raw-VCF import for both TSV and VCF exports.
        final_mt_cp = _append_path(args.checkpoint_run_dir, 'final_sites.mt')
        print(f'Checkpointing candidate sites to {final_mt_cp}', flush=True)
        mt = mt.checkpoint(final_mt_cp)

    # TSV must reflect the final variant-level filters as well.
    final_sites = mt.rows().select().key_by('locus', 'alleles')
    candidate = candidate.filter(hl.is_defined(final_sites[candidate.locus, candidate.alleles]))
    candidate = candidate.annotate(
        contig=candidate.locus.contig,
        pos=candidate.locus.position,
        ref=candidate.alleles[0],
        alt=candidate.alleles[1],
    )
    candidate = candidate.key_by('contig', 'pos', 'ref', 'alt', 'child_id')
    candidate = candidate.select(
        'family_id', 'father_id', 'mother_id',
        'child_GT', 'father_GT', 'mother_GT',
        'child_DP', 'father_DP', 'mother_DP',
        'child_GQ', 'father_GQ', 'mother_GQ',
        'child_AB', 'father_AB', 'mother_AB',
        'child_AD_ref', 'child_AD_alt',
        'father_AD_ref', 'father_AD_alt',
        'mother_AD_ref', 'mother_AD_alt',
        'PREQC_AC', 'QUAL', 'GLNEXUS_AQ',
    )
    if not args.no_checkpoints:
        candidate_cp = _append_path(args.checkpoint_run_dir, 'candidate_trios.ht')
        print(f'Checkpointing candidate trio table to {candidate_cp}', flush=True)
        candidate = candidate.checkpoint(candidate_cp)

    # Counts refer to the final TSV/VCF after *all* current QC filters.
    # The candidate TSV, not the VCF, is the primary downstream artifact.
    candidate_events = candidate.count()
    candidate_sites = mt.count_rows()
    print(f'Candidate sites: {candidate_sites:,}; candidate trio-events: {candidate_events:,}', flush=True)
    candidate.export(args.stage_tsv)
    if not hfs.exists(args.stage_tsv):
        raise RuntimeError('Candidate table export did not create its output')

    if not args.no_vcf:
        preqc_ac = mt._preqc_ac
        qc = mt.variant_qc
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                AC=hl.array([qc.AC[1]]),
                AN=qc.AN,
                AF=hl.array([hl.or_missing(qc.AN > 0, hl.float64(qc.AC[1]) / qc.AN)]),
                AQ=hl.array([mt.info.AQ[mt.a_index - 1]]),
                PREQC_AC=preqc_ac,
                N_CANDIDATE_TRIOS=mt._n_candidate_trios,
            )
        )
        # Drop all non-VCF row fields and keep common entry fields; avoid
        # exporting nonstandard GLnexus RNC arrays, which describe the raw
        # unsplit call state and no longer match the split GT.
        mt = mt.select_rows('rsid', 'qual', 'filters', 'info')
        mt = mt.select_entries('GT', 'DP', 'AD', 'GQ', 'PL')
        metadata = {
            'info': {
                'AC': {'Number': 'A', 'Type': 'Integer', 'Description': 'ALT count among retained qualifying trio genotypes'},
                'AN': {'Number': '1', 'Type': 'Integer', 'Description': 'Called allele number among retained qualifying trio genotypes'},
                'AF': {'Number': 'A', 'Type': 'Float', 'Description': 'ALT frequency among retained qualifying trio genotypes'},
                'AQ': {'Number': 'A', 'Type': 'Integer', 'Description': 'GLnexus full-cohort allele quality, unchanged'},
                'PREQC_AC': {'Number': '1', 'Type': 'Integer', 'Description': 'Trio-cohort alternate allele count before genotype QC'},
                'N_CANDIDATE_TRIOS': {'Number': '1', 'Type': 'Integer', 'Description': 'Number of qualifying trios at this allele'},
            },
            'format': {
                'GT': {'Number': '1', 'Type': 'String', 'Description': 'Biallelic downcoded genotype'},
                'DP': {'Number': '1', 'Type': 'Integer', 'Description': 'Original genotype read depth'},
                'AD': {'Number': 'R', 'Type': 'Integer', 'Description': 'Biallelic downcoded allele depths'},
                'GQ': {'Number': '1', 'Type': 'Integer', 'Description': 'Original GLnexus genotype quality'},
                'PL': {'Number': 'G', 'Type': 'Integer', 'Description': 'Downcoded genotype likelihoods'},
            },
        }
        hl.export_vcf(mt, args.stage_vcf, metadata=metadata)
        if not hfs.exists(args.stage_vcf):
            raise RuntimeError('Candidate VCF export did not create its output')
    # Export to run-unique staging paths before publishing any final artifact.
    # The JSON marker is written LAST, and is the only success signal for
    # downstream Nextflow stages on eventually/non-atomically copied stores.
    artifacts = [(args.stage_tsv, args.final_tsv)]
    if not args.no_vcf:
        artifacts.append((args.stage_vcf, args.final_vcf))
    manifest = {
        'status': 'complete',
        'run_id': args.run_id,
        'chunk_id': args.chunk_id,
        'finished_utc': datetime.now(timezone.utc).isoformat(),
        'input_vcf': prep['input_vcf'],
        'input_version': prep['input_version'],
        'prepared_matrix': args.matrix,
        'preparation_manifest': args.prep_manifest,
        'pedigree': args.pedigree,
        'lcr_bed': None if args.skip_lcr else args.lcr_bed,
        'lcr_skipped': args.skip_lcr,
        'reference_genome': 'GRCh38',
        'hail_version': hl.__version__,
        'spark_java_heap_max_gib': heap_gb,
        'input_variant_records': input_records,
        'input_samples': len(vcf_samples),
        'pedigree_members': len(ped_samples),
        'pedigree_trios': len(trios),
        'candidate_sites': candidate_sites,
        'candidate_trio_events': candidate_events,
        'filters': {'min_qual_inclusive': args.min_qual,
                    'max_ac_exclusive': args.max_ac_exclusive},
        'outputs': [final for _, final in artifacts],
    }
    with hfs.open(args.stage_manifest, 'w') as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write('\n')
    _publish_results(artifacts, args.stage_manifest, args.final_manifest, args.overwrite,
                     obsolete_outputs=[args.final_vcf] if args.no_vcf else ())
    print(f'Published primary TSV: {args.final_tsv}', flush=True)
    if not args.no_vcf:
        print(f'Published validation VCF: {args.final_vcf}', flush=True)
    print(f'Completed: {args.final_manifest}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise
