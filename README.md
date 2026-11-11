# Reusable GLnexus trio Mendelian code-2 screen (draft)

This draft separates a **durable per-chunk MatrixTable** from **analysis-specific filtering**. It preserves the student code's core autosomal Mendelian code-2 eligibility criteria and output column schema. This is **not** `hl.mendel_errors()` and **not** a validated de novo caller.

## Contents

- `scripts/01_prepare_matrix.py` — import GLnexus pVCF, restrict to autosomes, split multiallelic variants with original GQ, keep SNVs and indels, persist `.mt` and JSON manifest.
- `scripts/02_mendelian_screen.py` — read prepared `.mt`; apply pedigree restriction, AC, LCR, QUAL, MONOALLELIC and trio-specific genotype QC; emit candidate TSV, optional VCF and completion manifest.
- `scripts/student_original.py` — unchanged student's original code (line-number prefixes removed), for auditing and matched-run comparison.
- `main.nf` — initial DSL2 Seqera/AWS Batch orchestration example, not yet production-validated.

## Important safety constraints

1. Provide **unique versioned output paths** for each pVCF and preprocessing version. The script refuses mismatched or partially prepared datasets; it does **not verify** that a user-supplied `input_version` actually matches the current source object. Collect S3 VersionIds or checksums yourself and enforce immutability outside this script.
2. Two jobs must never prepare the same MatrixTable path concurrently. Pre-create a deduplicated chunk manifest and schedule one writer per MatrixTable. S3 directory writes are not atomic; interrupted writes require manual repair.
3. The preparation step assumes GLnexus `GT`, `DP`, `AD`, `GQ`, `PL` and `INFO/AQ`, GRCh38 `chr` contigs, and GT as the first FORMAT tag. Confirm the reduced SPARK pVCFs actually preserve these fields before using them.
4. The original numerical filtering rules are preserved. Original GQ is restored after `split_multi_hts`; **PL and AD remain downcoded**, so validate multiallelic cases separately.
5. Run the old and new implementations on the **same small pVCF chunk** and compare sorted TSV (including IDs, counts, AB, GQ, AQ, AC) and VCF before scaling. Beware shuffled row order and non-deterministic compression bytes.
6. `--no-vcf` is the suggested default for fast first-pass TSV analysis. `--lcr-bed` is required unless intentionally skipping the filter. Always record the original BED asset identity.
7. Keep preparation outputs indefinitely as long as they are needed; clean up per-run analysis checkpoints. The main `.mt` is not under the temporary checkpoint path.
8. These scripts have passed **Python syntax/static tests only**, not a Hail runtime integration test. Container/JVM/S3 permission configuration and performance must be tested on EC2/AWS Batch first.
9. Files are read directly from S3 by Hail. On AWS Batch, configure executor IAM permissions, container image, S3 workDir and scratch storage. Nextflow's `-resume` is complementary to, not a substitute for, persistent MatrixTables.
10. **There is a known scientific caveat:** the original code's final `variant_qc.gq_stats.mean` is evaluated after masking entries for qualifying trios. If this filter changes candidate outputs unexpectedly, compare its exact interpretation rather than silently removing it.

## Standalone preparation

```bash
python scripts/01_prepare_matrix.py \
  --vcf 's3a://BUCKET/input/chr22_chunk.vcf.gz' \
  --matrix 's3a://BUCKET/hail/v1/chr22_chunk.mt' \
  --manifest 's3a://BUCKET/hail/v1/chr22_chunk.prepare.json' \
  --input-version 'actual-S3-version-or-checksum' \
  --local-cores 8 --tmp-dir 's3a://BUCKET/hail/tmp'
```

Rerunning with exactly the same input identity and manifest reuses the MatrixTable without importing the VCF. If `--input-version` is inaccurate or the VCF was overwritten in place, this safety assumption fails; never treat a random label as an input checksum.

## Standalone analysis

```bash
python scripts/02_mendelian_screen.py \
  --matrix 's3a://BUCKET/hail/v1/chr22_chunk.mt' \
  --prep-manifest 's3a://BUCKET/hail/v1/chr22_chunk.prepare.json' \
  --pedigree 's3a://BUCKET/metadata/pedigrees.tsv' \
  --lcr-bed 's3a://BUCKET/reference/LCR-hs38.bed' \
  --out-prefix 's3a://BUCKET/results/chr22_chunk' \
  --local-cores 8 --tmp-dir 's3a://BUCKET/hail/tmp' --no-vcf
```

## Nextflow (after configuring AWS Batch + a validated container)

Populate `chunks.example.tsv` with one row for each unique nonoverlapping chunk; columns are `chunk_id`, `vcf`, `input_version`, `matrix`, `manifest`. S3 paths and container tag below are placeholders:

```bash
nextflow run main.nf -profile aws \
  --chunks chunks.tsv \
  --pedigree 's3://BUCKET/metadata/pedigrees.tsv' \
  --lcr_bed 's3://BUCKET/reference/LCR-hs38.bed' \
  --results 's3a://BUCKET/results/run_001' \
  --hail_tmp 's3a://BUCKET/hail/tmp' \
  --hail_container 'YOUR_TESTED_HAIL_IMAGE@sha256:...' \
  -work-dir 's3://BUCKET/nextflow-work'
```

The CLI flags `--lcr_bed`, `--hail_tmp` and `--hail_container` use Nextflow parameter names, not Python argument names. The `aws` profile does not include a guessed Batch queue; supply that via Seqera.

## Baseline filters carried across

- Selected samples: all pedigree trio members must exist; extraneous VCF samples removed before AC.
- GRCh38 autosomal SNPs/indels only, split using Hail HTS semantics, original GQ restored.
- Pre-genotype-QC trio-cohort AC strictly below 20 (default).
- LCR excluded unless `--skip-lcr`; `QUAL >= 150` (default); `MONOALLELIC` excluded.
- Both parents hom-ref and child heterozygous at same biallelic allele.
- All trio DP 10..200; parental GQ >=30; child GQ >=99; parental AB <=0.05.
- Child AB 0.22..0.78 SNVs, 0.20..0.80 indels; parental DP >=16 for indels.
- Final retained-trio genotype GQ mean >=50 and post-mask alternate-allele count >0.

## Still to benchmark

- On a small chunk, compare old versus new `candidate_trios.tsv.bgz` rows, and optionally the VCF genotype fields and INFO annotations.
- Confirm outputs can be exported to the intended S3 URI under the Hail version/container used in Seqera.
- Profile initial import, split, materialisation, candidate aggregation and export; quantify S3 storage costs before genome-wide expansion.
