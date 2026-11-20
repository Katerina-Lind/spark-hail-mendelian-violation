#!/usr/bin/env python3
"""Combine small per-chunk genotype exports, calculate KING-robust and IBS0.

KING between-family estimator matches the formula documented for hl.king().
This avoids constructing an enormous concatenated MatrixTable just to compare six people.
"""
import argparse
import csv
import gzip
import io
import itertools
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from urllib.parse import urlparse

import boto3

EXPECTED = {
    ('SDSM-V8W', 'SDSM-V8Y'): 'father-child',
    ('SDSM-V8W', 'SDSM-V8V'): 'mother-child',
    ('SDSM-2FNA', 'SDSM-2GIZ'): 'father-child',
    ('SDSM-2FNA', 'SDSM-2FNB'): 'mother-child',
}


def s3_parts(uri):
    p = urlparse(uri)
    if p.scheme != 's3' or not p.netloc or not p.path.strip('/'):
        raise ValueError('Expected S3 URI: ' + uri)
    return p.netloc, p.path.lstrip('/')


def get_bytes(client, uri):
    b, k = s3_parts(uri)
    return client.get_object(Bucket=b, Key=k)['Body'].read()


def put_new(client, uri, data, mime):
    bucket, key = s3_parts(uri)
    client.put_object(Bucket=bucket, Key=key, Body=data, ContentType=mime,
                      IfNoneMatch='*')


def main():
    p = argparse.ArgumentParser()
    for x in ('manifest', 'samples', 'results'):
        p.add_argument('--' + x, required=True)
    a = p.parse_args()
    samples = list(dict.fromkeys(s.strip() for s in a.samples.split(',') if s.strip()))
    if len(samples) != 6:
        raise ValueError('Expected six samples')
    s3 = boto3.client('s3', region_name='us-east-1')
    rows = list(csv.DictReader(io.StringIO(get_bytes(s3, a.manifest).decode()), delimiter='\t'))
    if len(rows) != 214:
        raise ValueError(f'Expected 214 chunk rows, got {len(rows)}')
    prefix = a.results.rstrip('/')
    # No old results should exist, including a partial previous execution.
    for suffix in ('kinship_pairs.tsv', 'kinship_summary.json', 'kinship.completed.json'):
        b, k = s3_parts(prefix + '/' + suffix)
        if s3.list_objects_v2(Bucket=b, Prefix=k, MaxKeys=1).get('KeyCount', 0):
            raise FileExistsError('Output exists: ' + k)

    by_site = {}
    sample_sites = Counter()
    chunk_loci = Counter()
    for row in rows:
        chunk = row['chunk_id']
        marker_uri = prefix + '/' + chunk + '.kinship_markers.completed.json'
        meta = json.loads(get_bytes(s3, marker_uri))
        if (meta.get('status') != 'complete' or meta.get('chunk_id') != chunk
                or meta.get('matrix') != row['matrix'] or meta.get('samples') != samples):
            raise ValueError(f'Chunk provenance mismatch: {chunk}')
        marker_file = prefix + '/' + chunk + '.kinship_markers.tsv.bgz'
        if meta.get('output') != marker_file:
            raise ValueError(f'Marker path mismatch: {chunk}')
        payload = get_bytes(s3, marker_file)
        with gzip.GzipFile(fileobj=io.BytesIO(payload)) as g:
            reader = csv.DictReader(io.TextIOWrapper(g, encoding='utf-8'), delimiter='\t')
            for r in reader:
                if r['sample_id'] not in samples:
                    raise ValueError('Unexpected sample in ' + chunk)
                if r['chromosome'] != row['chromosome']:
                    raise ValueError('Chromosome mismatch in ' + chunk)
                dosage = int(r['dosage'])
                if dosage not in (0, 1, 2):
                    raise ValueError('Non-diploid dosage')
                key = (r['chromosome'], int(r['pos']), r['ref'], r['alt'])
                if key not in by_site:
                    by_site[key] = {}
                    chunk_loci[chunk] += 1
                if r['sample_id'] in by_site[key]:
                    raise ValueError(f'Duplicate genotype for {key} and {r["sample_id"]}')
                by_site[key][r['sample_id']] = dosage
                sample_sites[r['sample_id']] += 1
    if not by_site:
        raise ValueError('Zero informative SNPs')
    missing = sorted(set(samples) - set(sample_sites))
    if missing:
        raise ValueError('Samples with zero called markers: ' + str(missing))

    all_chroms = sorted({r['chromosome'] for r in rows}, key=lambda x: int(x.replace('chr', '')))
    pairs = []
    for chrom in ['all'] + all_chroms:
        sites = [dosages for key, dosages in by_site.items()
                 if chrom == 'all' or key[0] == chrom]
        for a_id, b_id in itertools.combinations(sorted(samples), 2):
            both = [(s[a_id], s[b_id]) for s in sites if a_id in s and b_id in s]
            n = len(both)
            n_het_a = sum(x == 1 for x, _ in both)
            n_het_b = sum(y == 1 for _, y in both)
            het_both = sum(x == y == 1 for x, y in both)
            ibs0 = sum(abs(x-y) == 2 for x, y in both)
            den = n_het_a + n_het_b
            min_het = min(n_het_a, n_het_b)
            # Hail hl.king implements the between-family KING-robust estimator.
            phi_between = (0.5 + (2*het_both - 4*ibs0 - den) / (4*min_het)) if min_het else float('nan')
            phi_within = (het_both - 2*ibs0) / den if den else float('nan')
            role = EXPECTED.get((a_id, b_id), EXPECTED.get((b_id, a_id), 'other'))
            pairs.append({
                'chromosome': chrom, 'sample_a': a_id, 'sample_b': b_id,
                'recorded_relationship': role,
                'n_joint_called': n,
                'n_het_a': n_het_a, 'n_het_b': n_het_b,
                'n_het_both': het_both, 'n_ibs0': ibs0,
                'ibs0_fraction': ibs0 / n if n else float('nan'),
                'phi_king_between': phi_between, 'phi_king_within': phi_within,
            })
    fields = list(pairs[0])
    sink = io.StringIO()
    writer = csv.DictWriter(sink, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(pairs)
    summary = {
        'status': 'complete', 'method': 'KING between-family estimator (same equation as Hail hl.king)',
        'matrix_chunks': len(rows), 'chromosomes': all_chroms,
        'distinct_snp_sites': len(by_site),
        'samples': samples, 'called_snp_entries_per_sample': dict(sample_sites),
        'marker_count_per_chunk': dict(chunk_loci),
        'pair_count_per_partition': 15,
        'limitations': [
            'Only chr1, chr2, chr22 analysed, not genome-wide',
            'LD pruning occurs within each chunk, not across chunk boundaries',
            'Genotypes distorted by contamination can bias kinship and IBS0',
            'SNPs were position-thinned before full-cohort QC/LD pruning',
            'Six selected samples only: cannot identify unknown sample swaps elsewhere in cohort',
        ],
        'finished_utc': datetime.now(timezone.utc).isoformat(),
    }
    put_new(s3, prefix + '/kinship_pairs.tsv', sink.getvalue().encode(), 'text/tab-separated-values')
    put_new(s3, prefix + '/kinship_summary.json', json.dumps(summary, indent=2).encode(), 'application/json')
    put_new(s3, prefix + '/kinship.completed.json', json.dumps({
        'status': 'complete', 'pairs_file': prefix+'/kinship_pairs.tsv',
        'summary_file': prefix+'/kinship_summary.json',
        'distinct_snp_sites': len(by_site), 'finished_utc': summary['finished_utc']
    }, indent=2).encode(), 'application/json')
    print('Completed:', len(by_site), 'SNPs, 15 pairs per chromosome + combined', flush=True)
    print('Output:', prefix + '/kinship_pairs.tsv', flush=True)


if __name__ == '__main__':
    main()
