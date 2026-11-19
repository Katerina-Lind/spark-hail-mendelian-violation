#!/usr/bin/env python3
"""Download per-chunk heterozygous SNV AB exports from S3 and plot sample-normalized histograms."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import boto3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prefix', required=True, help='s3://.../results/chr22_het_ab_v1')
    p.add_argument('--outdir', default='chr22_het_ab_analysis')
    p.add_argument('--expected-chunks',type=int,default=17)
    a = p.parse_args()
    if not a.prefix.startswith('s3://'):
        p.error('S3 prefix required')
    bucket, prefix = a.prefix[5:].split('/',1)
    prefix = prefix.rstrip('/') + '/'
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    s3 = boto3.client('s3',region_name='us-east-1')
    objects = []
    for page in s3.get_paginator('list_objects_v2').paginate(Bucket=bucket,Prefix=prefix):
        objects.extend(item['Key'] for item in page.get('Contents',[]) if item['Key'].endswith('.het_snp_ab.tsv.bgz'))
    objects.sort()
    if len(objects) != a.expected_chunks:
        raise RuntimeError(f'Expected {a.expected_chunks} AB files, found {len(objects)}. Check the Seqera run and S3 prefix')
    names = sorted(k.rsplit('/',1)[-1].replace('.het_snp_ab.tsv.bgz','') for k in objects)
    for chunk in names:
        key = prefix + chunk + '.het_snp_ab.completed.json'
        data = json.loads(s3.get_object(Bucket=bucket,Key=key)['Body'].read())
        if data.get('status') != 'complete' or data.get('chunk_id') != chunk:
            raise RuntimeError(f'Invalid completion marker for {chunk}')
    values = defaultdict(list)
    for i, key in enumerate(objects, 1):
        local = outdir / key.rsplit('/',1)[-1]
        if not local.exists():
            s3.download_file(bucket,key,str(local))
        with gzip.open(local,'rt') as f:
            reader = csv.DictReader(f,delimiter='\t')
            for row in reader:
                ab = float(row['AB'])
                if 0 <= ab <= 1:
                    values[row['sample_id']].append(ab)
                else:
                    raise RuntimeError(f'Unexpected AB {ab} at {key}')
        print(f'{i}/{len(objects)}: {local.name}',flush=True)
    if not values:
        raise RuntimeError('No AB observations found')
    with (outdir / 'ab_summary.tsv').open('w',newline='') as f:
        writer = csv.writer(f,delimiter='\t')
        writer.writerow(['sample_id','n_het_snvs','median_ab','mean_ab','pct_ab_020_030','pct_ab_040_060'])
        for sample in sorted(values):
            a_values = np.asarray(values[sample],dtype=np.float32)
            writer.writerow([sample,len(a_values),round(float(np.median(a_values)),5),
                round(float(a_values.mean()),5),
                round(100*float(np.mean((a_values >= .2)&(a_values < .3))),2),
                round(100*float(np.mean((a_values >= .4)&(a_values <= .6))),2)])
    fig, ax = plt.subplots(figsize=(10,5.8))
    bins = np.linspace(0,1,51)
    for sample in sorted(values):
        ax.hist(values[sample],bins=bins,weights=np.full(len(values[sample]),100/len(values[sample])),
                histtype='step',linewidth=1.8,label=f'{sample} (n={len(values[sample]):,})')
    ax.set(xlabel='Alternate allele balance (AD_alt / (AD_ref + AD_alt))',
           ylabel='Percent of heterozygous SNP calls per bin (width 0.02)',
           title='SPARK chr22: heterozygous SNP allele balance by sample')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / 'chr22_het_snp_ab_histogram.png',dpi=180)
    plt.close(fig)
    print('\nSummary:', outdir / 'ab_summary.tsv')
    print('Histogram:', outdir / 'chr22_het_snp_ab_histogram.png')

if __name__ == '__main__':
    main()
