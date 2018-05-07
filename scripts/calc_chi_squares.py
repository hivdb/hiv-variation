#! /usr/bin/env python
import os
import csv
import click
from collections import defaultdict

import numpy as np
from scipy.stats import chi2_contingency

GENE_CHOICES = ('PR', 'RT', 'IN')
NAIVE = True
TREATED = False


@click.command()
@click.option('-i', '--indir', type=str, default='../local',
              help='input source directory')
@click.option('-o', '--outdir', type=str, default='../data',
              help='output target directory')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
@click.argument('subtype', type=str, required=False)
def calc_chi_squares(indir, outdir, gene, subtype):
    if gene in ('PR', 'RT'):
        raise NotImplementedError
    indir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), indir
    ))
    outdir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), outdir
    ))
    os.makedirs(outdir, exist_ok=True)
    subtype_text = subtype or 'all'
    file_codons = os.path.join(indir, '{}_codons.txt'.format(gene.lower()))
    file_chi2 = os.path.join(
        outdir, '{}_{}_chi2.txt'.format(gene.lower(), subtype_text.lower()))
    counter = defaultdict(lambda: {
        NAIVE: set(),
        TREATED: set()
    })
    with open(file_codons) as file_codons:
        reader = csv.DictReader(file_codons, delimiter='\t')
        for codon in reader:
            if subtype and codon['Subtype'] != subtype:
                continue
            aa = codon['AA']
            if aa in '.X-' or len(aa) > 1:
                # skip unsequenced or mixtures
                continue
            pos = int(codon['Position'])
            is_naive = int(codon['NumDrugs']) == 0
            ptid = int(codon['PtID'])
            counter[(pos, aa)][is_naive].add(ptid)
    site_counts = {}
    total_counts = defaultdict(lambda: {
        NAIVE: 0, TREATED: 0
    })
    for (pos, aa), site_counter in counter.items():
        naive_count = len(site_counter[NAIVE])
        treated_count = len(site_counter[TREATED])

        site_counts[(pos, aa)] = {
            NAIVE: naive_count,
            TREATED: treated_count
        }
        total_counts[pos][NAIVE] += naive_count
        total_counts[pos][TREATED] += treated_count

    with open(file_chi2, 'w') as file_chi2:
        writer = csv.writer(file_chi2, delimiter='\t')
        writer.writerow([
            'Position', 'AA',
            '# Naive Positive', '# Naive Patients',
            '# Treated Positive', '# Treated Patients',
            'P Value'
        ])
        for (pos, aa), count in sorted(site_counts.items()):
            naive_total = total_counts[pos][NAIVE]
            naive_pos = count[NAIVE]
            naive_neg = naive_total - naive_pos
            treated_total = total_counts[pos][TREATED]
            treated_pos = count[TREATED]
            treated_neg = treated_total - treated_pos
            obs = np.array([
                [naive_pos, naive_neg],
                [treated_pos, treated_neg]
            ])
            if naive_neg == treated_neg == 0:
                p = 1.0
            else:
                _, p, _, _ = chi2_contingency(obs)
            writer.writerow([
                pos, aa,
                naive_pos, naive_total,
                treated_pos, treated_total, p])


if __name__ == '__main__':
    calc_chi_squares()
