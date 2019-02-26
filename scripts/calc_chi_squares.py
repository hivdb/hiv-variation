#! /usr/bin/env python
import csv
import json
import click

import numpy as np
from scipy.stats import chi2_contingency

from drmlookup import build_drmlookup, build_algdrmlookup_with_numalgs

GENE_CHOICES = ('PR', 'RT', 'IN')
NAIVE = True
TREATED = False
SIGNIFICANCE_LEVEL = 0.01
MAX_NAIVE_PCNT = 0.05
MIN_FOLD_CHANGE = 2
DRMLOOKUP = build_drmlookup()
ALGDRMLOOKUP = build_algdrmlookup_with_numalgs()

IMPORTANT_POSITIONS = {
    'PR': [],
    'RT': [],
    'IN': {66, 92, 118, 138, 140, 143, 147, 148, 155, 263}
}


@click.command()
@click.option('-i', '--prevalence-file', type=click.File('r'),
              help='input prevalence source')
@click.option('-o', '--output-file', type=click.File('w'),
              help='output target TSV')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
def calc_chi_squares(prevalence_file, output_file, gene):
    if gene in ('PR', 'RT'):
        raise NotImplementedError
    prevalence_data = json.load(prevalence_file)

    drmlookup = DRMLOOKUP[gene]
    algdrmlookup = ALGDRMLOOKUP[gene]
    impposlookup = IMPORTANT_POSITIONS[gene]
    writer = csv.DictWriter(
        output_file,
        [
            'Position', 'AA',
            '# Naive Positive', '% Naive Positive', '# Naive Patients',
            '# Treated Positive', '% Treated Positive', '# Treated Patients',
            'P Value', 'Fold Change', 'Is Major DRM', 'Is Accessory DRM',
            'Is DRM', '# Algs', 'Is Important Position', 'Is Significant'
        ], delimiter='\t')
    writer.writeheader()
    rows = {}
    for item in prevalence_data:
        if item['subtype'] != 'All':
            # we only want all subtypes
            continue
        pos = item['position']
        aa = item['aa']
        if (pos, aa) not in rows:
            is_major = (pos, aa, True) in drmlookup
            is_acc = (pos, aa, False) in drmlookup
            rows[(pos, aa)] = {
                'Position': pos,
                'AA': aa,
                'Is Major DRM': is_major,
                'Is Accessory DRM': is_acc,
                'Is DRM': is_major or is_acc,
                '# Algs': algdrmlookup.get((pos, aa), 0),
                'Is Important Position': pos in impposlookup,
            }
        row = rows[(pos, aa)]
        rx = item['rx_type']
        count = item['count']
        total = item['total']
        pcnt = item['percent']
        if rx == 'naive':
            row['# Naive Positive'] = count
            row['% Naive Positive'] = pcnt
            row['# Naive Patients'] = total
        else:
            row['# Treated Positive'] = count
            row['% Treated Positive'] = pcnt
            row['# Treated Patients'] = total
    for row in rows.values():
        naive_pos = row['# Naive Positive']
        naive_neg = row['# Naive Patients'] - naive_pos
        treated_pos = row['# Treated Positive']
        treated_neg = row['# Treated Patients'] - treated_pos
        obs = np.array([
            [naive_pos, naive_neg],
            [treated_pos, treated_neg]
        ])
        if naive_pos == treated_pos == 0:
            continue
        elif naive_neg == treated_neg == 0:
            p = 1.0
        else:
            _, p, _, _ = chi2_contingency(obs)
        fold_change = 1e4
        naive_pos_pcnt = row['% Naive Positive']
        treated_pos_pcnt = row['% Treated Positive']
        if row['% Naive Positive'] > 0:
            fold_change = (treated_pos_pcnt / naive_pos_pcnt)
        is_sig = row['# Algs'] > 1 or (p < SIGNIFICANCE_LEVEL and
                                       naive_pos_pcnt < MAX_NAIVE_PCNT and
                                       fold_change > MIN_FOLD_CHANGE and
                                       naive_pos_pcnt < treated_pos_pcnt)
        row['P Value'] = p
        row['Fold Change'] = fold_change
        row['Is Significant'] = is_sig
        writer.writerow(row)


if __name__ == '__main__':
    calc_chi_squares()
