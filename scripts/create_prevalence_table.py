#! /usr/bin/env python
import csv
import json
import click

from drmlookup import build_algdrmlookup_with_numalgs

import numpy as np
from scipy.stats import chi2_contingency


GENE_CHOICES = ('PR', 'RT', 'IN')
SIGNIFICANCE_LEVEL = 0.01
MIN_TREATED_CASES = 3
MAX_NAIVE_PCNT = 0.005
MIN_FOLD_CHANGE = 2
ALGDRMLOOKUP = build_algdrmlookup_with_numalgs()

MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG', 'D', 'F', 'G']


@click.command()
@click.option('-i', '--prevalence-file', type=click.File('r'),
              help='input prevalence source')
@click.option('-o', '--output-file', type=click.File('w'),
              help='output target TSV')
@click.option('--num-algs-range', type=int, nargs=2,
              default=(1, 4), show_default=True,
              help='specify the # Algs range of including mutations')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
def create_prevalence_table(
        prevalence_file, output_file, num_algs_range, gene):
    if gene in ('PR', 'RT'):
        raise NotImplementedError
    prevalence_data = json.load(prevalence_file)
    algdrmlookup = ALGDRMLOOKUP[gene]
    header = [
        'Position', 'AA',
        '# Naive (All)',
        '# Naive Cases (All)',
        'Naive Prev (All)',
        '# Treated (All)',
        '# Treated Cases (All)',
        'Treated Prev (All)'
    ]
    for subtype in MAJOR_SUBTYPES + ['Others']:
        header.extend([
            'Naive Prev ({})'.format(subtype),
            '# Naive ({})'.format(subtype),
        ])
    header.extend([
        'Max Naive Total',
        'Max Naive Cases',
        'Max Naive Prev',
        'Max Naive Subtype',
        '# Algs',
        'P Value',
        'Fold Change',
        'Is Significant',
    ])
    writer = csv.DictWriter(
        output_file, header, extrasaction='ignore', delimiter='\t')
    writer.writeheader()
    rows = {}
    for item in prevalence_data:
        pos = item['position']
        aa = item['aa']
        num_algs = algdrmlookup.get((pos, aa), 0)
        if num_algs < num_algs_range[0] or num_algs > num_algs_range[1]:
            continue
        if (pos, aa) not in rows:
            rows[(pos, aa)] = {
                'Position': pos,
                'AA': aa,
                'Max Naive Prev': 0,
                'Max Naive Total': 0,
                'Max Naive Subtype': '-',
                '# Algs': num_algs,
            }
        row = rows[(pos, aa)]
        rx = item['rx_type']
        subtype = item['subtype']
        count = item['count']
        total = item['total']
        pcnt = item['percent']
        if subtype in ['All', 'Others'] + MAJOR_SUBTYPES:
            if rx == 'naive':
                row['# Naive Cases ({})'.format(subtype)] = count
                row['Naive Prev ({})'.format(subtype)] = pcnt
                row['# Naive ({})'.format(subtype)] = total
            if rx == 'art':
                row['# Treated Cases ({})'.format(subtype)] = count
                row['Treated Prev ({})'.format(subtype)] = pcnt
                row['# Treated ({})'.format(subtype)] = total
        if subtype not in ('All', 'Others', 'Unknown') and rx == 'naive':
            if total < 200:
                # an arbitrary threshold
                continue
            if pcnt > row['Max Naive Prev']:
                row['Max Naive Cases'] = count
                row['Max Naive Prev'] = pcnt
                row['Max Naive Total'] = total
                row['Max Naive Subtype'] = subtype
    for row in rows.values():
        naive_pos = row['# Naive Cases (All)']
        naive_neg = row['# Naive (All)'] - naive_pos
        treated_pos = row['# Treated Cases (All)']
        treated_neg = row['# Treated (All)'] - treated_pos
        obs = np.array([
            [naive_pos, naive_neg],
            [treated_pos, treated_neg]
        ])
        if naive_pos == treated_pos == 0:
            p = 1.0
        elif naive_neg == treated_neg == 0:
            p = 1.0
        else:
            _, p, _, _ = chi2_contingency(obs)
        fold_change = 1e4
        naive_pos_pcnt = row['Naive Prev (All)']
        treated_pos_pcnt = row['Treated Prev (All)']
        if naive_pos_pcnt > 0:
            fold_change = (treated_pos_pcnt / naive_pos_pcnt)
        is_sig = row['# Algs'] > 1 or (treated_pos >= MIN_TREATED_CASES and
                                       naive_pos_pcnt < MAX_NAIVE_PCNT and
                                       fold_change >= MIN_FOLD_CHANGE)
        row['P Value'] = p
        row['Fold Change'] = fold_change
        row['Is Significant'] = is_sig
        writer.writerow(row)


if __name__ == '__main__':
    create_prevalence_table()
