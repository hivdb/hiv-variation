#! /usr/bin/env python
import csv
import json
import click

from drmlookup import build_algdrmlookup

GENE_CHOICES = ('PR', 'RT', 'IN')
NAIVE = True
TREATED = False
SIGNIFICANCE_LEVEL = 0.01
MAX_NAIVE_PCNT = 0.05
MIN_FOLD_CHANGE = 2
ALGDRMLOOKUP = build_algdrmlookup(1)

MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG']


@click.command()
@click.option('-i', '--prevalence-file', type=click.File('r'),
              help='input prevalence source')
@click.option('-o', '--output-file', type=click.File('w'),
              help='output target TSV')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
def create_prevalence_table(prevalence_file, output_file, gene):
    if gene in ('PR', 'RT'):
        raise NotImplementedError
    prevalence_data = json.load(prevalence_file)
    algdrmlookup = ALGDRMLOOKUP[gene]
    header = [
        'Position', 'AA',
        '# Naive Positive (All)',
        '% Naive Positive (All)',
        '# Naive Patients (All)',
        '# Treated Positive (All)',
        '% Treated Positive (All)',
        '# Treated Patients (All)'
    ]
    for subtype in MAJOR_SUBTYPES + ['Others']:
        header.extend([
            '% Naive Positive ({})'.format(subtype),
            '# Naive Patients ({})'.format(subtype),
        ])
    header.extend([
        '% Max Naive Positive',
        '# Max Naive Patients',
        '# Max Naive Subtype',
    ])
    writer = csv.DictWriter(
        output_file, header, extrasaction='ignore', delimiter='\t')
    writer.writeheader()
    rows = {}
    for item in prevalence_data:
        pos = item['position']
        aa = item['aa']
        if (pos, aa) not in algdrmlookup:
            continue
        if (pos, aa) not in rows:
            rows[(pos, aa)] = {
                'Position': pos,
                'AA': aa,
                '% Max Naive Positive': 0,
                '# Max Naive Patients': 0,
                '# Max Naive Subtype': '-',
            }
        row = rows[(pos, aa)]
        rx = item['rx_type']
        subtype = item['subtype']
        count = item['count']
        total = item['total']
        pcnt = item['percent']
        if subtype in ['All', 'Others'] + MAJOR_SUBTYPES:
            if rx == 'naive':
                row['# Naive Positive ({})'.format(subtype)] = count
                row['% Naive Positive ({})'.format(subtype)] = pcnt
                row['# Naive Patients ({})'.format(subtype)] = total
            if rx == 'art':
                row['# Treated Positive ({})'.format(subtype)] = count
                row['% Treated Positive ({})'.format(subtype)] = pcnt
                row['# Treated Patients ({})'.format(subtype)] = total
        if subtype not in ('All', 'Others') and rx == 'naive':
            if total < 100:
                # an arbitrary threshold
                continue
            if pcnt > row['% Max Naive Positive']:
                row['# Max Naive Positive'] = count
                row['% Max Naive Positive'] = pcnt
                row['# Max Naive Patients'] = total
                row['# Max Naive Subtype'] = subtype
    writer.writerows(rows.values())


if __name__ == '__main__':
    create_prevalence_table()
