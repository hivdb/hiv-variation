#! /usr/bin/env python
import csv
import json
import click

import numpy as np
from scipy.stats import chi2_contingency


@click.command()
@click.option('-i', '--prevalence-file', type=click.File('r'),
              help='input prevalence source')
@click.option('-o', '--output-file', type=click.File('w'),
              help='output target TSV')
@click.option('--subtypes', multiple=True, type=str,
              default=('A', 'B'), show_default=True,
              help='calc for these subtypes')
def calc_chi_squares(prevalence_file, output_file, subtypes):
    prevalence_data = json.load(prevalence_file)

    writer = csv.DictWriter(
        output_file, ['Gene', 'Position', 'AA', 'P Value'],
        delimiter='\t')
    writer.writeheader()
    rows = {}
    for item in prevalence_data:
        gene = item['gene']
        pos = item['position']
        aa = item['aa']
        rx = item['rx_type']
        subtype = item['subtype']
        if subtype not in subtypes or rx not in ('art', 'naive'):
            continue
        if (gene, pos, aa) not in rows:
            rows[(gene, pos, aa)] = {}
        row = rows[(gene, pos, aa)]
        row[(rx, subtype, 0)] = item['count']
        row[(rx, subtype, 1)] = item['total'] - item['count']
    for (gene, pos, aa), row in rows.items():
        obs0 = []
        postotal = 0
        negtotal = 0
        for subtype in subtypes:
            obs1 = []
            for rx in ('art', 'naive'):
                pos = row.get((rx, subtype, 0), 0)
                neg = row.get((rx, subtype, 1), 0)
                obs1.append([pos, neg])
                postotal += pos
                negtotal += neg
            obs0.append(obs1)
        obs = np.array(obs0)
        if postotal == 0:
            continue
        elif negtotal == 0:
            p = 1.0
        else:
            _, p, _, _ = chi2_contingency(obs)
        writer.writerow({
            'Gene': gene,
            'Position': pos,
            'AA': aa,
            'P Value': p
        })


if __name__ == '__main__':
    calc_chi_squares()
