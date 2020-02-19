#! /usr/bin/env python
import csv
import json
import click
from itertools import groupby

from hivfacts import hivsdrm

import numpy as np
from scipy.stats import fisher_exact


GENE_CHOICES = ('PR', 'RT', 'IN')
SIGNIFICANCE_LEVEL = 0.01
MIN_TREATED_CASES = 3
MAX_NAIVE_PCNT = 0.005
MIN_FOLD_CHANGE = 2
SDRM_LOOKUP = hivsdrm.HIVSDRM()


def aapcnt_sortkey(item):
    return item['gene'], item['drugclass'], item['position'], item['aa']


@click.command()
@click.argument('drm_feature_csv', type=click.File('r'))
@click.argument('extended_csv', type=click.File('w'))
@click.option('--aapcnt', type=click.File('r'),
              help='input prevalence source')
@click.option('--major-subtypes', multiple=True, type=str,
              default=('A', 'B', 'C', 'CRF01_AE', 'CRF02_AG', 'D', 'F', 'G'),
              show_default=True, help='stat for these subtypes')
def extend_drm_feature(drm_feature_csv, extended_csv, aapcnt, major_subtypes):
    rows = list(csv.reader(drm_feature_csv))
    orig_header = rows.pop(0)
    prevalence_data = json.load(aapcnt)
    prevalence_data = sorted(prevalence_data, key=aapcnt_sortkey)
    prevalence_data = {
        mut: list(items)
        for mut, items in groupby(prevalence_data, aapcnt_sortkey)
    }
    header = [
        *orig_header,
        'SDRM',
        '# Naive (All)',
        '# Naive Cases (All)',
        '# Naive w/ SDRM',
        '# Naive w/o SDRM',
        'Naive Prev (All)',
        '# Treated (All)',
        '# Treated Cases (All)',
        '# Treated w/ SDRM',
        '# Treated w/o SDRM',
        'Treated Prev (All)'
    ]
    for subtype in list(major_subtypes) + ['Others']:
        header.extend([
            'Naive Prev ({})'.format(subtype),
            '# Naive ({})'.format(subtype),
        ])
    header.extend([
        'Max Naive Total',
        'Max Naive Cases',
        'Max Naive Prev',
        'Max Naive Subtype',
    ])
    header.extend([
        'P Value',
        'Fold Change',
    ])
    writer = csv.DictWriter(
        extended_csv, header, extrasaction='ignore')
    writer.writeheader()
    for row in rows:
        gene, pos, aa, drugclass = tuple(row[:4])
        pos = int(pos)
        mutkey = (gene, drugclass, pos, aa)
        # normalize indels
        row[2] = row[2].replace('-', 'del').replace('_', 'ins')
        row = dict(zip(orig_header, row))
        row.update({
            'Max Naive Prev': '0%',
            'Max Naive Total': 0,
            'Max Naive Subtype': '-',
        })
        row['SDRM'] = 1 if SDRM_LOOKUP.is_sdrm(gene, pos, aa) else ''
        items = prevalence_data[mutkey]
        for item in items:
            rx = item['rx_type']
            subtype = item['subtype']
            count = item['count']
            total = item['total']
            pcnt = item['percent']
            if rx == 'naive' and subtype == 'All':
                row['# Naive w/o SDRM'] = item['without_sdrm']
                row['# Naive w/ SDRM'] = item['with_sdrm']
            if rx == 'art' and subtype == 'All':
                row['# Treated w/o SDRM'] = item['without_sdrm']
                row['# Treated w/ SDRM'] = item['with_sdrm']
            if subtype in ['All', 'Others'] + list(major_subtypes):
                if rx == 'naive':
                    row['# Naive Cases ({})'.format(subtype)] = count
                    row['Naive Prev ({})'.format(subtype)] = \
                        '{}%'.format(pcnt * 100)
                    row['# Naive ({})'.format(subtype)] = total
                if rx == 'art':
                    row['# Treated Cases ({})'.format(subtype)] = count
                    row['Treated Prev ({})'.format(subtype)] = \
                        '{}%'.format(pcnt * 100)
                    row['# Treated ({})'.format(subtype)] = total
            if subtype not in ('All', 'Others', 'Unknown') and rx == 'naive':
                if total < 200:
                    # an arbitrary threshold
                    continue
                if pcnt > float(row['Max Naive Prev'][:-1]) / 100:
                    row['Max Naive Cases'] = count
                    row['Max Naive Prev'] = '{}%'.format(pcnt * 100)
                    row['Max Naive Total'] = total
                    row['Max Naive Subtype'] = subtype
        naive_pos = row['# Naive Cases (All)']
        naive_neg = row['# Naive (All)'] - naive_pos
        treated_pos = row['# Treated Cases (All)']
        treated_neg = row['# Treated (All)'] - treated_pos
        obs = np.array([
            [naive_pos, naive_neg],
            [treated_pos, treated_neg]
        ])
        try:
            _, p = fisher_exact(obs)
            # _, p, _, _ = chi2_contingency(obs)
        except ValueError:
            p = 1.0
        fold_change = 1e2
        naive_pos_pcnt = float(row['Naive Prev (All)'][:-1]) / 100
        treated_pos_pcnt = float(row['Treated Prev (All)'][:-1]) / 100
        if naive_pos_pcnt > 0:
            fold_change = (treated_pos_pcnt / naive_pos_pcnt)
        row['P Value'] = p
        row['Fold Change'] = fold_change
        writer.writerow(row)


if __name__ == '__main__':
    extend_drm_feature()
