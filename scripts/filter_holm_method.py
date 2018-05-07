#! /usr/bin/env python
import os
import csv
import click


GENE_CHOICES = ('PR', 'RT', 'IN')
NAIVE_NONPOLYMORPHIC_THRESHOLD = 0.005
SIGNIFICANCE_LEVEL = 0.01
INTVALUES = ['Position', '# Naive Positive', '# Naive Patients',
             '# Treated Positive', '# Treated Patients']
FLOATVALUES = ['P Value']


@click.command()
@click.option('-i', '--indir', type=str, default='../data',
              help='input source directory')
@click.option('-o', '--outdir', type=str, default='../data',
              help='output target directory')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
def filter_holm_method(indir, outdir, gene):
    if gene in ('PR', 'RT'):
        raise NotImplementedError
    indir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), indir
    ))
    outdir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), outdir
    ))
    os.makedirs(outdir, exist_ok=True)
    file_chi2 = os.path.join(indir, 'in_chi2.txt')
    candidates = []
    with open(file_chi2) as file_chi2:
        reader = csv.DictReader(file_chi2, delimiter='\t')
        for row in reader:
            for intvalue in INTVALUES:
                row[intvalue] = int(row[intvalue])
            for floatvalue in FLOATVALUES:
                row[floatvalue] = float(row[floatvalue])
            naive_pcnt = row['# Naive Positive'] / row['# Naive Patients']
            if naive_pcnt >= NAIVE_NONPOLYMORPHIC_THRESHOLD:
                continue
            treated_pcnt = \
                row['# Treated Positive'] / row['# Treated Patients']
            row['% Naive Positive'] = '=TEXT({}, "0.00%")'.format(naive_pcnt)
            row['% Treated Positive'] = \
                '=TEXT({}, "0.00%")'.format(treated_pcnt)
            candidates.append(row)
    candidates = sorted(candidates, key=lambda c: c['P Value'])
    file_holm = os.path.join(outdir, 'in_holm_filtered.txt')
    with open(file_holm, 'w') as file_holm:
        writer = csv.DictWriter(
            file_holm,
            ['Position', 'AA', '% Naive Positive', '# Naive Positive',
             '# Naive Patients', '% Treated Positive', '# Treated Positive',
             '# Treated Patients', 'P Value'],
            delimiter='\t')
        writer.writeheader()
        for rank, cand in zip(range(len(candidates), 0, -1), candidates):
            pvalue = cand['P Value']
            if pvalue <= SIGNIFICANCE_LEVEL / rank:
                writer.writerow(cand)


if __name__ == '__main__':
    filter_holm_method()
