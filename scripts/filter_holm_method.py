#! /usr/bin/env python
import os
import csv
import click


GENE_CHOICES = ('PR', 'RT', 'IN')
NAIVE_NONPOLYMORPHIC_THRESHOLD = 0.005
TREATED_NAIVE_RATIO_THRESHOLD = 5
SIGNIFICANCE_LEVEL = 0.01
INTVALUES = ['Position', '# Naive Positive', '# Naive Patients',
             '# Treated Positive', '# Treated Patients']
FLOATVALUES = ['% Naive Positive', '% Treated Positive',
               'P Value', 'Fold Change']


@click.command()
@click.option('-i', '--indir', type=str, default='../data',
              help='input source directory')
@click.option('-o', '--outdir', type=str, default='../data',
              help='output target directory')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
@click.argument('subtype', type=str, required=False)
def filter_holm_method(indir, outdir, gene, subtype):
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
    file_chi2 = os.path.join(
        indir, '{}_{}_chi2.txt'.format(gene.lower(), subtype_text.lower()))
    candidates = []
    with open(file_chi2) as file_chi2:
        reader = csv.DictReader(file_chi2, delimiter='\t')
        for row in reader:
            for intvalue in INTVALUES:
                row[intvalue] = int(row[intvalue])
            for floatvalue in FLOATVALUES:
                row[floatvalue] = float(row[floatvalue])
            naive_pcnt = row['% Naive Positive']
            if naive_pcnt >= NAIVE_NONPOLYMORPHIC_THRESHOLD:
                continue
            treated_pcnt = row['% Treated Positive']
            if naive_pcnt > 0 and \
                    treated_pcnt / naive_pcnt < TREATED_NAIVE_RATIO_THRESHOLD:
                continue
            candidates.append(row)
    candidates = sorted(candidates, key=lambda c: c['P Value'])
    file_holm = os.path.join(
        outdir, '{}_{}_tsm.txt'.format(gene.lower(), subtype_text.lower()))
    file_nontsm_drm = os.path.join(
        outdir, '{}_{}_nontsm_drm.txt'.format(gene.lower(),
                                              subtype_text.lower()))
    with open(file_holm, 'w') as file_holm, \
            open(file_nontsm_drm, 'w') as file_nontsm_drm:
        writer = csv.DictWriter(
            file_holm,
            ['Position', 'AA', '% Naive Positive', '# Naive Positive',
             '# Naive Patients', '% Treated Positive', '# Treated Positive',
             '# Treated Patients', 'P Value', 'Rank', 'P Value (Corrected)',
             'Fold Change', 'Is Major DRM', 'Is Accessory DRM', 'Is DRM',
             'Is Important Position'],
            delimiter='\t')
        drm_writer = csv.DictWriter(
            file_nontsm_drm,
            ['Position', 'AA', '% Naive Positive', '# Naive Positive',
             '# Naive Patients', '% Treated Positive', '# Treated Positive',
             '# Treated Patients', 'P Value', 'Rank', 'P Value (Corrected)',
             'Fold Change', 'Is Major DRM', 'Is Accessory DRM', 'Is DRM',
             'Is Important Position'],
            delimiter='\t')
        writer.writeheader()
        drm_writer.writeheader()
        significants = []
        nontsm_drms = []
        broken = False
        for rank, cand in zip(range(len(candidates), 0, -1), candidates):
            pvalue = cand['P Value']
            cand['Rank'] = rank
            cand['P Value (Corrected)'] = pvalue * rank
            if not broken and pvalue * rank <= SIGNIFICANCE_LEVEL:
                significants.append(cand)
            else:
                broken = True
                if cand['Is DRM'].upper() == 'TRUE':
                    nontsm_drms.append(cand)
        significants = sorted(significants,
                              key=lambda c: (c['Position'], c['AA']))
        writer.writerows(significants)
        nontsm_drms = sorted(nontsm_drms,
                             key=lambda c: (c['Position'], c['AA']))
        drm_writer.writerows(nontsm_drms)


if __name__ == '__main__':
    filter_holm_method()
