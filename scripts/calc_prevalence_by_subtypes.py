#! /usr/bin/env python
import os
import csv
import click
from collections import defaultdict

from drmlookup import build_algdrmlookup

GENE_CHOICES = ('PR', 'RT', 'IN')
NAIVE = True
TREATED = False
SIGNIFICANCE_LEVEL = 0.01
MAX_NAIVE_PCNT = 0.05
MIN_FOLD_CHANGE = 2
ALGDRMLOOKUP = build_algdrmlookup()

MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG']


@click.command()
@click.option('-i', '--indir', type=str, default='../local',
              help='input source directory')
@click.option('-o', '--outdir', type=str, default='../data',
              help='output target directory')
@click.argument('gene', required=True,
                type=click.Choice(GENE_CHOICES))
def calc_prevalence_by_subtypes(indir, outdir, gene):
    if gene in ('PR', 'RT'):
        raise NotImplementedError
    indir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), indir
    ))
    outdir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), outdir
    ))
    os.makedirs(outdir, exist_ok=True)
    file_codons = os.path.join(
        indir, '{}_codons.txt'.format(gene.lower()))
    file_output = os.path.join(
        outdir, '{}_prevalence.txt'.format(gene.lower()))

    counter = defaultdict(
        lambda: defaultdict(lambda: {
            NAIVE: set(),
            TREATED: set()
        })
    )

    with open(file_codons) as file_codons:
        reader = csv.DictReader(file_codons, delimiter='\t')
        for codon in reader:
            subtype = codon['Subtype']
            is_major_subtype = subtype in MAJOR_SUBTYPES
            aas = codon['AA']
            if aas in '.X-' or len(aas) > 2:
                # skip unsequenced or highly ambiguous mixtures
                continue
            pos = int(codon['Position'])
            is_naive = int(codon['NumDrugs']) == 0
            ptid = int(codon['PtID'])
            for aa in aas:
                counter[(pos, aa)]['.'][is_naive].add(ptid)
                counter[(pos, aa)][subtype][is_naive].add(ptid)
                if not is_major_subtype:
                    # one mutation by patient
                    counter[(pos, aa)]['Others'][is_naive].add(ptid)
    site_counts = defaultdict(lambda: defaultdict(
        lambda: {NAIVE: 0, TREATED: 0}
    ))
    total_counts = defaultdict(lambda: defaultdict(
        lambda: {NAIVE: 0, TREATED: 0}
    ))
    for (pos, aa), subtype_counter in counter.items():
        for subtype, site_counter in subtype_counter.items():
            naive_count = len(site_counter[NAIVE])
            treated_count = len(site_counter[TREATED])
            site_counts[(pos, aa)][subtype] = {
                NAIVE: naive_count,
                TREATED: treated_count
            }
            total_counts[pos][subtype][NAIVE] += naive_count
            total_counts[pos][subtype][TREATED] += treated_count

    algdrmlookup = ALGDRMLOOKUP[gene]
    with open(file_output, 'w') as file_output:
        writer = csv.writer(file_output, delimiter='\t')
        header = [
            'Position', 'AA',
            '# Naive Positive', '% Naive Positive', '# Naive Patients',
            '# Treated Positive', '% Treated Positive', '# Treated Patients'
        ]
        for subtype in MAJOR_SUBTYPES + ['Others']:
            header.extend([
                '# Naive Positive ({})'.format(subtype),
                '% Naive Positive ({})'.format(subtype),
                '# Naive Patients ({})'.format(subtype),
                '# Treated Positive ({})'.format(subtype),
                '% Treated Positive ({})'.format(subtype),
                '# Treated Patients ({})'.format(subtype),
            ])
        header.extend([
            '# Max Naive Positive',
            '% Max Naive Positive',
            '# Max Naive Patients',
            '# Max Naive Subtype',
        ])
        writer.writerow(header)
        for (pos, aa), subtype_counts in sorted(site_counts.items()):
            if (pos, aa) not in algdrmlookup:
                continue
            row = [pos, aa]
            for subtype in ['.'] + MAJOR_SUBTYPES + ['Others']:
                count = subtype_counts[subtype]
                naive_total = total_counts[pos][subtype][NAIVE]
                naive_pos = count[NAIVE]
                treated_total = total_counts[pos][subtype][TREATED]
                treated_pos = count[TREATED]
                naive_pos_pcnt = naive_pos / naive_total
                treated_pos_pcnt = treated_pos / treated_total

                row.extend([
                    naive_pos, naive_pos_pcnt, naive_total,
                    treated_pos, treated_pos_pcnt, treated_total,
                ])
            max_naive = [0, 0, 0, '-']
            for subtype, count in subtype_counts.items():
                if subtype == 'Others':
                    continue
                naive_total = total_counts[pos][subtype][NAIVE]
                naive_pos = count[NAIVE]
                if naive_total < 100:
                    # an arbitrary threshold
                    continue
                naive_pos_pcnt = naive_pos / naive_total
                if naive_pos_pcnt > max_naive[1]:
                    max_naive = [naive_pos, naive_pos_pcnt,
                                 naive_total, subtype]
            row.extend(max_naive)
            writer.writerow(row)


if __name__ == '__main__':
    calc_prevalence_by_subtypes()
