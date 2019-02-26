#! /usr/bin/env python
import re
import csv
import click
import requests

from collections import Counter
from drmlookup import build_algdrmlookup_with_numalgs

GENOPHENO_DATASET_URL = {
    'PI': ('https://hivdb.stanford.edu/'
           'download/GenoPhenoDatasets/PI_DataSet.Full.txt'),
    'NRTI': ('https://hivdb.stanford.edu/'
             'download/GenoPhenoDatasets/NRTI_DataSet.Full.txt'),
    'NNRTI': ('https://hivdb.stanford.edu/'
              'download/GenoPhenoDatasets/NNRTI_DataSet.Full.txt'),
    'INSTI': ('https://hivdb.stanford.edu/'
              'download/GenoPhenoDatasets/INI_DataSet.Full.txt')
}

ALGDRMLOOKUP = build_algdrmlookup_with_numalgs()
DRUG_CLASS_GENE_MAP = {
    'PI': 'PR',
    'NRTI': 'RT',
    'NNRTI': 'RT',
    'INSTI': 'IN'
}


def load_phenodata(drug_class, filter_func=lambda row: True):
    if drug_class not in GENOPHENO_DATASET_URL:
        raise NotImplementedError
    resp = requests.get(GENOPHENO_DATASET_URL[drug_class])
    reader = csv.DictReader(
        resp.iter_lines(decode_unicode=True), delimiter='\t')
    for row in reader:
        if filter_func(row):
            yield row


def create_filter_func(method, all_methods):

    if all_methods:
        return lambda row: True

    method = set(method)

    def filter_func(row):
        if row['Method'] in method:
            return True
        return False

    return filter_func


def iter_mutations(phenodata):
    """Iterate every mutation from giving genopheno dataset"""
    key_pattern = re.compile(r'^P\d+$')
    for row in phenodata:
        for key, aas in row.items():
            if aas in ('-', '.') or not key_pattern.match(key):
                continue
            aas = aas.replace('#', 'i').replace('~', 'd')
            pos = int(key[1:])
            for aa in aas:
                yield pos, aa


@click.command()
@click.option('-o', '--output-file', type=click.File('w'), default='-',
              show_default=True, help='output target TSV')
@click.option('--method', type=str, multiple=True,
              default=('PhenoSense', 'HeLaCD4-ViiV'),
              show_default=True, help='include methods only')
@click.option('--all-methods', is_flag=True,
              help='include all methods')
@click.argument('drug_class', required=True,
                type=click.Choice(GENOPHENO_DATASET_URL.keys()))
def stat_phenomuts(output_file, method, all_methods, drug_class):
    filter_func = create_filter_func(method, all_methods)
    phenodata = load_phenodata(drug_class, filter_func)
    mutations = iter_mutations(phenodata)
    result = Counter(mutations)
    result = sorted(result.items(), key=lambda r: (r[0][0], -r[1]))
    writer = csv.writer(output_file, delimiter='\t')
    writer.writerow(['# Filter conditions:'])
    writer.writerow(['# - Method IN ({})'.format(', '.join(method))])
    writer.writerow([])
    writer.writerow(['Position', 'AA', '# Algs', 'Count'])
    gene = DRUG_CLASS_GENE_MAP[drug_class]
    algdrmlookup = ALGDRMLOOKUP[gene]
    for (pos, aa), count in result:
        numalgs = algdrmlookup.get((pos, aa), 0)
        writer.writerow([pos, aa, numalgs, count])


if __name__ == '__main__':
    stat_phenomuts()
