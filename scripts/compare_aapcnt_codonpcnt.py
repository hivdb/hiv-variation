#! /usr/bin/env python

import csv
import click
from hivfacts import HIVAAPcnt, HIVCodonPcnt

from codonutils import translate_codon, get_codons


def codon2aa(codon):
    if codon == 'del':
        return '-'
    elif codon == 'ins':
        return '_'
    else:
        return translate_codon(codon)


@click.command()
@click.option('--subtype', type=str, default='all',
              help='specify an HIV subtype/group')
@click.option('--rx-type', type=click.Choice(['art', 'naive', 'all']),
              default='all', help='specify treatment type')
@click.argument('output', type=click.File('w'))
def main(subtype, rx_type, output):
    aapcnt_reader = HIVAAPcnt(rx_type, subtype)
    codonpcnt_reader = HIVCodonPcnt(rx_type, subtype)
    writer = csv.writer(output)
    writer.writerow(['Gene', 'Position', 'AA', 'Codon',
                     'NumPossibleCodons',
                     'AAPcnt', 'CodonPcnt'])

    num_possible_codons = {
        '-': 1,
        '_': 1
    }

    for codon_obj in codonpcnt_reader.get():
        gene = codon_obj['gene']
        position = codon_obj['position']
        codon = codon_obj['codon']
        aa = codon2aa(codon)
        aa_obj = aapcnt_reader.get(gene, position, aa)
        aa_pcnt = aa_obj['percent']
        cd_pcnt = codon_obj['percent']
        if aa not in num_possible_codons:
            num_possible_codons[aa] = len(get_codons(aa))
        if cd_pcnt < 0.0001:
            writer.writerow([
                gene, position, aa, codon,
                num_possible_codons[aa],
                '{:.10f}%'.format(100 * aa_pcnt),
                '{:.10f}%'.format(100 * cd_pcnt)
            ])


if __name__ == '__main__':
    main()
