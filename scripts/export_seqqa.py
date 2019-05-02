#! /usr/bin/env python
import csv
import json
from collections import defaultdict

import click
from hivdbql import app

db = app.db
models = app.models

Isolate = models.Isolate
Host = models.Host
Species = models.Species
ClinicalIsolate = models.ClinicalIsolate
Subtype = models.Subtype
Sequence = models.Sequence
Reference = models.Reference

UNUSUAL_AAPCNT_THRESHOLD = {
    'PR': 0.002,  # 0.2%
    'RT': 0.002,  # 0.2%
    'IN': 0.005   # 0.5%
}
GENES = ('PR', 'RT', 'IN')
DRUG_CLASSES = ('PI', 'NRTI', 'NNRTI', 'RTI', 'INSTI')
DRUG_CLASS_GENE_MAP = {
    'PI': 'PR',
    'NRTI': 'RT',
    'NNRTI': 'RT',
    'RTI': 'RT',
    'INSTI': 'IN'
}
MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG', 'D', 'F', 'G']
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY_-*'
# UNUSUAL_CUTOFF = 0.0001  # 1 in 10,000 or 0.01%
CSV_HEADER = [
    'IsolateID',
    'Gene',
    '# Unusuals',
    'Unusuals',
    '# APOBECs',
    'APOBECs'
]
QUERY_CHUNK_SIZE = 500

CRITERIA_CHOICES = {
    'HIV1_ONLY': Isolate._species.has(Species.species == 'HIV1'),
    'HIV2_ONLY': Isolate._species.has(Species.species == 'HIV2'),
    'PLASMA_ONLY': Isolate.clinical_isolate.has(
        ClinicalIsolate.source == 'Plasma'),
    'NO_CLONES': Isolate.clinical_isolate.has(
        ClinicalIsolate.clone_method == 'None'),
    'NO_QA_ISSUES': ~Isolate._filter.has(),
    'GENBANK_ONLY': Isolate.sequences.any(
        Sequence.accession.isnot(None) &
        (Sequence.accession != '')
    ),
    'NO_PARTIAL_MUTS': Isolate.sequences.any(
        Sequence.sequence_type == 'PartialMutationList'
    ),
}


def build_consensus_lookup(aapcnt_data):
    table = defaultdict(lambda: (None, -1.))
    for aapcnt in aapcnt_data:
        if aapcnt['subtype'] != 'All' or aapcnt['rx_type'] != 'all':
            continue
        gene = aapcnt['gene']
        pos = aapcnt['position']
        table[(gene, pos)] = max(
            table[(gene, pos)],
            (aapcnt['aa'], aapcnt['percent']),
            key=lambda o: o[1])
    return table


def unusual_mutation_lookup(aapcnt_data):
    table = {}
    for aapcnt in aapcnt_data:
        if aapcnt['subtype'] != 'All' or aapcnt['rx_type'] != 'all':
            continue
        gene = aapcnt['gene']
        pcnt = aapcnt['percent']
        aa = aapcnt['aa']
        if aa != '*' and pcnt > UNUSUAL_AAPCNT_THRESHOLD[gene]:
            continue
        # TODO: HIV2 only
        if gene == 'RT' and pcnt > 240:
            continue
        if gene == 'IN' and pcnt > 270:
            continue
        table[(gene, aapcnt['position'], aa)] = pcnt
    return table


def apobec_mutation_lookup(apobec_json):
    apobec_data = json.load(apobec_json)
    table = set()
    for apobec in apobec_data:
        table.add((apobec['gene'], apobec['position'], apobec['aa']))
    return table


def iter_isolates(drugclass, criteria, is_hiv2):
    print('Processing {} isolates...'
          .format(drugclass))
    gene = DRUG_CLASS_GENE_MAP[drugclass]
    if is_hiv2:
        criteria += ('HIV2_ONLY',)
    else:
        criteria += ('HIV1_ONLY',)
    conds = [CRITERIA_CHOICES[crkey] for crkey in criteria]
    query = (
        Isolate.query
        .filter(
            Isolate.gene == gene,
            Isolate.isolate_type == 'Clinical',
            Isolate._host.has(Host.host == 'Human'),
            *conds
        )
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.insertions))
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.mixtures))
    )
    if not is_hiv2:
        # for old HIV-2 isolate, there's no subtype table record
        query = query.filter(
            Isolate._subtype.has(Subtype.subtype.notin_(
                ['O', 'N', 'P', 'CPZ']
            ))
        )

    total = query.count()
    query = query.order_by(Isolate.id)
    for offset in range(0, total, QUERY_CHUNK_SIZE):
        print('  {}/{} isolates...'.format(offset, total), end='\r')
        yield from query.limit(QUERY_CHUNK_SIZE).offset(offset)
    print('  {0} isolates...                   '.format(total))


def run_seqqa(drugclass, criteria, is_hiv2,
              cons_lookup, uum_lookup, apm_lookup):
    for isolate in iter_isolates(drugclass, criteria, is_hiv2):
        gene = isolate.gene
        # this method returns consensus or single sequence
        sequence = isolate.get_or_create_consensus()
        unusuals = []
        apobecs = []
        for pos, aas in sequence.aas:
            cons = cons_lookup[(gene, pos)][0]
            if '_' in aas:
                aas = '_'
            if len(aas) > 4:
                continue
            for aa in aas:
                key = (gene, pos, aa)
                if key in uum_lookup:
                    pcnt = uum_lookup[key]
                    unusuals.append('{}{}{} ({:.2f}%)'
                                    .format(cons, pos, aa, pcnt * 100))
                if key in apm_lookup:
                    apobecs.append('{}{}{}'.format(cons, pos, aa))
        yield {
            'IsolateID': isolate.id,
            'Gene': gene,
            '# Unusuals': len(unusuals),
            'Unusuals': ', '.join(unusuals),
            '# APOBECs': len(apobecs),
            'APOBECs': ', '.join(apobecs),
        }


@click.command()
@click.option('--aapcnt-json', type=click.File('r'), required=True)
@click.option('--apobec-json', type=click.File('r'), required=True)
@click.option('--filter', type=click.Choice(CRITERIA_CHOICES.keys()),
              multiple=True, default=('NO_CLONES', 'NO_QA_ISSUES'),
              show_default=True, help='specify filter criteria')
@click.option('--hiv2', is_flag=True, help='create table for HIV-2 sequences')
@click.argument('output_file', type=click.File('w'), default='-')
def export_seqqa(aapcnt_json, apobec_json, output_file, filter, hiv2):
    result = []
    aapcnt_data = json.load(aapcnt_json)
    cons_lookup = build_consensus_lookup(aapcnt_data)
    uum_lookup = unusual_mutation_lookup(aapcnt_data)
    apm_lookup = apobec_mutation_lookup(apobec_json)
    with app.app_context():
        for dc in ('PI', 'RTI', 'INSTI'):
            result.extend(
                run_seqqa(dc, filter, hiv2, cons_lookup,
                          uum_lookup, apm_lookup))
    writer = csv.DictWriter(output_file, CSV_HEADER)
    writer.writeheader()
    writer.writerows(result)


if __name__ == '__main__':
    export_seqqa()
