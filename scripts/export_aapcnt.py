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

GENES = ('PR', 'RT', 'IN')
DRUG_CLASSES = ('PI', 'NRTI', 'NNRTI', 'RTI', 'INSTI')
DRUG_CLASS_GENE_MAP = {
    'PI': 'PR',
    'NRTI': 'RT',
    'NNRTI': 'RT',
    'RTI': 'RT',
    'INSTI': 'IN'
}
DRUG_CLASS_RX_QUERIES = {
    'naive': {
        'PI': (models.PRTotalRx.num_pis == 0,
               models.PRTotalRx.cmp_num_pis == '=',
               models.PRTotalRx.pi_unknown.is_(False)),
        'NRTI': (models.RTTotalRx.num_nrtis == 0,
                 models.RTTotalRx.cmp_num_nrtis == '=',
                 models.RTTotalRx.nrti_unknown.is_(False)),
        'NNRTI': (models.RTTotalRx.num_nnrtis == 0,
                  models.RTTotalRx.cmp_num_nnrtis == '=',
                  models.RTTotalRx.nnrti_unknown.is_(False)),
        'RTI': (models.RTTotalRx.num_nrtis == 0,
                models.RTTotalRx.cmp_num_nrtis == '=',
                models.RTTotalRx.nrti_unknown.is_(False),
                models.RTTotalRx.num_nnrtis == 0,
                models.RTTotalRx.cmp_num_nnrtis == '=',
                models.RTTotalRx.nnrti_unknown.is_(False)),
        'INSTI': (models.INTotalRx.num_iis == 0,
                  models.INTotalRx.cmp_num_iis == '=',
                  models.INTotalRx.ii_unknown.is_(False)),
    },
    'art': {
        'PI': ((models.PRTotalRx.num_pis > 0) |
               (models.PRTotalRx.pi_list.like('%PI%'))),
        'NRTI': ((models.RTTotalRx.num_nrtis > 0) |
                 (models.RTTotalRx.nrti_list.like('%NRTI%'))),
        'NNRTI': ((models.RTTotalRx.num_nnrtis > 0) |
                  (models.RTTotalRx.nnrti_list.like('%NNRTI%'))),
        'RTI': ((models.RTTotalRx.num_nrtis > 0) |
                (models.RTTotalRx.nrti_list.like('%NRTI%')) |
                (models.RTTotalRx.num_nnrtis > 0) |
                (models.RTTotalRx.nnrti_list.like('%NNRTI%'))),
        'INSTI': ((models.INTotalRx.num_iis > 0) |
                  (models.INTotalRx.ii_list.like('%INI%')),),
    },
}
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY_-*'
# UNUSUAL_CUTOFF = 0.0001  # 1 in 10,000 or 0.01%
CSV_HEADER = [
    'gene',
    'position',
    'aa',
    'percent',
    'count',
    'total'
]
QUERY_CHUNK_SIZE = 500

CRITERIA_CHOICES = {
    'PLASMA_ONLY': Isolate.clinical_isolate.has(
        ClinicalIsolate.source == 'Plasma'),
    'NO_CLONES': Isolate.clinical_isolate.has(
        ClinicalIsolate.clone_method == 'None'),
    'NO_QA_ISSUES': ~Isolate._filter.has(),
    'GENBANK_ONLY': Isolate.references.any(
        Reference.medline_id.isnot(None) &
        (Reference.medline_id != '')
    ),
    'NO_PARTIAL_MUTS': Isolate.sequences.any(
        Sequence.sequence_type == 'PartialMutationList'
    ),
}


def iter_isolates(drugclass, subtype, rx_type, criteria):
    print('Processing {} isolates ({}/{})...'
          .format(drugclass, subtype, rx_type))
    gene = DRUG_CLASS_GENE_MAP[drugclass]

    conds = [CRITERIA_CHOICES[crkey] for crkey in criteria]
    query = (
        Isolate.query
        .filter(
            Isolate.gene == gene,
            Isolate.isolate_type == 'Clinical',
            Isolate._host.has(Host.host == 'Human'),
            Isolate._species.has(Species.species == 'HIV1'),
            *conds
        )
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.insertions))
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.mixtures))
    )

    if subtype == 'All':
        query = query.filter(
            Isolate._subtype.has(Subtype.subtype.notin_(
                ['NA', 'U', 'O', 'N', 'P', 'CPZ', 'Unknown']
            ))
        )
    else:
        query = query.filter(
            Isolate._subtype.has(Subtype.subtype == subtype)
        )
    if rx_type != 'all':
        rx_query = DRUG_CLASS_RX_QUERIES[rx_type][drugclass]
        attrname = '_{}_total_rx'.format(gene.lower())
        query = query.filter(Isolate.clinical_isolate.has(
            getattr(ClinicalIsolate, attrname).has(db.and_(*rx_query))
        ))

    total = query.count()
    query = query.order_by(Isolate.id)
    for offset in range(0, total, QUERY_CHUNK_SIZE):
        print('  {}/{}...'.format(offset, total), end='\r')
        yield from query.limit(QUERY_CHUNK_SIZE).offset(offset)
    print('  {0}/{0}...'.format(total))


def stat_mutations(drugclass, subtype, rx_type, criteria):
    counter = defaultdict(set)
    total_counter = defaultdict(set)
    for isolate in iter_isolates(drugclass, subtype, rx_type, criteria):
        ptid = isolate.patient_id
        # if isolate.filter_reason:
        #     print('Warning: Isolate {} has filter_reason'.format(isolate.id))
        # this method returns consensus or single sequence
        sequence = isolate.get_or_create_consensus()
        for pos, aa in sequence.aas:
            if '_' in aa:
                aa = '_'
            elif len(aa) > 1:
                # ignore mixtures
                continue
            counter[(pos, aa)].add(ptid)
            total_counter[pos].add((ptid, aa))
    gene = DRUG_CLASS_GENE_MAP[drugclass]
    for pos in total_counter:
        total = len(total_counter[pos])
        for aa in AMINO_ACIDS:
            count = len(counter[(pos, aa)])
            yield {
                'gene': gene,
                'drugclass': drugclass,
                'subtype': subtype,
                'rx_type': rx_type,
                'position': pos,
                'aa': aa,
                'percent': count / total,
                'count': count,
                'total': total
            }


@click.command()
@click.argument('drugclass', nargs=-1, required=True,
                type=click.Choice(DRUG_CLASSES))
@click.argument('output_file', type=click.File('w'))
@click.option('--subtype', type=str, multiple=True,
              default=('All', 'A', 'B', 'C', 'CRF01_AE', 'CRF02_AG'),
              show_default=True, help='specify an HIV subtype/group')
@click.option('--rx-type', type=click.Choice(('naive', 'art', 'all')),
              multiple=True, default=('naive', 'art', 'all'),
              show_default=True, help='specify treatment type')
@click.option('--filter', type=click.Choice(CRITERIA_CHOICES.keys()),
              multiple=True, default=(), show_default=True,
              help='specify filter criteria')
@click.option('--format', type=click.Choice(['json', 'csv']),
              default='json', help='output format')
@click.option('--ugly-json', is_flag=True,
              help='output compressed (ugly) JSON instead of a pretty one')
def export_aapcnt(drugclass, output_file, subtype,
                  rx_type, filter, format, ugly_json):
    if drugclass != ('INSTI', ):
        raise NotImplementedError
    result = []
    with app.app_context():
        for dc in drugclass:
            for st in subtype:
                for rt in rx_type:
                    result.extend(stat_mutations(dc, st, rt, filter))
    if format == 'json':
        indent = None if ugly_json else '  '
        json.dump(result, output_file, indent=indent)
    elif format == 'csv':
        writer = csv.DictWriter(output_file, CSV_HEADER)
        writer.writeheader()
        writer.writerows(result)


if __name__ == '__main__':
    export_aapcnt()
