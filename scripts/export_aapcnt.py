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
MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG']
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
    'GENBANK_ONLY': Isolate.sequences.any(
        Sequence.accession.isnot(None) &
        (Sequence.accession != '')
    ),
    'NO_PARTIAL_MUTS': Isolate.sequences.any(
        Sequence.sequence_type == 'PartialMutationList'
    ),
}


def iter_isolates(drugclass, rx_type, criteria):
    print('Processing {} isolates ({})...'
          .format(drugclass, rx_type))
    gene = DRUG_CLASS_GENE_MAP[drugclass]

    conds = [CRITERIA_CHOICES[crkey] for crkey in criteria]
    query = (
        Isolate.query
        .filter(
            Isolate.gene == gene,
            Isolate.isolate_type == 'Clinical',
            Isolate._host.has(Host.host == 'Human'),
            Isolate._species.has(Species.species == 'HIV1'),
            Isolate._subtype.has(Subtype.subtype.notin_(
                ['O', 'N', 'P', 'CPZ']
            )),
            *conds
        )
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.insertions))
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.mixtures))
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
        print('  {}/{} isolates...'.format(offset, total), end='\r')
        yield from query.limit(QUERY_CHUNK_SIZE).offset(offset)
    print('  {0} isolates...                   '.format(total))


def stat_mutations(drugclass, rx_type, criteria, allow_mixture):
    counter = defaultdict(lambda: defaultdict(set))
    total_counter = defaultdict(lambda: defaultdict(set))
    rx_counter = defaultdict(set)
    ptids = set()
    subtype_counter = defaultdict(set)
    for isolate in iter_isolates(drugclass, rx_type, criteria):
        ptid = isolate.patient_id
        ptids.add(ptid)
        # this method returns consensus or single sequence
        sequence = isolate.get_or_create_consensus()
        subtype = isolate.subtype
        if subtype in MAJOR_SUBTYPES + ['D', 'F', 'G']:
            subtype_counter[subtype].add(ptid)
        else:
            subtype_counter['Others'].add(ptid)
        if rx_type == 'art':
            rx = isolate.total_rx
            if 'RAL' in rx.ii_list:
                rx_counter['RAL'].add(ptid)
            elif 'EVG' in rx.ii_list:
                rx_counter['EVG'].add(ptid)
            elif 'DTG' in rx.ii_list:
                rx_counter['DTG'].add(ptid)
            elif 'INI' in rx.ii_list:
                rx_counter['INI'].add(ptid)
        for pos, aas in sequence.aas:
            if '_' in aas:
                aas = '_'
            elif not allow_mixture and len(aas) > 1:
                # ignore mixtures
                continue
            for aa in aas:
                counter[(pos, aa)][subtype].add(ptid)
                total_counter[pos][subtype].add((ptid, aa))
                counter[(pos, aa)]['All'].add(ptid)
                total_counter[pos]['All'].add((ptid, aa))
                if subtype not in MAJOR_SUBTYPES:
                    counter[(pos, aa)]['Others'].add(ptid)
                    total_counter[pos]['Others'].add((ptid, aa))
    total = len(ptids)
    print('  {0} patients:'.format(total))
    inothersubtypes = set()
    for subtype in MAJOR_SUBTYPES + ['D', 'F', 'G', 'Others']:
        print('    Subtype {0}: {1} patients'
              .format(subtype, len(subtype_counter[subtype])))
        dups = inothersubtypes & subtype_counter[subtype]
        if dups:
            print(
                '    Following patients has multiple subtypes: {0}'
                .format(', '.join(str(d) for d in dups)))
        inothersubtypes |= subtype_counter[subtype]
    if rx_type == 'art':
        ralset = rx_counter['RAL']
        evgset = rx_counter['EVG']
        dtgset = rx_counter['DTG']
        iniset = rx_counter['INI']
        ralonly = len(ralset - evgset - dtgset - iniset)
        print('    RAL only: {0}/{1}'.format(ralonly, total))
        evgonly = len(evgset - ralset - dtgset - iniset)
        print('    EVG only: {0}/{1}'.format(evgonly, total))
        dtgonly = len(dtgset - ralset - evgset - iniset)
        print('    DTG only: {0}/{1}'.format(dtgonly, total))
        print('    Remains: {0}/{1}'.format(
            total - ralonly - evgonly - dtgonly, total))
    gene = DRUG_CLASS_GENE_MAP[drugclass]
    for pos in total_counter:
        totals = total_counter[pos]
        for aa in AMINO_ACIDS:
            counts = counter[(pos, aa)]
            for subtype in totals:
                total = len(totals[subtype])
                count = len(counts[subtype])
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
@click.option('--filter', type=click.Choice(CRITERIA_CHOICES.keys()),
              multiple=True, default=('NO_CLONES', 'NO_QA_ISSUES'),
              show_default=True, help='specify filter criteria')
@click.option('--allow-mixture', is_flag=True,
              help='the result should include mixtures')
@click.option('--format', type=click.Choice(['json', 'csv']),
              default='json', help='output format')
@click.option('--ugly-json', is_flag=True,
              help='output compressed (ugly) JSON instead of a pretty one')
def export_aapcnt(drugclass, output_file, filter,
                  allow_mixture, format, ugly_json):
    if drugclass != ('INSTI', ):
        raise NotImplementedError
    result = []
    with app.app_context():
        for dc in drugclass:
            for rt in ('naive', 'art'):
                result.extend(stat_mutations(dc, rt, filter, allow_mixture))
    if format == 'json':
        indent = None if ugly_json else '  '
        json.dump(result, output_file, indent=indent)
    elif format == 'csv':
        writer = csv.DictWriter(output_file, CSV_HEADER)
        writer.writeheader()
        writer.writerows(result)


if __name__ == '__main__':
    export_aapcnt()
