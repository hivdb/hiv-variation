#! /usr/bin/env python
import csv
import json
from collections import defaultdict, Counter

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
RefLink = models.RefLink

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
               (models.PRTotalRx.pi_list.like('%PI%')),),
        'NRTI': ((models.RTTotalRx.num_nrtis > 0) |
                 (models.RTTotalRx.nrti_list.like('%RTI%')),),
        'NNRTI': ((models.RTTotalRx.num_nnrtis > 0) |
                  (models.RTTotalRx.nnrti_list.like('%RTI%')),),
        'RTI': ((models.RTTotalRx.num_nrtis > 0) |
                (models.RTTotalRx.nrti_list.like('%RTI%')) |
                (models.RTTotalRx.num_nnrtis > 0) |
                (models.RTTotalRx.nnrti_list.like('%RTI%')),),
        'INSTI': ((models.INTotalRx.num_iis > 0) |
                  (models.INTotalRx.ii_list.like('%INI%')),),
    },
}
MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG', 'D', 'F', 'G']
MAJOR_SUBTYPES_HIV2 = ['A', 'B']
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY_-'
# UNUSUAL_CUTOFF = 0.0001  # 1 in 10,000 or 0.01%
CSV_HEADER = [
    'gene',
    'position',
    'subtype',
    'aa',
    'rx_type',
    'percent',
    'count',
    'total'
]
QUERY_CHUNK_SIZE = 500

CRITERIA_CHOICES = {
    'HIV1_ONLY': Isolate._species.has(Species.species == 'HIV1'),
    'HIV2_ONLY': Isolate._species.has(Species.species == 'HIV2'),
    'PLASMA_ONLY': Isolate.clinical_isolate.has(
        ClinicalIsolate.source == 'Plasma'),
    'SANGER_ONLY': Isolate.clinical_isolate.has(
        ClinicalIsolate.seq_method == 'Dideoxy'),
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
    'PUBLISHED': Isolate.references.any(
        RefLink.reference.has(Reference.published.is_(True))
    ),
}


def load_subtype_lookup(fp):
    rows = csv.DictReader(fp)
    result = {}
    for row in rows:
        result[row['Accession']] = row['Subtype']
    return result


def iter_isolates(drugclass, rx_type, criteria, is_hiv2):
    print('Processing {} isolates ({})...'
          .format(drugclass, rx_type))
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
            # Exclude Damond 2005
            ~Isolate.references.any(RefLink.reference_id == 2350),
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


def stat_mutations(drugclass, rx_type, criteria, is_hiv2,
                   denominator, allow_mixture, seq_subtypes):
    denom_use_patients = denominator == 'patients'
    counter = defaultdict(lambda: defaultdict(set))
    total_counter = defaultdict(lambda: defaultdict(set))
    rx_counter = defaultdict(set)
    ptids = set()
    subtype_counter = Counter()
    major_subtypes = MAJOR_SUBTYPES_HIV2 if is_hiv2 else MAJOR_SUBTYPES
    gene = DRUG_CLASS_GENE_MAP[drugclass]
    for isolate in iter_isolates(drugclass, rx_type, criteria, is_hiv2):
        ptid = isolate.patient_id
        ptids.add(ptid)
        # this method returns consensus or single sequence
        sequence = isolate.get_or_create_consensus()
        subtype_key = sequence.accession or 'PT{}'.format(isolate.patient_id)
        if sequence.sequence_type == 'Consensus':
            for seq in isolate.sequences:
                if seq.sequence_type == 'Sequence':
                    subtype_key = seq.accession
                    break
        subtype = seq_subtypes.get(subtype_key, isolate.subtype)
        subtype = subtype and subtype.replace('Group ', '')
        if subtype in major_subtypes:
            subtype_counter[subtype] += 1
        else:
            subtype_counter['Others'] += 1
        # TODO: print number of patients for each drug
        if not is_hiv2 and rx_type == 'art':
            rx = isolate.total_rx
            if 'RAL' in rx.ii_list:
                rx_counter['RAL'].add(ptid)
            if 'EVG' in rx.ii_list:
                rx_counter['EVG'].add(ptid)
            if 'DTG' in rx.ii_list:
                rx_counter['DTG'].add(ptid)
            if 'INI' in rx.ii_list:
                rx_counter['INI'].add(ptid)
        for pos, aas in sequence.aas:
            if '_' in aas:
                aas = '_'
            elif len(aas) > 4:
                # ignore any high ambiguous position
                continue
            elif not allow_mixture and len(aas) > 1:
                # ignore mixtures
                continue
            for aa in aas:
                if aa == '*':
                    continue
                totalkey = ptid if denom_use_patients else (ptid, aa)
                counter[(pos, aa)][subtype].add(ptid)
                total_counter[pos][subtype].add(totalkey)
                counter[(pos, aa)]['All'].add(ptid)
                total_counter[pos]['All'].add(totalkey)
                if subtype not in major_subtypes:
                    counter[(pos, aa)]['Others'].add(ptid)
                    total_counter[pos]['Others'].add(totalkey)
    total = len(ptids)
    print('  {0} {1}:'.format(total, denominator))
    for subtype in major_subtypes + ['Others']:
        print('    Subtype {0}: {1} isolates'
              .format(subtype, subtype_counter[subtype]))
    if not is_hiv2 and rx_type == 'art':
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
        knownrx = len(ralset | evgset | dtgset)
        print('    Combined: {0}/{1}'.format(
            knownrx - ralonly - evgonly - dtgonly, total))
        unknownrx = len(iniset - ralset - evgset - dtgset)
        print('    Unknown: {0}/{1}'.format(
            unknownrx, total))
    for pos in sorted(total_counter.keys()):
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
@click.option('--no-filter', is_flag=True,
              help='Don\'t apply any extra filter to the query')
@click.option('--allow-mixture', is_flag=True,
              help='the result should include mixtures')
@click.option('--denominator', default='patients',
              type=click.Choice(['patients', 'patient-variants']),
              help='denominator when calculating the prevalence')
@click.option('--subtype-lookup', type=click.File('r'))
@click.option('--hiv2', is_flag=True, help='create table for HIV-2 sequences')
@click.option('--format', type=click.Choice(['json', 'csv']),
              default='json', help='output format')
@click.option('--ugly-json', is_flag=True,
              help='output compressed (ugly) JSON instead of a pretty one')
def export_aapcnt(drugclass, output_file, filter, no_filter,
                  allow_mixture, denominator, subtype_lookup,
                  hiv2, format, ugly_json):
    result = []
    if no_filter:
        filter = []
    seq_subtypes = {}
    if subtype_lookup:
        seq_subtypes = load_subtype_lookup(subtype_lookup)
    with app.app_context():
        for dc in drugclass:
            for rt in ('all', 'naive', 'art'):
                result.extend(
                    stat_mutations(dc, rt, filter, hiv2, denominator,
                                   allow_mixture, seq_subtypes))
    if format == 'json':
        indent = None if ugly_json else '  '
        json.dump(result, output_file, indent=indent)
    elif format == 'csv':
        writer = csv.DictWriter(output_file, CSV_HEADER, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(result)


if __name__ == '__main__':
    export_aapcnt()
