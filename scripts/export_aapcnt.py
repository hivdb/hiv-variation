#! /usr/bin/env python
import csv
import json
from itertools import groupby
from collections import defaultdict, Counter

import click
from hivdbql import app
from hivfacts import hivsdrm

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
GENE_RANGES = {
    'PR': (1, 100),
    'RT': (1, 561),
    'IN': (1, 289)
}
DRUG_CLASSES = ('PI', 'NRTI', 'NNRTI', 'NRTI2', 'NNRTI2', 'RTI', 'INSTI')
DRUG_CLASS_GENE_MAP = {
    'PI': 'PR',
    'NRTI': 'RT',
    'NNRTI': 'RT',
    'NRTI2': 'RT',
    'NNRTI2': 'RT',
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
        # NRTI2 & NNRTI2 have RTI-naive and NRTI/NNRTI-treated
        'NRTI2': (models.RTTotalRx.num_nrtis == 0,
                  models.RTTotalRx.cmp_num_nrtis == '=',
                  models.RTTotalRx.nrti_unknown.is_(False),
                  models.RTTotalRx.num_nnrtis == 0,
                  models.RTTotalRx.cmp_num_nnrtis == '=',
                  models.RTTotalRx.nnrti_unknown.is_(False)),
        'NNRTI2': (models.RTTotalRx.num_nrtis == 0,
                   models.RTTotalRx.cmp_num_nrtis == '=',
                   models.RTTotalRx.nrti_unknown.is_(False),
                   models.RTTotalRx.num_nnrtis == 0,
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
        'NRTI2': ((models.RTTotalRx.num_nrtis > 0) |
                  (models.RTTotalRx.nrti_list.like('%RTI%')),),
        'NNRTI2': ((models.RTTotalRx.num_nnrtis > 0) |
                   (models.RTTotalRx.nnrti_list.like('%RTI%')),),
        'RTI': ((models.RTTotalRx.num_nrtis > 0) |
                (models.RTTotalRx.nrti_list.like('%RTI%')) |
                (models.RTTotalRx.num_nnrtis > 0) |
                (models.RTTotalRx.nnrti_list.like('%RTI%')),),
        'INSTI': ((models.INTotalRx.num_iis > 0) |
                  (models.INTotalRx.ii_list.like('%INI%')),),
    },
    'truvada': {
        'PI': (False,),
        'NRTI': (models.RTTotalRx.num_nrtis == 2,
                 models.RTTotalRx.cmp_num_nrtis == '=',
                 models.RTTotalRx.drug_tdf == 'Yes',
                 (models.RTTotalRx.drug_ftc == 'Yes') |
                 (models.RTTotalRx.drug_3tc == 'Yes')),
        'NNRTI': (False,),
        'NRTI2': (models.RTTotalRx.num_nrtis == 2,
                  models.RTTotalRx.cmp_num_nrtis == '=',
                  models.RTTotalRx.drug_tdf == 'Yes',
                  (models.RTTotalRx.drug_ftc == 'Yes') |
                  (models.RTTotalRx.drug_3tc == 'Yes')),
        'NNRTI2': (False,),
        'RTI': (False,),
        'INSTI': (False,),
    },
}
MAJOR_SUBTYPES = ['A', 'B', 'C', 'CRF01_AE', 'CRF02_AG', 'D', 'F', 'G']
OTHER_SUBTYPES = [
    'CRF07_BC', 'CRF12_BF', 'CRF06_cpx', 'CRF08_BC', 'CRF35_AD',
    'CRF18_cpx', 'CRF45_cpx', 'CRF63_02A1'
]
MAJOR_SUBTYPES_HIV2 = ['HIV-2 A', 'HIV-2 B']
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
QUERY_CHUNK_SIZE = 20000

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
    'BEFORE_2009': Isolate.isolate_date < '2009-01-01',
    'AFTER_2009': Isolate.isolate_date >= '2009-01-01',
}

SDRM_LOOKUP = hivsdrm.HIVSDRM()


def load_subtype_lookup(fp):
    rows = csv.DictReader(fp)
    result = {}
    for row in rows:
        result[row['Accession']] = row['Subtype']
    return result


def iter_isolates(drugclass, rx_type, criteria, is_hiv2, exclude_ptids):
    click.echo('Processing {} isolates ({})...'
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
            # ~Isolate.references.any(RefLink.reference_id == 2350),
            *conds
        )
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.insertions))
        .options(db.selectinload(Isolate.sequences)
                 .selectinload(Sequence.mixtures))
    )
    if exclude_ptids:
        query = query.filter(
            Isolate.patient_id.notin_(exclude_ptids)
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
    query = query.order_by(Isolate.patient_id, Isolate.id)
    for offset in range(0, total, QUERY_CHUNK_SIZE):
        click.echo('  {}/{} isolates...\r'.format(offset, total), nl=False)
        yield from query.limit(QUERY_CHUNK_SIZE).offset(offset)
    click.echo('  {0} isolates...                   '.format(total))


def get_sdrms(isolates):
    sdrms = set()
    for isolate in isolates:
        gene = isolate.gene
        sequence = isolate.get_or_create_consensus()
        for pos, aas in sequence.aas:
            if '_' in aas:
                aas = '_'
            if len(aas) > 4:
                continue
            sdrms |= SDRM_LOOKUP.get_sdrms(gene, pos, aas)
    return sdrms


def find_exclude_ptids(rx, drugclasses, criteria, is_hiv2,
                       max_sdrm_per_naive_person):
    sdrms_per_patient = defaultdict(set)
    for drugclass in drugclasses:
        isolates_by_ptid = groupby(
            iter_isolates(drugclass, rx, criteria, is_hiv2, None),
            lambda isolate: isolate.patient_id
        )
        for ptid, isolates in isolates_by_ptid:
            isolates = list(isolates)
            sdrms_per_patient[ptid] |= get_sdrms(isolates)
    exclude_ptids = []
    for ptid, sdrms in sdrms_per_patient.items():
        if len(sdrms) > max_sdrm_per_naive_person:
            exclude_ptids.append(ptid)
    return exclude_ptids


def stat_mutations(drugclass, rx_type, criteria, is_hiv2,
                   denominator, allow_mixture, seq_subtypes,
                   exclude_ptids):
    denom_use_patients = denominator == 'patients'
    counter = defaultdict(lambda: defaultdict(set))
    total_counter = defaultdict(lambda: defaultdict(set))
    # rx_counter = defaultdict(set)
    ptids = set()
    subtype_counter = Counter()
    major_subtypes = MAJOR_SUBTYPES_HIV2 if is_hiv2 else MAJOR_SUBTYPES
    other_subtypes = [] if is_hiv2 else OTHER_SUBTYPES
    report_subtypes = set(major_subtypes + other_subtypes)
    gene = DRUG_CLASS_GENE_MAP[drugclass]
    isolates_by_ptid = groupby(
        iter_isolates(drugclass, rx_type, criteria, is_hiv2, exclude_ptids),
        lambda isolate: isolate.patient_id
    )
    for ptid, isolates in isolates_by_ptid:
        ptids.add(ptid)
        isolates = list(isolates)
        sdrms = get_sdrms(isolates)
        for isolate in isolates:
            # this method returns consensus or single sequence
            sequence = isolate.get_or_create_consensus()
            subtype_key = (sequence.accession or
                           'PT{}'.format(isolate.patient_id))
            if sequence.sequence_type == 'Consensus':
                for seq in isolate.sequences:
                    if seq.sequence_type == 'Sequence':
                        subtype_key = seq.accession
                        break
            subtype = seq_subtypes.get(subtype_key, isolate.subtype)
            subtype = subtype and subtype.replace('Group ', '')
            if subtype in report_subtypes:
                subtype_counter[subtype] += 1
            else:
                subtype_counter[subtype] += 1
                subtype_counter['Others'] += 1
            # # TODO: print number of patients for each drug
            # if not is_hiv2 and rx_type == 'art':
            #     rx = isolate.total_rx
            #     if 'RAL' in rx.ii_list:
            #         rx_counter['RAL'].add(ptid)
            #     if 'EVG' in rx.ii_list:
            #         rx_counter['EVG'].add(ptid)
            #     if 'DTG' in rx.ii_list:
            #         rx_counter['DTG'].add(ptid)
            #     if 'INI' in rx.ii_list:
            #         rx_counter['INI'].add(ptid)
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
                    counter[(pos, aa)]['All'].add(ptid)
                    total_counter[pos]['All'].add(totalkey)
                    if len(sdrms - {(gene, pos, aa)}) > 0:
                        counter[(pos, aa)]['WithSDRM'].add(ptid)
                    else:
                        counter[(pos, aa)]['WithoutSDRM'].add(ptid)
                    if subtype in report_subtypes:
                        counter[(pos, aa)][subtype].add(ptid)
                        total_counter[pos][subtype].add(totalkey)
                    else:
                        counter[(pos, aa)]['Others'].add(ptid)
                        total_counter[pos]['Others'].add(totalkey)
    # if rx_type == 'naive':
    #     for key, value in counter.items():
    #         if value['WithSDRM']:
    #             click.echo("Naive With SDRMs: {}, {}"
    #                        .format(repr(key), repr(value['WithSDRM'])))
    total = len(ptids)
    click.echo('  {0} {1}:'.format(total, denominator))
    for subtype in report_subtypes:
        click.echo('    Subtype {0}: {1} isolates'
                   .format(subtype, subtype_counter[subtype]))
    # if not is_hiv2 and rx_type == 'art':
    #     ralset = rx_counter['RAL']
    #     evgset = rx_counter['EVG']
    #     dtgset = rx_counter['DTG']
    #     iniset = rx_counter['INI']
    #     ralonly = len(ralset - evgset - dtgset - iniset)
    #     print('    RAL only: {0}/{1}'.format(ralonly, total))
    #     evgonly = len(evgset - ralset - dtgset - iniset)
    #     print('    EVG only: {0}/{1}'.format(evgonly, total))
    #     dtgonly = len(dtgset - ralset - evgset - iniset)
    #     print('    DTG only: {0}/{1}'.format(dtgonly, total))
    #     knownrx = len(ralset | evgset | dtgset)
    #     print('    Combined: {0}/{1}'.format(
    #         knownrx - ralonly - evgonly - dtgonly, total))
    #     unknownrx = len(iniset - ralset - evgset - dtgset)
    #     print('    Unknown: {0}/{1}'.format(
    #         unknownrx, total))
    for pos in range(*GENE_RANGES[gene]):
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
                    'total': total,
                    'with_sdrm': len(counts['WithSDRM']),
                    'without_sdrm': len(counts['WithoutSDRM']),
                }


@click.command()
@click.argument('drugclass', nargs=-1, required=True,
                type=click.Choice(DRUG_CLASSES))
@click.argument('output_file', type=click.File('w'))
@click.option('--species', type=str, default='HIV1',
              help='specify an HIV species')
@click.option('--filter', type=click.Choice(CRITERIA_CHOICES.keys()),
              multiple=True, default=('NO_CLONES', 'NO_QA_ISSUES'),
              show_default=True, help='specify filter criteria')
@click.option('--no-filter', is_flag=True,
              help='Don\'t apply any extra filter to the query')
@click.option('--rx-type', type=click.Choice(['art', 'naive',
                                              'all', 'truvada']),
              multiple=True, default=['art', 'naive', 'all'],
              help='specify treatment type')
@click.option('--allow-mixture', is_flag=True,
              help='the result should include mixtures')
@click.option('--denominator', default='patients',
              type=click.Choice(['patients', 'patient-variants']),
              help='denominator when calculating the prevalence')
@click.option('--subtype-lookup', type=click.File('r'),
              help=('use extra subtype lookup table '
                    'instead of database subtype'))
@click.option('--max-sdrm-per-naive-person', type=int,
              help='Maximum allowed SDRMs found in naive patient')
@click.option('--format', type=click.Choice(['json', 'csv']),
              default='json', help='output format')
@click.option('--ugly-json', is_flag=True,
              help='output compressed (ugly) JSON instead of a pretty one')
def export_aapcnt(drugclass, output_file, species, filter, no_filter,
                  rx_type, allow_mixture, denominator, subtype_lookup,
                  max_sdrm_per_naive_person, format, ugly_json):
    result = []
    if no_filter:
        filter = []
    seq_subtypes = {}
    if subtype_lookup:
        seq_subtypes = load_subtype_lookup(subtype_lookup)
    hiv2 = species == 'HIV2'
    with app.app_context():
        exclude_ptids = []
        if 'naive' in rx_type and max_sdrm_per_naive_person is not None:
            exclude_ptids = find_exclude_ptids(
                'naive', drugclass, filter, hiv2, max_sdrm_per_naive_person)
        for rx in rx_type:
            for dc in drugclass:
                result.extend(
                    stat_mutations(dc, rx, filter, hiv2, denominator,
                                   allow_mixture, seq_subtypes,
                                   exclude_ptids if rx == 'naive' else []))
    if format == 'json':
        indent = None if ugly_json else '  '
        json.dump(result, output_file, indent=indent)
    elif format == 'csv':
        writer = csv.DictWriter(output_file, CSV_HEADER, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(result)


if __name__ == '__main__':
    export_aapcnt()
