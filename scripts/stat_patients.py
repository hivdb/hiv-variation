#! /usr/bin/env python
import csv
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


def stat_patients(drugclass, rx_type, criteria, is_hiv2,
                  exclude_ptids):
    # gene = DRUG_CLASS_GENE_MAP[drugclass]
    patient_counter = defaultdict(set)
    isolate_counter = Counter()
    major_subtypes = MAJOR_SUBTYPES_HIV2 if is_hiv2 else MAJOR_SUBTYPES
    other_subtypes = [] if is_hiv2 else OTHER_SUBTYPES
    report_subtypes = set(major_subtypes + other_subtypes)

    isolates_by_ptid = groupby(
        iter_isolates(drugclass, rx_type, criteria, is_hiv2, exclude_ptids),
        lambda isolate: isolate.patient_id
    )
    for ptid, isolates in isolates_by_ptid:
        patient_counter['All'].add(ptid)
        for isolate in isolates:
            subtype = isolate.subtype
            if subtype in report_subtypes:
                patient_counter[subtype].add(ptid)
                isolate_counter[subtype] += 1
            else:
                patient_counter['Other'].add(ptid)
                isolate_counter['Other'] += 1
            isolate_counter['All'] += 1
    for subtype in ['All'] + major_subtypes + other_subtypes + ['Other']:
        yield {
            'Category': 'Patient',
            'Subtype': subtype,
            'DrugClass': drugclass,
            'RxType': rx_type,
            'Count':  len(patient_counter[subtype])
        }

    for subtype in ['All'] + major_subtypes + other_subtypes + ['Other']:
        yield {
            'Category': 'Isolate',
            'Subtype': subtype,
            'DrugClass': drugclass,
            'RxType': rx_type,
            'Count':  isolate_counter[subtype]
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
@click.option('--subtype-lookup', type=click.File('r'),
              help=('use extra subtype lookup table '
                    'instead of database subtype'))
@click.option('--max-sdrm-per-naive-person', type=int,
              help='Maximum allowed SDRMs found in naive patient')
def export_aapcnt(drugclass, output_file, species, filter, no_filter,
                  rx_type, subtype_lookup,
                  max_sdrm_per_naive_person):
    if no_filter:
        filter = []
    hiv2 = species == 'HIV2'
    with app.app_context():
        writer = csv.DictWriter(
            output_file,
            ['Category', 'Subtype', 'DrugClass', 'RxType', 'Count']
        )
        exclude_ptids = []
        if 'naive' in rx_type and max_sdrm_per_naive_person is not None:
            exclude_ptids = find_exclude_ptids(
                'naive', drugclass, filter, hiv2, max_sdrm_per_naive_person)
        for rx in rx_type:
            for dc in drugclass:
                writer.writerows(
                    stat_patients(dc, rx, filter, hiv2,
                                  exclude_ptids if rx == 'naive' else []))


if __name__ == '__main__':
    export_aapcnt()
