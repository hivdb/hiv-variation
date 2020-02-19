#! /usr/bin/env python
import sys
import csv
from itertools import groupby

from hivdbql import app
from hivdbql.utils import dbutils

from hypermut2 import hypermut

from hivdbql.utils.codonutils import translate_codon

db = app.db
models = app.models

FIELDS = [
    'PtID',
    'IsolateDate',
    'Genes',
    'AccessionID',
    'Source',
    'SeqMethod',
    'Hypermut PositiveAPOBECs',
    'Hypermut PotentialAPOBECs',
    'Hypermut PositiveControls',
    'Hypermut PotentialControls',
    'Hypermut Ratio (A/B)/(C/D)',
    'Adjusted Hypermut Ratio (+1)',
    'Hypermut P-value (Fisher Exact)',
    'HIVDB APOBEC Count',
    'Stop Codon Count',
    'HIVDB APOBECs'
]

GENE_NACONS = [
    'CCTCAAATCACTCTTTGGCAGCGACCCCTTGTCACAATAAAAATAGGGGGACAGCTAAGGGAAGCTCTATTAG'
    'ATACAGGAGCAGATGATACAGTATTAGAAGAAATAAATTTGCCAGGAAAATGGAAACCAAAAATGATAGGGGG'
    'AATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAAATACTCATAGAAATTTGTGGAAAAAAGGCTATAGGT'
    'ACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATATGTTGACTCAGATTGGTTGCACTTTAA'
    'ATTTT',
    'CCAATTAGTCCTATTGAAACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAGGTTAAACAATGGC'
    'CATTGACAGAAGAAAAAATAAAAGCATTAACAGAAATTTGTACAGAAATGGAAAAGGAAGGAAAAATTTCAAA'
    'AATTGGGCCTGAAAATCCATACAACACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAA'
    'TTAGTAGATTTCAGAGAACTCAATAAAAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCAG'
    'CAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGGGATGCATATTTTTCAGTTCCTTTAGATGA'
    'AGACTTCAGGAAGTATACTGCATTCACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTAC'
    'AATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAGAGTAGCATGACAAAAATCTTAGAGCCCT'
    'TTAGAACAAAAAATCCAGAAATAGTTATCTACCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAAT'
    'AGGGCAGCATAGAGCAAAAATAGAGGAGTTAAGAGAACATCTATTGAGGTGGGGATTTACCACACCAGACAAA'
    'AAACATCAGAAAGAACCTCCATTTCTTTGGATGGGATATGAACTCCATCCTGACAAATGGACAGTACAGCCTA'
    'TACAGCTGCCAGAAAAAGACAGCTGGACTGTCAATGATATACAGAAGTTAGTGGGAAAACTAAATTGGGCAAG'
    'TCAGATTTATCCAGGGATTAAAGTAAAGCAACTATGTAAACTCCTTAGGGGAGCCAAAGCACTAACAGACATA'
    'GTACCACTGACTGAAGAAGCAGAATTAGAATTGGCAGAGAACAGGGAGATTCTAAAAGAACCAGTACATGGAG'
    'TATATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAACAAGGGCAAGACCAATGGACATATCAAAT'
    'TTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAGTATGCAAAAATGAGGTCTGCCCACACTAATGATGTA'
    'AAACAATTAACAGAAGCAGTGCAAAAAATAGCCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTA'
    'GACTACCCATACAAAAAGAAACATGGGAGACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTG'
    'GGAGTTTGTCAATACCCCTCCTCTAGTAAAATTATGGTACCAGTTAGAAAAAGAACCCATAGTAGGAGCAGAA'
    'ACTTTCTATGTAGATGGGGCAGCTAATAGGGAGACTAAACTAGGAAAAGCAGGATATGTTACTGACAGAGGAA'
    'GACAAAAAGTTGTTTCCCTAACTGAAACAACAAATCAGAAGACTGAATTACAAGCAATTCATCTAGCTTTGCA'
    'GGATTCAGGATCAGAAGTAAACATAGTAACAGACTCACAGTATGCATTAGGAATCATTCAAGCACAACCAGAT'
    'AAGAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTACCTGTCATGGG'
    'TACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTTCTGGAATCAGGAAAGTGCT'
    'A',
    'TTTTTAGATGGGATAGATAAGGCTCAAGAAGAACATGAAAAATATCACAGCAATTGGAGAGCAATGGCTAGTG'
    'ATTTTAATCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCTAGCTGTGATAAATGTCAGCTAAAAGGGGAAGC'
    'CATGCATGGACAAGTAGACTGTAGTCCAGGGATATGGCAATTAGATTGTACACATTTAGAAGGAAAAGTTATC'
    'CTGGTAGCAGTCCATGTAGCCAGTGGCTATATAGAAGCAGAAGTTATCCCAGCAGAAACAGGACAGGAAACAG'
    'CATACTTTATATTAAAATTAGCAGGAAGATGGCCAGTAAAAGTAATACATACAGACAATGGCAGCAATTTCAC'
    'CAGTGCTGCGGTTAAGGCAGCCTGTTGGTGGGCAGGTATCCAGCAGGAATTTGGAATTCCCTACAATCCCCAA'
    'AGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATCATAGGGCAGGTAAGAGATCAAGCTGAAC'
    'ACCTTAAGACAGCAGTACAAATGGCAGTATTCATTCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAG'
    'TGCAGGGGAAAGAATAATAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAA'
    'ATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAGACCCAATTTGGAAAGGACCAGCAAAACTACTCTGGA'
    'AAGGTGAAGGGGCAGTAGTAATACAAGACAATAGTGAAATAAAGGTAGTACCAAGAAGAAAAGCAAAGATCAT'
    'TAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAGGTAGACAGGATGAGGAT'
]

GENE_SIZES = [99, 560, 288]

GENE_ORDER = {
    'PR': 0,
    'RT': 1,
    'IN': 2
}

DRUG_CLASSES = ['PI', 'RTI', 'INSTI']

ACTIVE_SITES = [
    ('PR', 25, 'N'),
    ('RT', 110, 'N'),
    ('RT', 185, 'N'),
    ('RT', 186, 'N'),
    ('IN', 64, 'N'),
    ('IN', 116, 'N'),
    ('IN', 152, 'K')
]


APOBEC_MUTATION_MAP = models.LUAPOBEC.apobec_lookup_table()


def iter_isolates():
    Isolate = models.Isolate
    query = Isolate.make_query(
        'HIV1', None, None, [
            'NO_CLONES',
            'PLASMA_OR_PBMC',
            'SANGER_ONLY',
            'NASEQ_AVAILABLE'
        ])
    query = query.options(db.selectinload('sequences'))
    query = query.options(db.joinedload('clinical_isolate'))
    isolates = dbutils.chunk_query(
        query, (
            Isolate.patient_id,
            Isolate.isolate_date,
            Isolate.gene
        ), chunksize=20000,
        on_progress=(lambda o, t:
                     print('{0}/{1} isolates...'.format(o, t), end='\r')),
        on_finish=(lambda t:
                   print('{0} isolates.                        '.format(t)))
    )
    isolates = groupby(
        isolates, lambda iso: (iso.patient_id, iso.isolate_date))
    for (ptid, date), isolates in isolates:
        isolates = sorted(isolates, key=lambda iso: GENE_ORDER[iso.gene])
        genes = ''.join(iso.gene for iso in isolates)
        if genes not in ('PRRT', 'PRRTIN', 'IN'):
            continue
        source = ','.join(sorted({
            iso.clinical_isolate.source
            for iso in isolates
        }))
        if ',' in source:
            continue
        yield isolates, genes, source


def compare(isolates):
    all_naseq = ''
    all_refseq = ''
    hivdb_result = [0, 0, []]
    for isolate in isolates:
        seq = isolate.get_or_create_consensus()

        gene = isolate.gene
        geneidx = GENE_ORDER[gene]
        gene_size = GENE_SIZES[geneidx]

        naseq = seq.naseq
        first_aa = seq.first_aa
        last_aa = seq.last_aa
        if naseq is None:
            print(seq)
        naseq = ('...' * (first_aa - 1) +
                 naseq +
                 '...' * (gene_size - last_aa))

        refseq = GENE_NACONS[geneidx]

        all_naseq += naseq
        all_refseq += refseq

        for pos, aas in seq.aas:
            if len(aas) > 4:
                # ignore highly ambiguous sites
                continue
            refaa = translate_codon(refseq[pos * 3 - 3:pos * 3])
            if any((gene, pos, aa) in APOBEC_MUTATION_MAP
                   for aa in aas):
                hivdb_result[0] += 1
                hivdb_result[2].append('{}{}{}'.format(refaa, pos, aas))
            if '*' in aas:
                hivdb_result[1] += 1
            if any((gene, pos, aa) in ACTIVE_SITES for aa in aas):
                hivdb_result[1] += 1

    hypermut_result = hypermut(all_refseq, all_naseq)
    return hypermut_result, hivdb_result


def compare_all():
    for isolates, genes, source in iter_isolates():
        hypermut_result, hivdb_result = compare(isolates)
        isolate0 = isolates[0]
        yield {
            'PtID': isolate0.patient_id,
            'IsolateDate': isolate0.isolate_date,
            'Genes': genes,
            'AccessionID': ','.join(sorted({
                iso.get_or_create_consensus().accession
                for iso in isolates
                if iso.get_or_create_consensus().accession
            })),
            'Source': source,
            'SeqMethod': ','.join(sorted({
                iso.clinical_isolate.seq_method
                for iso in isolates
            })),
            'Hypermut PositiveAPOBECs': hypermut_result[0],
            'Hypermut PotentialAPOBECs': hypermut_result[1],
            'Hypermut PositiveControls': hypermut_result[2],
            'Hypermut PotentialControls': hypermut_result[3],
            'Hypermut Ratio (A/B)/(C/D)': hypermut_result[4],
            'Adjusted Hypermut Ratio (+1)': hypermut_result[5],
            'Hypermut P-value (Fisher Exact)': hypermut_result[6],
            'HIVDB APOBEC Count': hivdb_result[0],
            'Stop Codon Count': hivdb_result[1],
            'HIVDB APOBECs': ', '.join(hivdb_result[2])
        }


def main():
    if len(sys.argv) != 2:
        print('Usage: {} <ALL_OUTPUT>'.format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    allout = sys.argv[1]
    with open(allout, 'w') as allout:
        all_writer = csv.DictWriter(allout, FIELDS)
        all_writer.writeheader()
        for result in compare_all():
            all_writer.writerow(result)


if __name__ == '__main__':
    with app.app_context():
        main()
