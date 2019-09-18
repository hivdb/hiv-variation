#! /usr/bin/env python
import sys
import csv

from hivdbql import app
from hivdbql.utils import dbutils

from hypermut2 import hypermut

from hivdbql.utils.codonutils import translate_codon

db = app.db
models = app.models

FIELDS = [
    'Gene',
    'IsolateID',
    'PtID',
    'IsolateDate',
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
    'HIVDB 2019APOBEC Count',
    'Stop Codon Count',
    'HIVDB APOBECs',
    'HIVDB 2019APOBECs',
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


def apobec_mutation_map():
    apobecs = models.LUAPOBEC.query.filter(
        models.LUAPOBEC.is_apobec.is_(True)
    )
    return {(m.gene, m.pos, m.hm) for m in apobecs}


APOBEC_MUTATION_MAP = apobec_mutation_map()


def apobec_2019_mutation_map():
    with open('/Users/philip/Dropbox/APOBECs/APBOECs20190809.csv') as fp:
        reader = csv.reader(fp)
        return {(g, int(p), a) for g, p, a in reader}


APOBEC_2019_MUTATION_MAP = apobec_2019_mutation_map()


def iter_sequences(dc):
    query = models.Isolate.make_query('HIV1', dc, 'all', [])
    query = query.options(db.selectinload('sequences'))
    query = query.options(db.joinedload('clinical_isolate'))
    isolates = dbutils.chunk_query(
        query, models.Isolate.id, chunksize=20000,
        on_progress=(lambda o, t:
                     print('{0}/{1} isolates...'.format(o, t), end='\r')),
        on_finish=(lambda t:
                   print('{0} isolates.                        '.format(t)))
    )
    for isolate in isolates:
        for seq in isolate.sequences:
            if seq.sequence_type != 'Sequence':
                continue
            if seq.naseq is None:
                continue
            yield isolate, seq


def compare(refseq, seq, gene, gene_size):
    naseq = seq.naseq
    first_aa = seq.first_aa
    last_aa = seq.last_aa
    naseq = ('...' * (first_aa - 1) +
             naseq +
             '...' * (gene_size - last_aa))
    hypermut_result = hypermut(refseq, naseq)
    hivdb_result = [0, 0, []]
    hivdb2019_result = [0, 0, []]
    for pos, aas in seq.aas:
        if len(aas) > 4:
            # ignore highly ambiguous sites
            continue
        refaa = translate_codon(refseq[pos * 3 - 3:pos * 3])
        if any((gene, pos, aa) in APOBEC_MUTATION_MAP
               for aa in aas):
            hivdb_result[0] += 1
            hivdb_result[2].append('{}{}{}'.format(refaa, pos, aas))
        if any((gene, pos, aa) in APOBEC_2019_MUTATION_MAP
               for aa in aas):
            hivdb2019_result[0] += 1
            hivdb2019_result[2].append('{}{}{}'.format(refaa, pos, aas))
        if '*' in aas:
            hivdb_result[1] += 1
            hivdb2019_result[1] += 1
        if any((gene, pos, aa) in ACTIVE_SITES for aa in aas):
            hivdb_result[1] += 1
            hivdb2019_result[1] += 1
    return hypermut_result, hivdb_result, hivdb2019_result


def compare_all():
    for geneidx, gene in enumerate(('PR', 'RT', 'IN')):
        gene_size = GENE_SIZES[geneidx]
        refseq = GENE_NACONS[geneidx]
        dc = DRUG_CLASSES[geneidx]
        for isolate, seq in iter_sequences(dc):
            hypermut_result, hivdb_result, hivdb2019_result = \
                compare(refseq, seq, gene, gene_size)
            hasclin = isolate.clinical_isolate is not None
            yield {
                'Gene': gene,
                'IsolateID': isolate.id,
                'PtID': isolate.patient_id,
                'IsolateDate': isolate.isolate_date,
                'AccessionID': seq.accession,
                'Source': (isolate.clinical_isolate.source
                           if hasclin else None),
                'SeqMethod': (isolate.clinical_isolate.seq_method
                              if hasclin else None),
                'Hypermut PositiveAPOBECs': hypermut_result[0],
                'Hypermut PotentialAPOBECs': hypermut_result[1],
                'Hypermut PositiveControls': hypermut_result[2],
                'Hypermut PotentialControls': hypermut_result[3],
                'Hypermut Ratio (A/B)/(C/D)': hypermut_result[4],
                'Adjusted Hypermut Ratio (+1)': hypermut_result[5],
                'Hypermut P-value (Fisher Exact)': hypermut_result[6],
                'HIVDB APOBEC Count': hivdb_result[0],
                'HIVDB 2019APOBEC Count': hivdb2019_result[0],
                'Stop Codon Count': hivdb_result[1],
                'HIVDB APOBECs': ', '.join(hivdb_result[2]),
                'HIVDB 2019APOBECs': ', '.join(hivdb2019_result[2])
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
    main()
