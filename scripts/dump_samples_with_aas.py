#! /usr/bin/env python
import os
import csv
import click
import pymysql.cursors

from codonutils import translate_codon

GENE_CHOICES = ('PR', 'RT', 'IN')
GENE_RANGES = {
    'PR': list(range(1, 100)),
    'RT': list(range(1, 561)),
    'IN': list(range(1, 289)),
}
SQL_QUERY_INSTIS = """
SELECT i.PtID, i.IsolateID, i.NumIIs as NumDrugs, s.Subtype
FROM _INTotalRx i JOIN tblSubtypes s ON s.IsolateID=i.IsolateID
WHERE Unknown='No' ORDER BY IsolateID
"""
SQL_QUERY_SEQS = """
SELECT IsolateID, SequenceID, Firstaa, Lastaa, NASeq, AASeq
FROM tblSequences WHERE IsolateID IN %s ORDER BY IsolateID
"""


def iter_seqs(conn, query, isolate_ids, every=1000):
    isolate_ids = list(isolate_ids)
    with conn.cursor() as cursor:
        offset = 0
        total = len(isolate_ids)
        while offset <= total:
            cursor.execute(
                query, (isolate_ids[offset:offset + every],))
            yield from cursor
            offset += every


def iter_codons(seq, gene_range):
    firstaa = seq['Firstaa']
    lastaa = seq['Lastaa']
    naseq = seq['NASeq']
    aaseq = seq['AASeq']
    aaseq_only = not naseq
    for pos in gene_range:
        if pos < firstaa or pos > lastaa:
            yield pos, '...', '.'
            continue
        aaidx0 = pos - firstaa
        if aaseq_only:
            yield pos, '...', aaseq[aaidx0]
        else:
            naidx0 = aaidx0 * 3
            codon = naseq[naidx0:naidx0 + 3]
            codon = codon.replace('~', '-')
            aa = '.' if '.' in codon else translate_codon(codon)
            yield pos, codon, aa


@click.command()
@click.option('-h', '--host', type=str, default='127.0.0.1', help='MySQL host')
@click.option('-P', '--port', type=int, default=3308, help='MySQL port')
@click.option('-u', '--user', type=str, default='root', help='MySQL user')
@click.option('-p', '--password', type=str, help='MySQL password')
@click.option('-d', '--db', default='HIVDB2', help='MySQL database name')
@click.option('-o', '--outdir', type=str, default='../local',
              help='output target directory')
@click.argument('genes', nargs=-1, required=True,
                type=click.Choice(GENE_CHOICES))
def dump_samples_with_aas(host, port, user, password, db, outdir, genes):
    if 'PR' in genes or 'RT' in genes:
        raise NotImplementedError
    conn = pymysql.connect(
        host=host, port=port, user=user, passwd=password,
        db=db, cursorclass=pymysql.cursors.DictCursor)
    outdir = os.path.abspath(os.path.join(
        os.path.dirname(__file__), outdir
    ))
    os.makedirs(outdir, exist_ok=True)
    file_isolates = os.path.join(outdir, 'in_isolates.txt')
    file_codons = os.path.join(outdir, 'in_codons.txt')
    with open(file_isolates, 'w') as file_isolates, \
            open(file_codons, 'w') as file_codons, \
            conn.cursor() as cursor:
        writer_isolates = csv.DictWriter(
            file_isolates,
            ['PtID', 'IsolateID', 'NumDrugs', 'Subtype'],
            delimiter='\t')
        writer_codons = csv.writer(file_codons, delimiter='\t')
        cursor.execute(SQL_QUERY_INSTIS)
        isolates = cursor.fetchall()
        writer_isolates.writeheader()
        writer_isolates.writerows(isolates)
        isolates = {i['IsolateID']: i for i in isolates}
        isolate_ids = isolates.keys()
        writer_codons.writerow(['PtID', 'IsolateID', 'NumDrugs',
                                'Subtype', 'Position', 'Codon', 'AA'])
        gene_range = GENE_RANGES['IN']
        for one in iter_seqs(conn, SQL_QUERY_SEQS, isolate_ids):
            isolate = isolates[one['IsolateID']]
            for pos, codon, aa in iter_codons(one, gene_range):
                writer_codons.writerow([
                    isolate['PtID'], isolate['IsolateID'],
                    isolate['NumDrugs'], isolate['Subtype'],
                    pos, codon, aa
                ])


if __name__ == '__main__':
    dump_samples_with_aas()
