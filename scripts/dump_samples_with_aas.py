#! /usr/bin/env python
import csv
import click
import pymysql.cursors

from codonutils import translate_codon

GENE_CHOICES = ('PR', 'RT', 'IN')
SQL_QUERY_INSTIS = """
SELECT PtID, IsolateID, NumIIs as NumDrugs FROM _INTotalRx
WHERE Unknown='No' ORDER BY IsolateID
"""
SQL_QUERY_INSTI_CODONS = """
SELECT IsolateID, cd.* FROM tblSequences s, _INNATriplets cd
WHERE s.SequenceID=cd.SequenceID AND IsolateID IN %s
ORDER BY IsolateID
"""


def iter_codons(conn, query, isolate_ids, every=1000):
    isolate_ids = list(isolate_ids)
    with conn.cursor() as cursor:
        offset = 0
        total = len(isolate_ids)
        while offset <= total:
            cursor.execute(
                query, (isolate_ids[offset:offset + every],))
            yield from cursor
            offset += every


@click.command()
@click.option('-h', '--host', type=str, default='127.0.0.1', help='MySQL host')
@click.option('-P', '--port', type=int, default=3308, help='MySQL port')
@click.option('-u', '--user', type=str, default='root', help='MySQL user')
@click.option('-p', '--password', type=str, help='MySQL password')
@click.option('-d', '--db', default='HIVDB2', help='MySQL database name')
@click.option('-o', '--out', type=click.File('w'),
              default='-', help='output target file')
@click.argument('genes', nargs=-1, required=True,
                type=click.Choice(GENE_CHOICES))
def dump_samples_with_aas(host, port, user, password, db, out, genes):
    if 'PR' in genes or 'RT' in genes:
        raise NotImplementedError
    conn = pymysql.connect(
        host=host, port=port, user=user, passwd=password,
        db=db, cursorclass=pymysql.cursors.DictCursor)
    with conn.cursor() as cursor:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['PtID', 'IsolateID', 'NumDrugs',
                         'Position', 'Codon', 'AA'])
        cursor.execute(SQL_QUERY_INSTIS)
        isolates = cursor.fetchall()
        isolates = {i['IsolateID']: i for i in isolates}
        isolate_ids = isolates.keys()
        for one in iter_codons(conn, SQL_QUERY_INSTI_CODONS, isolate_ids):
            isolate = isolates[one['IsolateID']]
            for pos in range(1, 289):
                codon = one['P{}'.format(pos)]
                codon = codon.replace('~', '-')
                aa = '.' if '.' in codon else translate_codon(codon)
                writer.writerow([
                    isolate['PtID'], isolate['IsolateID'],
                    isolate['NumDrugs'], pos, codon, aa
                ])


if __name__ == '__main__':
    dump_samples_with_aas()
