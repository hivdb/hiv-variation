#! /usr/bin/env python
import csv
from hivdbql import app
from itertools import groupby

models = app.models
Isolate = models.Isolate
RefLink = models.RefLink

AA_RANGES = {
    'PR': (1, 100),
    'RT': (1, 241),
    'IN': (1, 271)
}


def create_gene_table(gene, drugclass, target_fp):
    query = Isolate.make_query(species='HIV2',
                               drugclass=drugclass,
                               rx_type='all',
                               criteria=['NO_QA_ISSUES'])
    query = query.filter(~Isolate.references.any(RefLink.reference_id == 2350))
    query = query.order_by(Isolate.patient_id, Isolate.isolate_date)
    patients = groupby(query, lambda i: i.patient_id)
    writer = csv.writer(target_fp, delimiter='\t')
    writer.writerow(
        ['PtID', 'IsolateID', 'IsolateDate', 'Rx', 'Genotype'] +
        ['P{}'.format(i) for i in range(*AA_RANGES[gene])])
    for ptid, isolates in patients:
        for isolate in isolates:
            if not isolate:
                continue
            aas = dict(isolate.get_or_create_consensus().aas)
            writer.writerow(
                [ptid, isolate.id,
                 isolate.isolate_date.strftime('%Y-%m-%d'),
                 ', '.join(isolate.drugs) or 'None',
                 isolate.subtype] +
                [aas.get(i, '.') for i in range(*AA_RANGES[gene])])


if __name__ == '__main__':
    with app.app_context():
        with open('/tmp/PRTable.txt', 'w') as fp:
            create_gene_table('PR', 'PI', fp)
        with open('/tmp/RTTable.txt', 'w') as fp:
            create_gene_table('RT', 'NRTI', fp)
        with open('/tmp/INTable.txt', 'w') as fp:
            create_gene_table('IN', 'INSTI', fp)
