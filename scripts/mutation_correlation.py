#! /usr/bin/env python

import re
import csv

import click
import numpy as np
from scipy.stats import spearmanr
from hivdbql import app
from hivdbql.utils import dbutils

db = app.db
models = app.models

MUTATION_PATTERN = re.compile(r'^[A-Z]?(\d+)([A-Z*_-]+)$')


def read_mutations(fp):
    mutations = set()
    for line in fp:
        line = line.strip()
        match = MUTATION_PATTERN.match(line)
        if line and match:
            pos, aa = match.groups()
            mutations.add((int(pos), aa))
    orderedmuts = sorted(mutations)
    return {m: orderedmuts.index(m) for m in orderedmuts}


def calc_spearman(both, m0only, m1only, none):
    dataset = (
        [(1, 1)] * both + [(1, 0)] * m0only +
        [(0, 1)] * m1only + [(0, 0)] * none
    )
    return spearmanr(dataset)


@click.command()
@click.argument('input_mutations_file', type=click.File('r'))
@click.argument('output_file', type=click.File('w'))
def mutation_corellation(input_mutations_file, output_file):
    mutations = read_mutations(input_mutations_file)
    mutationitems = sorted(mutations.items(), key=lambda i: i[1])
    nummuts = len(mutations)
    writer = csv.writer(output_file)

    matrix = np.zeros([nummuts, nummuts, 0b100], dtype=np.int64)
    writer.writerow(['MutX', 'MutY', '#XY', '#X',
                     '#Y', '#Null', 'Rho', 'P'])
    # query = models.Isolate.make_query(
    #     'HIV1', 'INSTI', 'all', ['NO_CLONES',
    #                              'NO_QA_ISSUES',
    #                              'PUBLISHED_ONLY'])
    query = (
        models.Patient.query
        .filter(models.Patient.isolates.any(db.and_(
            *models.Isolate.make_criteria(
                # TODO: should be able to defined by command line
                'HIV1', 'INSTI', 'all', ['NO_CLONES',
                                         'NO_QA_ISSUES',
                                         'PUBLISHED_ONLY'])
        )))
        .options(db.selectinload(models.Patient.isolates)
                 .selectinload(models.Isolate.sequences)
                 .selectinload(models.Sequence.insertions))
        .options(db.selectinload(models.Patient.isolates)
                 .selectinload(models.Isolate.sequences)
                 .selectinload(models.Sequence.mixtures))
    )

    patients = dbutils.chunk_query(
        query, models.Patient.id, chunksize=500,
        on_progress=(lambda o, t:
                     print('{0}/{1} patients...'.format(o, t), end='\r')),
        on_finish=(lambda t:
                   print('{0} patients.              '.format(t)))
    )
    patcount = 0
    seqcount = 0
    for patient in patients:
        patmatrix = np.zeros_like(matrix)
        patflag = False
        for isolate in patient.isolates:
            # TODO: should be able to defined by command line
            if isolate.gene != 'IN':
                continue
            seq = isolate.get_or_create_consensus()
            first_aa = seq.first_aa
            last_aa = seq.last_aa
            seqmuts = {m for m in seq.aas if m in mutations}
            if not seqmuts:
                continue
            seqcount += 1
            patflag = True
            for m0, m0idx in mutationitems:
                if m0[0] < first_aa or m0[0] > last_aa:
                    # disqualified because of out of range
                    continue
                for m1, m1idx in mutationitems[m0idx + 1:]:
                    if m1[0] < first_aa or m1[0] > last_aa:
                        # disqualified because of out of range
                        continue
                    hasm0 = m0 in seqmuts
                    hasm1 = m1 in seqmuts
                    if hasm0 and hasm1:
                        # contains both
                        patmatrix[m0idx, m1idx, 0b11] = 1
                    elif hasm0 and not hasm1:
                        # contains m0
                        patmatrix[m0idx, m1idx, 0b10] = 1
                    elif not hasm0 and hasm1:
                        # contains m1
                        patmatrix[m0idx, m1idx, 0b01] = 1
                    else:  # elif not hasm0 and not hasm1:
                        # contains none
                        patmatrix[m0idx, m1idx, 0b00] = 1
        matrix += patmatrix
        patcount += patflag
    print('{} patients ({} sequences) have at least one given mutation.'
          .format(patcount, seqcount))
    for m0, m0idx in mutationitems:
        for m1, m1idx in mutationitems[m0idx + 1:]:
            both = matrix[m0idx, m1idx, 0b11]
            m0only = matrix[m0idx, m1idx, 0b10]
            m1only = matrix[m0idx, m1idx, 0b01]
            none = matrix[m0idx, m1idx, 0b00]
            if both != 0 or m0only * m1only != 0:
                rho, p = calc_spearman(both, m0only, m1only, none)
            else:
                rho = p = ''
            writer.writerow([
                '{}{}'.format(*m0),
                '{}{}'.format(*m1),
                both, m0only, m1only, none, rho, p
            ])


if __name__ == '__main__':
    with app.app_context():
        mutation_corellation()
