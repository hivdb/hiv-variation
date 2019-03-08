#! /usr/bin/env python

import click
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
from sklearn.manifold import MDS
from adjustText import adjust_text


@click.command()
@click.argument('input-file', type=click.File('r'))
@click.argument('output-file', type=click.File('w'), default='-')
def main(input_file, output_file):
    data = pd.read_csv(input_file, delimiter='\t', index_col=0, dtype={
        'Mut': str, 'Algs': int, 'Freq': int, 'Poly': int,
        'Selection': int, 'Susc': int, 'Position': int
    })
    arr = np.array(data)
    mds = MDS(random_state=4)
    pos = mds.fit_transform(arr)
    plot.figure(figsize=(8, 8))
    texts = []
    plot.scatter(pos[:, 0], pos[:, 1], c='white')
    for (x, y), mut in zip(pos, data.index):
        x += (random.random() - .5) / 50
        y += (random.random() - .5) / 50
        texts.append(plot.text(x, y, mut, fontsize=9))
    print(pos)
    print(adjust_text(texts, lim=20, text_from_points=False))
    plot.savefig("plots/mds-drms.pdf")


if __name__ == '__main__':
    main()
