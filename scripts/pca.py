#! /usr/bin/env python

import csv
import click
import random

import numpy as np
import pandas as pd
import matplotlib.pyplot as plot
from sklearn.decomposition import PCA
from adjustText import adjust_text


def plot_pca_all(data, pca_arr):
    mutations = np.array(data.index)
    plot.figure(figsize=(40, 24))
    n = 1
    for i in range(5):
        for j in range(i + 1, 6):
            plot.subplot(3, 5, n)
            n += 1
            xs = pca_arr[:, i]
            ys = pca_arr[:, j]
            plot.scatter(xs, ys, s=1, c='white')
            texts = []
            for x, y, mut in zip(xs, ys, mutations):
                # adjust_text has a bug that text with exactly same x/y
                # sometime overlap together perfectly; by slightly tweaking
                # the coordinate can avoid the issue
                x += (random.random() - 0.5) / 50
                y += (random.random() - 0.5) / 50
                texts.append(plot.text(x, y, mut, fontsize=9))
            plot.xlabel('PC{}'.format(i + 1))
            plot.ylabel('PC{}'.format(j + 1))
            print(i, j, adjust_text(texts, lim=10, text_from_points=False))
    plot.savefig("plots/pca-drms.pdf")


def plot_pca(data, pca_arr):
    mutations = np.array(data.index)
    plot.figure(figsize=(8, 8))
    xs = pca_arr[:, 0]
    ys = pca_arr[:, 1]
    plot.scatter(xs, ys, s=1, c='white')
    texts = []
    for x, y, mut in zip(xs, ys, mutations):
        # adjust_text has a bug that text with exactly same x/y
        # sometime overlap together perfectly; by slightly tweaking
        # the coordinate can avoid the issue
        x += (random.random() - 0.5) / 50
        y += (random.random() - 0.5) / 50
        texts.append(plot.text(x, y, mut, fontsize=9))
    plot.xlabel('PC1')
    plot.ylabel('PC2')
    adjust_text(texts, lim=10, text_from_points=False)
    plot.savefig("plots/pca-drms-pc1-pc2.pdf")


@click.command()
@click.argument('input-file', type=click.File('r'))
@click.argument('output-file', type=click.File('w'), default='-')
@click.option('--plot', is_flag=True)
def main(input_file, output_file, plot):
    data = pd.read_csv(input_file, delimiter='\t', index_col=0, dtype={
        'Mut': str, 'Algs': int, 'Freq': int, 'Poly': int,
        'Selection': int, 'Susc': int, 'Position': int
    })
    arr = np.array(data)
    pca = PCA()
    arr = pca.fit_transform(arr)
    arr = np.array(pd.DataFrame(arr).mul(-1))
    writer = csv.writer(output_file, delimiter='\t')
    pcs = ['PC{}'.format(i + 1) for i in range(6)]
    writer.writerow(['', *pcs])
    writer.writerow(['Components'])
    for i, (feature, comp) in enumerate(
            zip(data.columns.values, pca.components_)):
        writer.writerow(['PC{}'.format(i + 1), *comp])
    writer.writerow([
        'Principal axes in feature space, representing the directions of '
        'maximum variance in the data. The components are sorted by '
        '"Explained Variance".'])
    writer.writerow([])
    writer.writerow(['Explained Variance',
                     *pca.explained_variance_])
    writer.writerow([
        'The amount of variance explained by each of the selected components.'
    ])
    writer.writerow([])
    writer.writerow(['Explained Variance Ratio',
                     *pca.explained_variance_ratio_])
    writer.writerow([
        'Percentage of variance explained by each of the selected components.'
    ])
    writer.writerow([])
    writer.writerow(['Singular Values', *pca.singular_values_])
    writer.writerow([
        'The singular values corresponding to each of the selected components.'
    ])
    writer.writerow([])
    writer.writerow(['Mean', *pca.mean_])
    writer.writerow([
        'Per-feature empirical mean, estimated from the training set.'
    ])
    writer.writerow([])
    writer.writerow(['Mutation'])
    mutations = np.array(data.index)
    for f, row in zip(mutations, arr):
        writer.writerow([f, *row])
    if plot:
        plot_pca(data, arr)


if __name__ == '__main__':
    main()
