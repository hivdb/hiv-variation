#! /usr/bin/env python

# Python 3.6+ is required to preserve dict key order

import csv
import click

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plot
# from sklearn.decomposition import PCA


@click.command()
@click.argument('input-file', type=click.File('r'))
@click.argument('output-file', type=click.File('w'), default='-')
def main(input_file, output_file):
    data = []
    reader = csv.reader(input_file, delimiter='\t')
    next(reader)  # skip title row
    muts = []
    for pos, aa, *features, _ in reader:
        muts.append('{}{}'.format(pos, aa))
        data.append([int(f) for f in features])
    arr = np.array(data)
    linked = linkage(arr)
    plot.figure(figsize=(6, 15))
    dendrogram(linked,
               orientation='left',
               labels=muts,
               distance_sort='descending',
               show_leaf_counts=True)
    currentaxis = plot.gca()
    for spine in currentaxis.spines.values():
        spine.set_visible(False)
    currentaxis.get_xaxis().set_visible(False)

    plot.savefig("plots/hierarcical-cluster.png", dpi=300)
    print(linked)

    # pca = PCA()
    # arr = pca.fit_transform(arr)
    # writer = csv.writer(output_file, delimiter='\t')
    # writer.writerow(
    #     ['Feature'] + ['PC{}'.format(i)
    #                    for i in range(1, 1 + len(feature_names))])
    # for f, row in zip(feature_names, arr):
    #     writer.writerow([f, *row])


if __name__ == '__main__':
    main()
