#! /usr/bin/env python

# Python 3.6+ is required to preserve dict key order

import csv
import click

import numpy as np
# import matplotlib.pyplot as plot
from sklearn.decomposition import PCA


@click.command()
@click.argument('input-file', type=click.File('r'))
@click.argument('output-file', type=click.File('w'), default='-')
def main(input_file, output_file):
    data = []
    reader = csv.reader(input_file, delimiter='\t')
    _, _, *feature_names, _ = next(reader)
    for pos, aa, *features, _ in reader:
        data.append([int(f) for f in features])
    arr = np.transpose(data)
    pca = PCA()
    arr = pca.fit_transform(arr)
    # plot.figure(figsize=(12, 12))
    writer = csv.writer(output_file, delimiter='\t')
    writer.writerow(
        ['Feature'] + ['PC{}'.format(i)
                       for i in range(1, 1 + len(feature_names))])
    for f, row in zip(feature_names, arr):
        writer.writerow([f, *row])

    # for (x, y), muts in grouped_muts.items():
    #     for mut in muts:
    #         writer.writerow([mut, x, y])
    #     muts = ','.join(muts)
    #     xoffset = .026 * len(muts)
    #     plot.annotate(muts, (x - xoffset, y + .05))
    # plot.xlim(-3, 3)
    # plot.ylim(-3, 3)
    # plot.xlabel('PC1')
    # plot.ylabel('PC2')
    # plot.savefig("plots/pca-drms.png", dpi=300)


if __name__ == '__main__':
    main()
