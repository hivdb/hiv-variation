#! /usr/bin/env python

import sys
import csv
import json

RX_TYPE = 'naive'
SUBTYPE = 'All'


def main():
    """Accept AAPcnt JSON file and output a csv table of coverage"""
    if len(sys.argv) != 3:
        print('Usage: {} <AAPCNT_JSON> <OUTPUT_CSV>'
              .format(sys.argv[0]), file=sys.stderr)
        exit(1)
    jsonfp, out = sys.argv[1:]
    with open(jsonfp) as jsonfp, open(out, 'w') as out:
        data = json.load(jsonfp)
        covs = {
            'PR': ['Coverage'] + [0] * 99,
            'RT': ['Coverage'] + [0] * 560,
            'IN': ['Coverage'] + [0] * 288
        }
        for row in data:
            if row['subtype'] == SUBTYPE and row['rx_type'] == RX_TYPE:
                pos0 = row['position']
                gene = row['gene']
                cov = row['total']
                covs[gene][pos0] = cov
        writer = csv.writer(out)
        writer.writerow(['PR'] + list(range(1, 100)))
        writer.writerow(covs['PR'])
        writer.writerow([])
        writer.writerow([])
        writer.writerow(['RT'] + list(range(1, 561)))
        writer.writerow(covs['RT'])
        writer.writerow([])
        writer.writerow([])
        writer.writerow(['IN'] + list(range(1, 289)))
        writer.writerow(covs['IN'])


if __name__ == '__main__':
    main()
