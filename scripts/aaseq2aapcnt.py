#! /usr/bin/env python
from __future__ import print_function

import sys
import json
from collections import Counter


def fasta_reader(filename):
    with open(filename) as fp:
        header = None
        seq = []
        for line in fp:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if seq:
                    yield header, ''.join(seq)
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            yield header, ''.join(seq)


def stat(gene, refheader, infile):
    sequences = list(fasta_reader(infile))
    counter = Counter()
    for header, seq in sequences:
        if header.endswith(refheader):
            refseq = seq
            break
    for header, seq in sequences:
        if header.endswith(refheader):
            continue
        pos = 0
        for r, a in zip(refseq, seq):
            if r != '-':
                pos += 1
                if a == '-':
                    continue
                counter[(pos, a)] += 1
                counter[(pos, '.')] += 1
    for (pos, aa), count in sorted(counter.items()):
        if aa == '.':
            continue
        total = counter[(pos, '.')]
        yield {
            'gene': gene,
            'position': pos,
            'subtype': 'All',
            'aa': aa,
            'rx_type': 'naive',
            'percent': count / total,
            'count': count,
            'total': total
        }


def main():
    if len(sys.argv) != 5:
        print('Usage: {} <GENE> <REF_HEADER> <INPUT FASTA> <OUTPUT JSON>'
              .format(sys.argv[0]), file=sys.stderr)
        exit(1)
    with open(sys.argv[4], 'w') as outfile:
        data = list(stat(sys.argv[1], sys.argv[2], sys.argv[3]))
        json.dump(data, outfile, indent=2)


if __name__ == '__main__':
    main()
