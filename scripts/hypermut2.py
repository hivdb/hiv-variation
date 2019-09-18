#! /usr/bin/env python
import re
import sys
from scipy.stats import fisher_exact

inf = float('inf')

IUPAC_CODES = {
    'R': r'[RAG]',
    'Y': r'[YCT]',
    'B': r'[BYSKCGT]',
    'D': r'[DRWKAGT]',
    'H': r'[HYWMACT]',
    'V': r'[VRSMAGC]',
    'N': r'[NRYBDHVWSKMACGT]',
    'W': r'[WAT]',
    'S': r'[SCG]',
    'K': r'[KGT]',
    'M': r'[MAC]',
}

VALID_NA_PATTERN = re.compile(r'[RYBDHVWSKMACGT]')


def expand_iupac(pattern):
    result = []
    for char in pattern:
        if char in IUPAC_CODES:
            result.append(IUPAC_CODES[char])
        else:
            result.append(char)
    return ''.join(result)


DEFAULT_PATTERNS = {
    'hypermut_from': re.compile(expand_iupac(r'G(?=RD)')),
    'hypermut_to': re.compile(expand_iupac(r'A(?=RD)')),
    'control_from': re.compile(expand_iupac(r'G(?=YN|RC)')),
    'control_to': re.compile(expand_iupac(r'A(?=YN|RC)')),
}


def fasta_reader(filename):
    with open(filename) as fp:
        header = None
        seq = []
        for line in fp:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if seq:
                    yield header, ''.join(seq).upper()
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            yield header, ''.join(seq).upper()


CACHE = {}


def cached_find_sites(seq, pattern, site_range):
    cachekey = (seq, pattern.pattern)
    if cachekey not in CACHE:
        CACHE[cachekey] = find_sites(seq, pattern, range(len(seq)))
    available_sites = set(CACHE[cachekey])
    site_range = set(site_range)
    return sorted(available_sites & site_range)


def find_sites(seq, pattern, site_range):
    sites = []
    site_range = set(site_range)
    min_site = min(site_range)
    max_site = max(site_range)
    if not site_range:
        return sites
    for match in pattern.finditer(seq[min_site:]):
        offset = min_site + match.span()[0]
        if offset not in site_range:
            continue
        if offset > max_site:
            break
        sites.append(offset)
    return sites


def get_comparable_sites(refseq, naseq):
    sites = []
    for offset, (ref, na) in enumerate(zip(refseq, naseq)):
        if not VALID_NA_PATTERN.match(ref):
            continue
        if not VALID_NA_PATTERN.match(na):
            continue
        sites.append(offset)
    return sites


def hypermut(refseq, naseq, patterns=DEFAULT_PATTERNS):
    comparable_sites = get_comparable_sites(refseq, naseq)
    potential_muts = cached_find_sites(
        refseq, patterns['hypermut_from'], comparable_sites)
    potential_ctrls = cached_find_sites(
        refseq, patterns['control_from'], comparable_sites)
    matched_muts = find_sites(
        naseq, patterns['hypermut_to'], potential_muts)
    matched_ctrls = find_sites(
        naseq, patterns['control_to'], potential_ctrls)
    num_potential_muts = len(potential_muts)
    num_matched_muts = len(matched_muts)
    num_potential_ctrls = len(potential_ctrls)
    num_matched_ctrls = len(matched_ctrls)
    try:
        oddsratio = (
            (num_matched_muts / num_potential_muts) /
            (num_matched_ctrls / num_potential_ctrls)
        )
    except ZeroDivisionError:
        oddsratio = inf
    adjusted_oddsratio = (
        ((num_matched_muts + 1) / (num_potential_muts + 1)) /
        ((num_matched_ctrls + 1) / (num_potential_ctrls + 1))
    )
    _, p = fisher_exact([
        [num_matched_muts, num_potential_muts - num_matched_muts],
        [num_matched_ctrls, num_potential_ctrls - num_matched_ctrls]
    ], 'greater')
    return (
        num_matched_muts,
        num_potential_muts,
        num_matched_ctrls,
        num_potential_ctrls,
        oddsratio,
        adjusted_oddsratio,
        p)


def main():
    if len(sys.argv) != 2:
        print('Usage: {} <FASTA_FILE>'.format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    fasta_filename = sys.argv[1]
    sequences = list(fasta_reader(fasta_filename))
    _, refseq = sequences.pop(0)
    for _, naseq in sequences:
        print(hypermut(refseq, naseq, DEFAULT_PATTERNS))


if __name__ == '__main__':
    main()
