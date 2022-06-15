#! /usr/bin/env python
"""A Python implementation of Hypermut2

See https://www.hiv.lanl.gov/cgi-bin/HYPERMUT/hypermut.cgi
"""

import re
import csv
import sys
from scipy.stats import fisher_exact  # type: ignore
from typing import TextIO, Literal, Dict, Tuple, Iterable, Set, List, Optional

import click  # type: ignore


inf: float = float('inf')

IUPAC_CODES: Dict[str, str] = {
    'R': r'[AG]',
    'Y': r'[CT]',
    'B': r'[^A]',
    'D': r'[^C]',
    'H': r'[^G]',
    'V': r'[^T]',
    'N': r'.',
    'W': r'[AT]',
    'S': r'[CG]',
    'K': r'[GT]',
    'M': r'[AC]',
}


def expand_iupac(pattern: str) -> str:
    result = []
    for char in pattern:
        if char in IUPAC_CODES:
            result.append(IUPAC_CODES[char])
        else:
            result.append(char)
    return ''.join(result)


def compile_mutation(mut: str) -> re.Pattern:
    return re.compile(expand_iupac(mut))


def compile_context(
    upstream: str,
    downstream: str
) -> re.Pattern:
    return re.compile(expand_iupac(
        '(?={}).(?={})'.format(
            upstream,
            downstream
        )
    ))


def fasta_reader(
    fasta_file: TextIO
) -> Iterable[Tuple[Optional[str], str]]:
    header: Optional[str] = None
    seq: List[str] = []
    for line in fasta_file:
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


CACHE: Dict[Tuple[str, str], Set[int]] = {}


def cached_find_sites(
    seq: str,
    pattern: re.Pattern
) -> Set[int]:
    cachekey: Tuple[str, str] = (seq, pattern.pattern)
    if cachekey not in CACHE:
        CACHE[cachekey] = find_sites(seq, pattern)
    return CACHE[cachekey]


def find_sites(
    seq: str,
    pattern: re.Pattern
) -> Set[int]:
    sites: Set[int] = set()
    for mat in pattern.finditer(seq):
        offset = mat.span()[0]
        sites.add(offset)
    return sites


def hypermut(
    refseq: str,
    naseq: str,
    patterns: Dict[str, re.Pattern],
    enforce_context: Literal['ref', 'both', 'query']
) -> Tuple[int, int, int, int, float, float, str, str, str, str]:

    potential_muts: Set[int] = set(range(min(len(refseq), len(naseq))))
    potential_ctrls: Set[int] = set(range(min(len(refseq), len(naseq))))

    if enforce_context in ('ref', 'both'):
        potential_muts &= cached_find_sites(
            refseq, patterns['hypermut_context'])
        potential_ctrls &= cached_find_sites(
            refseq, patterns['control_context'])
    if enforce_context in ('query', 'both'):
        potential_muts &= find_sites(
            naseq, patterns['hypermut_context'])
        potential_ctrls &= find_sites(
            naseq, patterns['control_context'])
    potential_muts &= cached_find_sites(
        refseq, patterns['hypermut_from'])
    potential_ctrls &= cached_find_sites(
        refseq, patterns['control_from'])

    matched_muts = find_sites(
        naseq, patterns['hypermut_to']) & potential_muts
    matched_ctrls = find_sites(
        naseq, patterns['control_to']) & potential_ctrls

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
    # adjusted_oddsratio = (
    #     ((num_matched_muts + 1) / (num_potential_muts + 1)) /
    #     ((num_matched_ctrls + 1) / (num_potential_ctrls + 1))
    # )
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
        # adjusted_oddsratio,
        p,
        ','.join([str(pos0 + 1) for pos0 in sorted(matched_muts)]),
        ','.join([str(pos0 + 1) for pos0 in sorted(potential_muts)]),
        ','.join([str(pos0 + 1) for pos0 in sorted(matched_ctrls)]),
        ','.join([str(pos0 + 1) for pos0 in sorted(potential_ctrls)])
    )


@click.command()
@click.argument('fasta_file', type=click.File('r'))
@click.option(
    '-e', '--enforce-context',
    type=click.Choice(['ref', 'both', 'query']),
    default='query',
    show_default=True,
    help=(
        'Should the downstream context be enforced to the reference, '
        'query, or both sequences.'
    )
)
@click.option(
    '-U', '--hypermut-upstream',
    default='',
    show_default=True,
    help='The upstream context of hypermut pattern.'
)
@click.option(
    '-u', '--control-upstream',
    default='',
    show_default=True,
    help='The upstream context of control pattern.'
)
@click.option(
    '-F', '--hypermut-from',
    default='G',
    show_default=True,
    help='The mutation from of hypermut pattern.'
)
@click.option(
    '-T', '--hypermut-to',
    default='A',
    show_default=True,
    help='The mutation to of hypermut pattern.'
)
@click.option(
    '-f', '--control-from',
    default='G',
    show_default=True,
    help='The mutation from of control pattern.'
)
@click.option(
    '-t', '--control-to',
    default='A',
    show_default=True,
    help='The mutation to of control pattern.'
)
@click.option(
    '-D', '--hypermut-downstream',
    default='RD',
    show_default=True,
    help='The downstream context of hypermut pattern.'
)
@click.option(
    '-d', '--control-downstream',
    default='YN|RC',
    show_default=True,
    help='The downstream context of control pattern.'
)
def main(
    fasta_file: TextIO,
    enforce_context: Literal['ref', 'both', 'query'],
    hypermut_upstream: str,
    control_upstream: str,
    hypermut_from: str,
    hypermut_to: str,
    control_from: str,
    control_to: str,
    hypermut_downstream: str,
    control_downstream: str
) -> None:
    sequences = list(fasta_reader(fasta_file))
    _, refseq = sequences.pop(0)

    patterns = {
        'hypermut_context': compile_context(
            hypermut_upstream,
            hypermut_downstream
        ),
        'control_context': compile_context(
            control_upstream,
            control_downstream
        ),
        'hypermut_from': compile_mutation(hypermut_from),
        'hypermut_to': compile_mutation(hypermut_to),
        'control_from': compile_mutation(control_from),
        'control_to': compile_mutation(control_to)
    }

    writer = csv.writer(sys.stdout, delimiter='\t')
    writer.writerow([
        'sequence',
        'num_matched_muts',
        'num_potential_muts',
        'num_matched_ctrls',
        'num_potential_ctrls',
        'oddsratio',
        # 'adjusted_oddsratio',
        'p',
        'matched_muts',
        'potential_muts',
        'matched_ctrls',
        'potential_ctrls'
    ])
    for name, naseq in sequences:
        writer.writerow([
            name,
            *hypermut(refseq, naseq, patterns, enforce_context)
        ])


if __name__ == '__main__':
    main()
