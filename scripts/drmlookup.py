import os
import json
from collections import defaultdict

from commons import get_scored_mutations

MUTSCORES_FILENAME = '../data/mut_scores.json'
MUTSCORES_FILENAME = os.path.abspath(os.path.join(
    os.path.dirname(__file__), MUTSCORES_FILENAME
))


def build_drmlookup_old(min_drm_score=0):
    with open(MUTSCORES_FILENAME) as fp:
        data = json.load(fp)
    result = defaultdict(set)
    for one in data:
        if one['score'] < min_drm_score:
            continue
        result[one['gene']].add((one['pos'], one['aa']))
    return result


def build_drmlookup(min_drm_score=0):
    drms = get_scored_mutations(True)
    result = defaultdict(set)
    for gene, pos, aa, is_major in drms:
        result[gene].add((pos, aa, is_major))
        result[gene].add((pos, aa))
    return result
