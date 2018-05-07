import os
import json
from collections import defaultdict

MUTSCORES_FILENAME = '../data/mut_scores.json'
MUTSCORES_FILENAME = os.path.abspath(os.path.join(
    os.path.dirname(__file__), MUTSCORES_FILENAME
))


def build_drmlookup(min_drm_score=0):
    with open(MUTSCORES_FILENAME) as fp:
        data = json.load(fp)
    result = defaultdict(set)
    for one in data:
        if one['score'] < min_drm_score:
            continue
        result[one['gene']].add((one['pos'], one['aa']))
    return result
