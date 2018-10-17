import re
import requests

URL_SCORED_MUTS = (
    'https://github.com/hivdb/sierra/raw/master/DrugResistance/src/main/'
    'resources/__cached_classes/edu.stanford.hivdb.drugresistance.database.'
    'MutationScores/mutScores.json'
)

URL_SCORED_COMBO_MUTS = (
    'https://github.com/hivdb/sierra/raw/master/DrugResistance/src/main/'
    'resources/__cached_classes/edu.stanford.hivdb.drugresistance.database.'
    'MutationComboScores/combinationScores.json'
)


def get_scored_mutations(includes_combo=False):
    resp = requests.get(URL_SCORED_MUTS)
    data = resp.json()
    result = {(one['gene'], one['pos'], one['aa'], True) for one in data}
    if includes_combo:
        resp = requests.get(URL_SCORED_COMBO_MUTS)
        data = resp.json()
        for one in data:
            gene = one['gene']
            muts = one['rule'].split('+')
            for mut in muts:
                pos, aas = re.match(r'^(\d+)(.+)$', mut).groups()
                for aa in aas:
                    result.add((gene, int(pos), aa, False))
    return result
