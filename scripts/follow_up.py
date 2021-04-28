import json
from rank_compounds import rank, get_all_smiles

def get_reduced_smiles_bits(smiles_bits, test_target):
    
    reduced = {}

    for t in smiles_bits:
        if t != test_target:
            reduced[t] = smiles_bits[t]
    
    return reduced

smiles_bits = json.load(open('data/datafiles/smiles_bits.json', 'r'))


for test_target in smiles_bits:

    print(test_target)

    smiles_bits_reduced = get_reduced_smiles_bits(smiles_bits, test_target)
    paired_compounds = []

    for target in smiles_bits_reduced:
        for smiles1 in smiles_bits_reduced[target]:
            for smiles2 in smiles_bits_reduced[target]:
                sim = len(set(smiles_bits_reduced[target][smiles1]).intersection(set(smiles_bits_reduced[target][smiles2])))/len(set(smiles_bits_reduced[target][smiles1]).union(set(smiles_bits_reduced[target][smiles2])))
                if smiles1 != smiles2 and sim > 0.0:
                    if sorted([smiles1, smiles2]) not in paired_compounds:
                        paired_compounds.append(sorted([smiles1, smiles2]))


    all_smiles = get_all_smiles(smiles_bits_reduced)
    comps, x = rank(smiles_bits_reduced, all_smiles, 250)

    total_information = []
    for c in smiles_bits[test_target]:
        total_information += smiles_bits[test_target][c]
    total_information = len(set(total_information))
    target_fraction = get_fraction(comps, smiles_bits, test_target, total_information)

    print('hits recovered in first 250:', len([i for i in smiles_bits[test_target] if i in comps])/len(smiles_bits[test_target]))
    print('information recovered in first 250:', target_fraction)

    print(len(smiles_bits[test_target]), len([i for i in smiles_bits[test_target] if i not in comps]))
    pred = 0
    red = 0
    for comp in comps:
        if comp in smiles_bits[test_target]:
            for pair in paired_compounds:
                if comp in pair:
                    rest = [i for i in pair if i != comp][0]
                    if rest not in comps and rest in smiles_bits[test_target]:
                        pred += 1
                    else:
                        red += 1
    print(test_target, pred, red)

