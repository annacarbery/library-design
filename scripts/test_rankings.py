import os
import json
from rank_compounds import get_next_best_compound, get_all_smiles, rank
import matplotlib.pyplot as plt

def get_reduced_smiles_bits(smiles_bits, test_target):
    
    reduced = {}

    for t in smiles_bits:
        if t != test_target:
            reduced[t] = smiles_bits[t]
    
    return reduced



smiles_bits = json.load(open('data/datafiles/smiles_bits.json', 'r'))

smol = []
for t in smiles_bits:
    if len(smiles_bits[t]) < 10:
        smol.append(t)

for t in smol:
    del smiles_bits[t]
# DSiP_smiles = all_smiles + ['C'] * 450
# random.shuffle(DSiP_smiles)

plt.figure(figsize=(12, 5))
plt.subplot(121)

fractions = []

for test_target in smiles_bits:

    smiles_bits_reduced = get_reduced_smiles_bits(smiles_bits, test_target)
    all_smiles = get_all_smiles(smiles_bits_reduced)
    comps, x = rank(smiles_bits_reduced, all_smiles, len(all_smiles))


    total_information = []
    for c in smiles_bits[test_target]:
        total_information += smiles_bits[test_target][c]
    total_information = len(set(total_information))

    target_fraction = [0]
    current_info = []

    for comp in comps:

        if comp in smiles_bits[test_target]:
            current_info += smiles_bits[test_target][comp]

        current_info = list(set(current_info))
        target_fraction.append(len(current_info)/total_information)

    fractions.append(target_fraction)

    print(test_target, len(comps), len(smiles_bits[test_target]), target_fraction[-1])
    plt.plot(range(len(target_fraction)), target_fraction, label=test_target)

    plt.xlabel('library size')
    plt.ylabel('information recovered')
    plt.tight_layout()
    plt.savefig('figures/fraction_test.png')
    

mean_fractions = []
fraction_lengths = [len(i) for i in fractions]
for i in range(min(fraction_lengths)):
    mean_fractions.append(sum([x[i] for x in fractions])/len(fractions))


plt.subplot(122)
plt.plot(range(len(mean_fractions)), mean_fractions)
plt.ylim(0, 1)
plt.xlabel('library size')
plt.ylabel('information recovered')
plt.tight_layout()
plt.savefig('figures/fraction_test.png')

