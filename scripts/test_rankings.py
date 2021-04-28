import os
import json
from rank_compounds import get_next_best_compound, get_all_smiles, rank
import matplotlib.pyplot as plt
import random


def get_reduced_smiles_bits(smiles_bits, test_target):
    
    reduced = {}

    for t in smiles_bits:
        if t != test_target:
            reduced[t] = smiles_bits[t]
    
    return reduced


def get_fraction(comps, smiles_bits, target, total_information):
    target_fraction = [0]
    current_info = []

    for comp in comps:

        if comp in smiles_bits[test_target]:
            current_info += smiles_bits[test_target][comp]

        current_info = list(set(current_info))
        target_fraction.append(len(current_info)/total_information)

    return target_fraction




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
random_fractions = []
comps_lists = []

for test_target in smiles_bits:

    smiles_bits_reduced = get_reduced_smiles_bits(smiles_bits, test_target)
    all_smiles = get_all_smiles(smiles_bits_reduced)
    comps, x = rank(smiles_bits_reduced, all_smiles, 250)
    comps_lists += comps


    total_information = []
    for c in smiles_bits[test_target]:
        total_information += smiles_bits[test_target][c]
    total_information = len(set(total_information))

    target_fraction = get_fraction(comps, smiles_bits, test_target, total_information)
    fractions.append(target_fraction)

    add_for_DSiP = 1007 - len(all_smiles)
    all_smiles += ['C'] * add_for_DSiP
    random.shuffle(all_smiles)
    random_fraction = get_fraction(all_smiles[:250], smiles_bits, test_target, total_information)
    random_fractions.append(random_fraction)


    print(test_target, len(comps), len(smiles_bits[test_target]), target_fraction[-1])
    plt.plot(range(len(target_fraction)), target_fraction, label=test_target)

    plt.xlabel('library size')
    plt.ylabel('information recovered')
    plt.ylim(0,1)
    plt.tight_layout()
    plt.savefig('figures/fraction_test.png')
    

mean_fractions = []
fraction_lengths = [len(i) for i in fractions]
for i in range(min(fraction_lengths)):
    mean_fractions.append(sum([x[i] for x in fractions])/len(fractions))

mean_random_fractions = []
random_fraction_lengths = [len(i) for i in random_fractions]
for i in range(min(random_fraction_lengths)):
    mean_random_fractions.append(sum([x[i] for x in random_fractions])/len(random_fractions))


plt.subplot(122)
plt.plot(range(len(mean_fractions)), mean_fractions, label='using ranked compounds')
plt.plot(range(len(mean_random_fractions)), mean_random_fractions, label='picking random compounds')
plt.ylim(0, 1)
plt.legend()
plt.xlabel('library size')
plt.ylabel('mean information recovered across all new targets')
plt.tight_layout()
plt.savefig('figures/fraction_test.png')


all_smiles = get_all_smiles(smiles_bits)
counts = [comps_lists.count(i) for i in all_smiles]
plt.close()
plt.figure(figsize=(8, 5))
plt.bar(list(set(counts)), [counts.count(i) for i in set(counts)])
plt.xlabel('number of times fragments in top 250')
plt.tight_layout()
plt.savefig('figures/compound_counts.png')


improvement = []
for i in range(len(fractions)):
    try:
        improvement.append((fractions[i][250]-random_fractions[i][250])/random_fractions[i][250])
    except:
        improvement.append(0)


plt.close()
plt.figure(figsize=(8, 5))
plt.barh([t for t in smiles_bits], improvement)
plt.xlabel('factor of information improvement in library size 250')
plt.ylabel('target')
plt.xlim(-2, 10)
plt.tight_layout()
plt.savefig('figures/improvement_by_target.png')