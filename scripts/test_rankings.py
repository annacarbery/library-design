import os
import json
from rank_compounds import get_next_best_compound, get_all_smiles, rank
import matplotlib.pyplot as plt
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover


def get_reduced_smiles_bits(smiles_bits, test_target):
    
    # get smiles bits without the target being tested 

    reduced = {}

    for t in smiles_bits:
        if t != test_target:
            reduced[t] = smiles_bits[t]
    
    return reduced


def get_fraction(comps, smiles_bits, target, total_information):

    # get the fraction of information recovered using a subset of DSiP compounds
    # returns list of fractional information

    target_fraction = [0]
    current_info = []

    for comp in comps:

        if comp in smiles_bits[test_target]:
            current_info += smiles_bits[test_target][comp]

        current_info = list(set(current_info))
        target_fraction.append(len(current_info)/total_information)

    return target_fraction


def ignore_targets(smiles_bits, target_tests, remove_less_than):

    smol = []
    for t in smiles_bits:

        if t in target_tests and target_tests[t] < remove_less_than:
            smol.append(t)
        elif t not in target_tests:
            smol.append(t)

    for t in smol:
        del smiles_bits[t]

    return smiles_bits


def get_total_information(smiles_bits, target):

    total_information = []
    for c in smiles_bits[target]:
        total_information += smiles_bits[target][c]
    total_information = len(set(total_information))
    return total_information


def multiple_random_runs(all_smiles, smiles_bits, test_target, total_information, lib_size=150, num_runs=1000):

    random_fractions_target = [0] * lib_size
    for i in range(num_runs):
        random.shuffle(all_smiles)
        random_run = get_fraction(all_smiles[:lib_size], smiles_bits, test_target, total_information)
        for l in range(lib_size):
            random_fractions_target[l] += random_run[l]
    random_fractions_target = [f/num_runs for f in random_fractions_target]
    random_fractions_target += [0]
    random_fractions_target = sorted(random_fractions_target)

    return random_fractions_target


def plot_target_fraction(target, fraction):

    plt.plot(range(len(fraction)), fraction, label=target)
    plt.xlabel('library size')
    plt.ylabel('information recovered')
    plt.ylim(0,1)
    plt.tight_layout()
    plt.savefig('figures/fraction_test.png')


def target_improvement(target_fraction, random_fractions_target):

    assert len(target_fraction) == len(random_fractions_target)

    if random_fractions_target[-1] > 0:
        imp = (target_fraction[-1]-random_fractions_target[-1])/random_fractions_target[-1]
    else:
        imp = 0
    
    return imp


def mean_across_targets(fraction):

    mean_fractions = []
    fraction_lengths = [len(i) for i in fractions]
    for i in range(min(fraction_lengths)):
        mean_fractions.append(sum([x[i] for x in fractions])/len(fractions))
    
    return mean_fractions


smiles_bits = json.load(open('data/datafiles/smiles_bits.json', 'r'))
frequent_comps = json.load(open('data/datafiles/frequently_tested_compounds.json', 'r'))
target_tests = json.load(open('data/datafiles/target_test_numbers.json', 'r'))
lib_size = 150

smiles_bits = ignore_targets(smiles_bits, target_tests, 250)

plt.figure(figsize=(12, 5))
plt.subplot(121)

fractions = []
random_fractions = []
comps_lists = []
improvement = []

for test_target in smiles_bits:

    # get all results except target being tested
    smiles_bits_reduced = get_reduced_smiles_bits(smiles_bits, test_target)
    all_smiles = [i for i in get_all_smiles(smiles_bits_reduced) if i in frequent_comps]

    # rank compounds based on past results
    comps, x = rank(smiles_bits_reduced, all_smiles, lib_size)
    comps_lists += comps

    # calculate total information from test target from full library
    total_information = get_total_information(smiles_bits, test_target)

    # get list of fraction of information recovered at each library size (ranked library)
    target_fraction = get_fraction(comps, smiles_bits, test_target, total_information)
    fractions.append(target_fraction)

    # add in dummy compounds for those we have never seen bind and shuffle new (larger) library
    add_for_DSiP = len(frequent_comps) - len(all_smiles)
    all_smiles += ['C'] * add_for_DSiP
    random.shuffle(all_smiles)

    # get list of fraction of information recovered at each library size (random library)
    random_fractions_target = multiple_random_runs(all_smiles, smiles_bits, test_target, total_information, num_runs=1000)
    random_fractions.append(random_fractions_target)

    # print and save improvement in information recovered from random and ranked sub-libaries
    print(test_target, len(comps), len(smiles_bits[test_target]), target_fraction[-1], random_fractions_target[-1])
    improvement.append(target_improvement(target_fraction, random_fractions_target))

    # plot information recovery as library size increases 
    plot_target_fraction(test_target, target_fraction)


mean_fractions = mean_across_targets(fractions)
mean_random_fractions = mean_across_targets(random_fractions)

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
plt.xlabel(f'number of times fragments in top {lib_size}')
plt.tight_layout()
plt.savefig('figures/compound_counts.png')


plt.close()
plt.figure(figsize=(8, 7))
plt.barh([t for t in smiles_bits], improvement)
plt.xlabel(f'factor of information improvement in DSiP-diverse {lib_size} vs average of 1,000 random runs')
plt.ylabel('target')
plt.tight_layout()
plt.savefig('figures/improvement_by_target.png')
