import os
import json
from rank_compounds import get_next_best_compound, get_all_smiles, rank, get_DSiP_smiles
import matplotlib.pyplot as plt
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
from structurally_diverse_library import MACCS_diverse

import numpy as np


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


def ignore_targets(smiles_bits, frequent_comps, target_screens, remove_less_than):

    # returns a dictionary of IFP bits for each compound in each target that have more than a threshold amount of frequent compounds tested

    smol = []
    for t in smiles_bits:
        if t in target_screens and len([i for i in target_screens[t] if i in frequent_comps]) < remove_less_than:
            smol.append(t)
        elif t not in target_screens:
            smol.append(t)
        elif len(smiles_bits[t]) < 10:
            smol.append(t)
            print(t, len(smiles_bits[t]))

    for t in smol:
        del smiles_bits[t]


    return smiles_bits


def get_total_information(smiles_bits, target):

    # returns number of unique IFP bits for target

    total_information = []
    for c in smiles_bits[target]:
        total_information += smiles_bits[target][c]

    total_information = len(set(total_information))

    return total_information


def multiple_random_runs(all_smiles, smiles_bits, test_target, total_information, lib_size=150, num_runs=1000):

    # generates results for screens on given target using random libraries

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

    # plots information recovery for each target

    plt.plot(range(len(fraction)), fraction, label=target)
    plt.xlabel('library size')
    plt.ylabel('information recovered')
    plt.ylim(0,1)
    plt.tight_layout()
    # plt.savefig('figures/fraction_test_300.png')


def target_improvement(target_fraction, random_fractions_target):

    # calculates improvement in information recovered between ranked library and random library (mean of many runs)

    assert len(target_fraction) == len(random_fractions_target)

    if random_fractions_target[-1] > 0:
        imp = (target_fraction[-1]-random_fractions_target[-1])/random_fractions_target[-1]
    else:
        imp = 0
    
    return imp


def mean_across_targets(fraction):

    # calculates the mean information recovered across all targets for a particular library type

    mean_fractions = []
    fraction_lengths = [len(i) for i in fraction]
    for i in range(min(fraction_lengths)):
        mean_fractions.append(sum([x[i] for x in fraction])/len(fraction))
    
    return mean_fractions



smiles_bits = json.load(open('data/datafiles/smiles_bits_clean.json', 'r')) # load in dictionary of targets with IFP bits attributed to each compound
frequent_comps = json.load(open('data/datafiles/frequently_tested_compounds.json', 'r')) # load in list of compounds that have been tested on more than 15 of our targets
lib_size = 110 # library size we are using to compare ranked and random libraries
runs = 10
target_screens = json.load(open('data/datafiles/target_full_screens.json', 'r')) # load in dictionary showing which compounds were tested on which targets

DSiP_smiles = get_DSiP_smiles()
smiles_bits = ignore_targets(smiles_bits, frequent_comps, target_screens, 500) # only take targets that have had more than 250 of our frequent compounds tested
print(len(set(frequent_comps)))
plt.figure(figsize=(12, 5))
plt.subplot(121)


ranked_fractions = []
random_fractions = []
diverse_fractions_cross = []
# diverse_fractions = []

comps_lists = []
improvement = []
improvement_diverse = []

bars = []

for test_target in smiles_bits:

    # get all results except target being tested
    smiles_bits_reduced = get_reduced_smiles_bits(smiles_bits, test_target)
    all_smiles = [i for i in get_all_smiles(smiles_bits_reduced) if Chem.MolToSmiles(Chem.MolFromSmiles(i)) in frequent_comps]

    # calculate total information from test target from full library
    total_information = get_total_information(smiles_bits, test_target)

    fractions = []
    for i in range(runs):
        random.shuffle(all_smiles)
        # rank compounds based on past results
        comps, x = rank(smiles_bits_reduced, all_smiles, lib_size)
        comps_lists += comps

        # get list of fraction of information recovered at each library size (ranked library)
        target_fraction = get_fraction(comps, smiles_bits, test_target, total_information)
        fractions.append(target_fraction)
    
    target_fractions = []
    for i in range(lib_size+1):
        fracs = [f[i] for f in fractions]
        target_fractions.append(sum(fracs)/runs)

    fractions.append(target_fractions)   


    # get compounds that were actually tested on this target
    all_smiles = [i for i in target_screens[test_target] if i in frequent_comps]

    # get list of fraction of information recovered at each library size (random library)
    random_fractions_target = multiple_random_runs(all_smiles, smiles_bits, test_target, total_information, lib_size=lib_size, num_runs=runs)
    random_fractions.append(random_fractions_target)

    # print and save improvement in information recovered from random and ranked sub-libaries
    # print(test_target, len([i for i in smiles_bits[test_target] if Chem.MolToSmiles(Chem.MolFromSmiles(i)) in frequent_comps]), target_fraction[-1], random_fractions_target[-1])
    improvement.append(target_improvement(target_fraction, random_fractions_target))

    # plot information recovery as library size increases 
    plot_target_fraction(test_target, target_fractions)

    ## !!! lots of runs for structurally diverse library!!!!
    many_diverse_fractions = []
    for i in range(runs):
        random.shuffle(all_smiles)
        div = MACCS_diverse(all_smiles, lib_size)
        many_diverse_fractions.append(get_fraction(div, smiles_bits, test_target, total_information))
    diverse_fractions = []
    for i in range(lib_size+1):
        fracs = [f[i] for f in many_diverse_fractions]
        diverse_fractions.append(sum(fracs)/runs)


    print(test_target, len([i for i in smiles_bits[test_target] if Chem.MolToSmiles(Chem.MolFromSmiles(i)) in frequent_comps]), target_fractions[-1], random_fractions_target[-1], diverse_fractions[-1])
    improvement_diverse.append(target_improvement(target_fractions, diverse_fractions))
    diverse_fractions_cross.append(diverse_fractions)
    bars.append([target_fractions[-1], random_fractions_target[-1], diverse_fractions[-1]])


mean_fractions = mean_across_targets(fractions)
mean_random_fractions = mean_across_targets(random_fractions)
mean_diverse_fractions = mean_across_targets(diverse_fractions_cross)


plt.subplot(122)
plt.plot(range(len(mean_fractions)), mean_fractions, label='using ranked compounds')
plt.plot(range(len(mean_random_fractions)), mean_random_fractions, label='picking random compounds')
plt.plot(range(len(mean_diverse_fractions)), mean_diverse_fractions, label='picking diverse compounds')
plt.ylim(0, 1)
plt.legend()
plt.xlabel('library size')
plt.ylabel('mean information recovered across all new targets')
plt.tight_layout()
# plt.savefig('figures/fraction_test_300.png')


# all_smiles = get_all_smiles(smiles_bits)
# counts = [comps_lists.count(i) for i in all_smiles]
# plt.close()
# plt.figure(figsize=(8, 5))
# plt.bar(list(set(counts)), [counts.count(i) for i in set(counts)])
# plt.xlabel(f'number of times fragments in top {lib_size}')
# plt.tight_layout()
# plt.savefig('figures/compound_counts.png')

plt.close()
plt.figure(figsize=(10, 8))
fig, ax = plt.subplots()
pos = np.arange(len(smiles_bits))
ax.barh(pos+0.15, improvement, height=0.25, label='vs random runs')
ax.barh(pos-0.15, improvement_diverse, height=0.25, label='vs MACCS-diverse libraries')
ax.set_yticks(pos)
ax.set_yticklabels([t for t in smiles_bits])
plt.xlabel(f'mean factor of information improvement in DSiP-diverse {lib_size} across 1000 runs')
plt.ylabel('target')
plt.legend()
plt.tight_layout()
# plt.savefig('figures/improvement_by_target_300.png')

bar_rank = [i[0] for i in bars]
bar_random = [i[1] for i in bars]
bar_diverse = [i[2] for i in bars]

plt.close()
plt.figure(figsize=(10, 8))
fig, ax = plt.subplots()
pos = np.arange(len(smiles_bits))
ax.bar(pos-0.25, bar_rank, width=0.25, label='top-ranked compounds')
ax.bar(pos, bar_random, width=0.25,label='random compounds')
ax.bar(pos+0.25, bar_diverse, width=0.25,label='diverse compounds')
ax.set_xticks(pos)
ax.set_xticklabels([t for t in smiles_bits])
plt.xticks(rotation=90)
plt.ylabel(f'information recovered with library sizes of {lib_size}')
plt.xlabel('target')
plt.legend()
plt.tight_layout()
plt.savefig('figures/results_bar.png')

