import os
import json
from rank_compounds import get_next_best_compound, get_all_smiles, rank, get_DSiP_smiles
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
    random_runs = []
    for i in range(runs):
        random.shuffle(all_smiles)
        fraction = get_fraction(all_smiles[:lib_size], smiles_bits, test_target, total_information)
        random_runs.append(fraction)

    target_fractions, target_std = mean_across_runs(random_runs)

    return target_fractions, target_std


def multiple_runs(all_smiles, smiles_bits, lib_size, test_target, total_information, runs=1000):
   
    ranked_runs = []

    for i in range(runs):

        random.shuffle(all_smiles)
        # rank compounds based on past results
        comps, x = rank(smiles_bits_reduced, all_smiles, lib_size)
        # comps_lists += comps

        # get list of fraction of information recovered at each library size (ranked library)
        target_fraction = get_fraction(comps, smiles_bits, test_target, total_information)
        ranked_runs.append(target_fraction)

    target_fractions, ranked_std = mean_across_runs(ranked_runs)

    return target_fractions, ranked_std


def multiple_diverse_runs(all_smiles, smiles_bits, lib_size, test_target, total_information, runs=1000):
    
    diverse_runs = []
    
    for i in range(runs):
        random.shuffle(all_smiles)
        div = MACCS_diverse(all_smiles, lib_size)
        diverse_runs.append(get_fraction(div, smiles_bits, test_target, total_information))

    target_fractions, target_std = mean_across_runs(diverse_runs)

    return target_fractions, target_std


def mean_across_runs(list_of_runs):

    # calculates mean and standard deviations of runs for a particular method

    mean_of_runs = []
    std_of_runs = []
    for i in range(len(list_of_runs[0])):
        fracs = [f[i] for f in list_of_runs]
        mean_of_runs.append(sum(fracs)/len(list_of_runs))
        std_of_runs.append(np.std(fracs))
    
    return mean_of_runs, std_of_runs
    

if __name__ == '__main__':

    smiles_bits = json.load(open('data/datafiles/smiles_bits_atomic.json', 'r')) # load in dictionary of targets with IFP bits attributed to each compound
    frequent_comps = json.load(open('data/datafiles/frequently_tested_compounds.json', 'r')) # load in list of compounds that have been tested on more than 15 of our targets
    target_screens = json.load(open('data/datafiles/target_full_screens.json', 'r')) # load in dictionary showing which compounds were tested on which targets

    lib_size = 300 # library size we are using to compare ranked and random libraries
    runs = 1000
    coverage = 500

    DSiP_smiles = get_DSiP_smiles()
    smiles_bits = ignore_targets(smiles_bits, frequent_comps, target_screens, coverage) # only take targets that have had more than 250 of our frequent compounds tested
    print(len(set(frequent_comps)))


    ranked_fractions = []
    ranked_stds = []

    random_fractions = []
    random_stds = []

    diverse_fractions = []
    diverse_stds = []

    for test_target in smiles_bits:

        # get all results except target being tested
        smiles_bits_reduced = get_reduced_smiles_bits(smiles_bits, test_target)
        all_smiles = [i for i in get_all_smiles(smiles_bits_reduced) if Chem.MolToSmiles(Chem.MolFromSmiles(i)) in frequent_comps and i in target_screens[test_target]]
        print(len(all_smiles))
        # calculate total information from test target from full library
        total_information = get_total_information(smiles_bits, test_target)

        # calculate information recovered and standard deviation over many runs for ranked library
        target_fractions, ranked_std = multiple_runs(all_smiles, smiles_bits, lib_size, test_target, total_information, runs)
        ranked_fractions.append(target_fractions) 
        ranked_stds.append(ranked_std)

        # get compounds that were actually tested on this target
        all_smiles = [i for i in target_screens[test_target] if i in frequent_comps]

        # get list of fraction of information recovered at each library size (random library)
        random_fractions_target, random_std = multiple_random_runs(all_smiles, smiles_bits, test_target, total_information, lib_size=lib_size, num_runs=runs)
        random_fractions.append(random_fractions_target)
        random_stds.append(random_std)

        # lots of runs for structurally diverse library
        diverse_runs, diverse_std = multiple_diverse_runs(all_smiles, smiles_bits, lib_size, test_target, total_information, runs)
        diverse_fractions.append(diverse_runs)
        diverse_stds.append(diverse_std)

        print(test_target, len([i for i in smiles_bits[test_target] if Chem.MolToSmiles(Chem.MolFromSmiles(i)) in frequent_comps]))
        print(target_fractions[-1], ranked_std[-1])
        print(random_fractions_target[-1], random_std[-1])
        print(diverse_runs[-1], diverse_std[-1])


    json.dump(ranked_fractions, open('data/outputs/ranked_fractions_atomic.json', 'w'))
    json.dump(random_fractions, open('data/outputs/random_fractions_atomic.json', 'w'))
    json.dump(diverse_fractions, open('data/outputs/diverse_fractions_atomic.json', 'w'))

    json.dump(ranked_stds, open('data/outputs/ranked_stds_atomic.json', 'w'))
    json.dump(random_stds, open('data/outputs/random_stds_atomic.json', 'w'))
    json.dump(diverse_stds, open('data/outputs/diverse_stds_atomic.json', 'w'))