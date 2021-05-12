import json
import matplotlib.pyplot as plt
import numpy as np 
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
import os

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

def get_DSiP_smiles():

    # using the sdf files provided by Enamine, the molecules are cleaned of salts and converted into SMILES strings

    DSiP = []

    for filename in os.listdir('data/DSiP'):
        DSiP += list(Chem.SDMolSupplier(f'data/DSiP/{filename}'))

    remover = SaltRemover.SaltRemover()
    DSiP = [remover.StripMol(mol) for mol in DSiP]
    DSiP_smiles = [Chem.MolToSmiles(i) for i in DSiP]
    DSiP_smiles = list(set(DSiP_smiles))

    return DSiP_smiles

def get_all_smiles(smiles_bits):

    # gets all smiles that have bound to proteins 

    all_smiles = []
    for i in smiles_bits:
        all_smiles += smiles_bits[i]

    all_smiles = list(set(all_smiles))

    return all_smiles


def get_next_best_compound(smiles_bits, current_info, current_compounds, all_compounds):

    # for each library size, selects new compound that gives most additional information
    
    new_compounds = [i for i in all_compounds if i not in current_compounds]

    best_new_info = 0
    best_new_compound = None
    best_new_bits = {}


    for compound in new_compounds:
        new_info = 0

        for target in smiles_bits:
            if compound in smiles_bits[target]:
                new_info += len([i for i in smiles_bits[target][compound] if i not in current_info[target]])

        if new_info > best_new_info:
            best_new_info = new_info
            best_new_compound = compound
            best_new_bits = {}
            for target in smiles_bits:
                if compound in smiles_bits[target]:
                    best_new_bits[target] = [i for i in smiles_bits[target][compound] if i not in current_info[target]]

        
    for target in best_new_bits:
        current_info[target] += best_new_bits[target]

    return current_info, best_new_compound        
        

def get_info_size(info_dict):

    # returns number of unique interactions across all targets 

    info_num = 0

    for target in info_dict:
        info_num += len(set(info_dict[target]))

    return info_num



def get_random_fraction(library, smiles_bits):

    random.shuffle(library)

    line_fraction = [0]

    for s in range(len(library)):
        info = []
        k = library[:s]
        for t in smiles_bits:
            info_targ = []
            for sm in smiles_bits[t]:
                if sm in k:
                    info_targ += smiles_bits[t][sm]

            info.append(len(set(info_targ)))
        
        line_fraction.append(sum(info))

    return line_fraction


def rank(smiles_bits, all_smiles, length):
    fraction = [0]
    comps = []

    current_info = {}
    for t in smiles_bits:
        current_info[t] = []

    for s in range(length):

        current_info, best_new_compound = get_next_best_compound(smiles_bits, current_info, comps, all_smiles)
        comps.append(best_new_compound)
        total_info_current = get_info_size(current_info)
        fraction.append(total_info_current)

    return comps, fraction



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
    # smiles_bits = json.load(open('data/datafiles/smiles_bits.json', 'r'))
    smiles_bits = json.load(open('data/datafiles/smiles_bits_clean.json', 'r')) # load in dictionary of targets with IFP bits attributed to each compound
    frequent_comps = json.load(open('data/datafiles/frequently_tested_compounds.json', 'r')) # load in list of compounds that have been tested on more than 15 of our targets
    target_screens = json.load(open('data/datafiles/target_full_screens.json', 'r')) # load in dictionary showing which compounds were tested on which targets


    lib_size = 300 # library size we are using to compare ranked and random libraries
    coverage = 500

    DSiP_smiles = get_DSiP_smiles()
    smiles_bits = ignore_targets(smiles_bits, frequent_comps, target_screens, coverage) # only take targets that have had more than 250 of our frequent compounds tested
    all_smiles = get_all_smiles(smiles_bits)
    # DSiP_smiles = all_smiles + ['C'] * 450
    # random.shuffle(DSiP_smiles)

    comps, fraction = rank(smiles_bits, all_smiles, len(DSiP_smiles))

    # json.dump(comps, open('data/outputs/ranked_compounds.json', 'w'))

    # random_full_fraction = get_random_fraction(DSiP_smiles, smiles_bits)

    random_smaller_fraction = []
    for i in range(1000):
        print(i)
        random_smaller_fraction.append(get_random_fraction(all_smiles, smiles_bits))

    random_smaller_fraction, std = mean_across_runs(random_smaller_fraction)
    # does_order_matter = get_random_fraction(comps, smiles_bits)

    fraction  = [i/max(fraction) for i in fraction]
    # random_full_fraction = [i/max(random_full_fraction) for i in random_full_fraction]
    random_smaller_fraction = [i/max(random_smaller_fraction) for i in random_smaller_fraction]
    # does_order_matter = [i/max(does_order_matter) for i in does_order_matter]


    # print(len(fraction), len(random_full_fraction), len(random_smaller_fraction))
    plt.close()
    plt.figure(figsize=(8,5))
    plt.plot(range(364), fraction[:364], label='one at a time, taking best fragment')
    # plt.plot(range(len(DSiP_smiles)+1), random_full_fraction, label='random order')
    plt.plot(range(364), random_smaller_fraction, label="random order with only compounds weve seen bind")
    # plt.plot(range(len(comps)+1), does_order_matter, label='random order with only compounds weve seen bind')
    plt.legend()
    plt.xlabel('library size')
    plt.ylabel('information recovered')
    plt.tight_layout()
    plt.savefig('figures/fraction.png')


    plt.close()

    print(len([i for i in all_smiles if i not in DSiP_smiles]))
    print(len(all_smiles))