import json
import matplotlib.pyplot as plt
import numpy as np
from test_rankings import ignore_targets
from matplotlib.font_manager import FontProperties

def plot_target_fraction(target, fraction):

    # plots information recovery for each target

    plt.plot(range(len(fraction)), fraction, label=target)
    plt.xlabel('library size')
    plt.ylabel('information recovered')
    plt.ylim(0,1)
    plt.tight_layout()
    plt.savefig('figures/fraction_test_atomic.png')


def mean_across_targets(fraction):

    # calculates the mean information recovered across all targets for a particular library type

    mean_fractions = []
    fraction_lengths = [len(i) for i in fraction]
    for i in range(min(fraction_lengths)):
        mean_fractions.append(sum([x[i] for x in fraction])/len(fraction))
    
    return mean_fractions


def target_improvement(target_fraction_list, random_fractions_target_list, lib_size):

    # calculates improvement in information recovered between ranked library and random library (mean of many runs)
    improvement = []
    
    for i in range(len(target_fraction_list)):

        target_fraction = target_fraction_list[i]
        random_fractions_target = random_fractions_target_list[i]

        assert len(target_fraction) == len(random_fractions_target)

        if random_fractions_target[lib_size] > 0:
            imp = (target_fraction[lib_size]-random_fractions_target[lib_size])/random_fractions_target[lib_size]
        else:
            imp = 0
        
        improvement.append(imp)

    return improvement

def plot_fraction_and_std(mean_fractions, mean_stds, label, color):
    plt.plot(range(len(mean_fractions)), mean_fractions, label=label)
    plt.fill_between(range(len(mean_fractions)), 
                [mean_fractions[i]-mean_stds[i] for i in range(len(mean_fractions))], 
                [mean_fractions[i]+mean_stds[i] for i in range(len(mean_fractions))],
                color=color, alpha=0.1)

def plot_fractions(ranked_fractions, ranked_stds, random_fractions, random_stds, diverse_fractions, diverse_stds):

    mean_ranked_fractions = mean_across_targets(ranked_fractions)
    mean_ranked_stds = mean_across_targets(ranked_stds)

    mean_random_fractions = mean_across_targets(random_fractions)
    mean_random_stds = mean_across_targets(random_stds)

    mean_diverse_fractions = mean_across_targets(diverse_fractions)
    mean_diverse_stds = mean_across_targets(diverse_stds)

    plt.figure(figsize=(12, 6))
    plt.subplot(121)
    for i in ranked_fractions:
        plot_target_fraction('target', i)

    plt.subplot(122)

    plot_fraction_and_std(mean_ranked_fractions, mean_ranked_stds, 'ranked compounds', '#1f77b4')
    plot_fraction_and_std(mean_random_fractions, mean_random_stds, 'random compounds', '#ff7f0e')
    plot_fraction_and_std(mean_diverse_fractions, mean_diverse_stds, 'diverse compounds', '#2ca02c')

    plt.ylim(0, 1)
    plt.legend()
    plt.xlabel('library size')
    plt.ylabel('mean information recovered across all new targets')
    plt.tight_layout()
    plt.savefig('figures/fraction_test_atomic.png')


ranked_fractions = json.load(open('data/outputs/ranked_fractions_atomic.json', 'r'))
ranked_stds = json.load(open('data/outputs/ranked_stds_atomic.json', 'r'))

random_fractions = json.load(open('data/outputs/random_fractions_atomic.json', 'r'))
random_stds = json.load(open('data/outputs/random_stds_atomic.json', 'r'))

diverse_fractions = json.load(open('data/outputs/diverse_fractions_atomic.json', 'r'))
diverse_stds = json.load(open('data/outputs/diverse_stds_atomic.json', 'r'))

frequent_comps = json.load(open('data/datafiles/frequently_tested_compounds.json', 'r'))
target_screens = json.load(open('data/datafiles/target_full_screens.json', 'r'))
smiles_bits = json.load(open('data/datafiles/smiles_bits_atomic.json', 'r'))
smiles_bits = ignore_targets(smiles_bits, frequent_comps, target_screens, 500)


plot_fractions(ranked_fractions, ranked_stds, random_fractions, random_stds, diverse_fractions, diverse_stds)

lib_size = 100

bar_rank = [i[lib_size] for i in ranked_fractions]
error_rank = [i[lib_size] for i in ranked_stds]
bar_random = [i[lib_size] for i in random_fractions]
error_random = [i[lib_size] for i in random_stds]
bar_diverse = [i[lib_size] for i in diverse_fractions]
error_diverse = [i[lib_size] for i in diverse_stds]

plt.close()
plt.figure(figsize=(10, 6))
fig, ax = plt.subplots()
pos = np.arange(len(smiles_bits))
ax.bar(pos-0.25, bar_rank, yerr=error_rank, width=0.25, label='top-ranked compounds', alpha=0.7)
ax.bar(pos, bar_random, yerr=error_random, width=0.25,label='random compounds', alpha=0.7)
ax.bar(pos+0.25, bar_diverse, yerr=error_diverse, width=0.25,label='diverse compounds', alpha=0.7)
ax.set_xticks(pos)
ax.set_xticklabels([t for t in smiles_bits])
plt.xticks(rotation=90)
plt.ylabel(f'information recovered with library sizes of {lib_size}')
plt.xlabel('target')
plt.ylim(0, 0.65)
fontP = FontProperties()
fontP.set_size('x-small')
plt.legend(prop=fontP)
plt.tight_layout()
plt.savefig('figures/results_bar_atomic.png')

improvement_random = target_improvement(ranked_fractions, random_fractions, lib_size)
improvement_diverse = target_improvement(ranked_fractions, diverse_fractions, lib_size)

plt.close()
plt.figure(figsize=(10, 8))
fig, ax = plt.subplots()
pos = np.arange(len(smiles_bits))
ax.barh(pos+0.15, improvement_random, height=0.25, label='vs random runs', color= '#ff7f0e', alpha=0.7)
ax.barh(pos-0.15, improvement_diverse, height=0.25, label='vs MACCS-diverse libraries', color='#2ca02c', alpha=0.7)
ax.set_yticks(pos)
ax.set_yticklabels([t for t in smiles_bits])
plt.xlabel(f'mean factor of information improvement in DSiP-diverse {lib_size}')
plt.ylabel('target')
plt.legend(prop=fontP)
plt.tight_layout()
plt.savefig('figures/improvement_by_target_atomic.png')