from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover, rdMolDescriptors, MACCSkeys
from rdkit import SimDivFilters
import numpy as np
import os
from rank_compounds import get_all_smiles, get_DSiP_smiles
import json


def dmat_sim(fps,ntopick):
    ds=[]
    for i in range(1,len(fps)):
         ds.extend(DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i],returnDistance=True))
    mmp =SimDivFilters.MaxMinPicker()
    ids=mmp.Pick(np.array(ds),len(fps),ntopick)
    return ids



smiles_bits = json.load(open('data/datafiles/smiles_bits.json', 'r'))
all_smiles = get_all_smiles(smiles_bits)
DSiP_smiles = json.load(open('data/datafiles/frequently_tested_compounds.json', 'r'))
ms = [Chem.MolFromSmiles(i) for i in DSiP_smiles]
# fps = [rdMolDescriptors.GetMorganFingerprintAsBitVect(m,2) for m in ms]
fps = [MACCSkeys.GenMACCSKeys(m) for m in ms]

dmat_ids=dmat_sim(fps,250)
MaxMinSmiles = [DSiP_smiles[i] for i in list(dmat_ids)]
MaxMinMols = [ms[i] for i in list(dmat_ids)]

MaxMinMACCS = [MACCSkeys.GenMACCSKeys(x) for x in MaxMinMols]
MACCS_sim = []

for i in range(len(MaxMinMACCS)):
    for j in range(i+1, len(MaxMinMACCS)):
        MACCS_sim.append(DataStructs.FingerprintSimilarity(MaxMinMACCS[i],MaxMinMACCS[j]))

import matplotlib.pyplot as plt

plt.hist(MACCS_sim)
plt.savefig('figures/MACCS_sim.png')

MaxMinLib = []
print(len(MaxMinSmiles))

all_fps = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(sm),2) for sm in all_smiles]
MaxMin_fps = [AllChem.GetMorganFingerprint(Chem.MolFromSmiles(sm),2) for sm in MaxMinSmiles]
for fp1 in range(len(all_fps)):
    for fp2 in range(len(MaxMin_fps)):
        mol_sim = DataStructs.DiceSimilarity(all_fps[fp1], MaxMin_fps[fp2])
        if mol_sim == 1.0:
            MaxMinLib.append(all_smiles[fp1])
            break

to_add = 250-len(MaxMinLib)
MaxMinLib += ['C'] * to_add

print(len(MaxMinLib))

json.dump(MaxMinLib, open('data/datafiles/MaxMinLib.json', 'w'))
json.dump(all_smiles, open('data/datafiles/all_smiles.json', 'w'))