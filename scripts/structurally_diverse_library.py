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


def MACCS_diverse(all_smiles, lib_size):
    ms = [Chem.MolFromSmiles(i) for i in all_smiles]
    fps = [MACCSkeys.GenMACCSKeys(m) for m in ms]
    dmat_ids=dmat_sim(fps,lib_size)
    MaxMinSmiles = [all_smiles[i] for i in list(dmat_ids)]
    # MaxMinMols = [ms[i] for i in list(dmat_ids)]

    return MaxMinSmiles
