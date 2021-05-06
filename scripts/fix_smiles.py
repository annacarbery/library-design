import json
from rdkit import Chem

smiles_bits = json.load(open('data/datafiles/smiles_bits.json', 'r')) # load in dictionary of targets with IFP bits attributed to each compound

smiles_bits_new = {}

for target in smiles_bits:
    smiles_bits_new[target] = {}
    for smiles in smiles_bits[target]:
        smiles_bits_new[target][Chem.MolToSmiles(Chem.MolFromSmiles(smiles))] = smiles_bits[target][smiles]


json.dump(smiles_bits_new, open('smiles_bits_clean.json', 'w'))