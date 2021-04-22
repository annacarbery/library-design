import os
from oddt.fingerprints_new import InteractionFingerprint
import oddt
import sys 
from pymol import cmd
import statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
import json
import argparse


def separate_files(filepath):

    # takes in a protein-ligand complex and separates files into protein and ligand pdb files

    cmd.reinitialize()
    cmd.load(filepath, 'complex')
    cmd.select('lig', 'resn LIG')
    cmd.save('lig.pdb', 'lig')
    cmd.extract('hets', 'complex and HETATM')
    cmd.save('prot.pdb', 'complex')


def get_IFP():

    # uses the previously split protein and ligand files and caculates binary protein-ligand interaction fingerprint

    lig = next(oddt.toolkit.readfile('pdb', 'lig.pdb'))
    prot = next(oddt.toolkit.readfile('pdb', 'prot.pdb'))

    prot.protein = True

    IFP = InteractionFingerprint(lig, prot, strict=False)

    return IFP


def get_DSiP_smiles():

    # using the sdf files provided by Enamine, the molecules are cleaned of salts and converted into SMILES strings

    DSiP = []

    for filename in os.listdir('data/DSiP'):
        DSiP += list(Chem.SDMolSupplier(f'data/DSiP/{filename}'))

    remover = SaltRemover.SaltRemover()
    DSiP = [remover.StripMol(mol) for mol in DSiP]
    DSiP_smiles = [Chem.MolToSmiles(i) for i in DSiP]

    return DSiP_smiles


def get_IFP_vectors(input_data):

    # all IFPs for a particular target are calculated

    ismiles = []
    ifrags = []
    ivecs = []

    for ligand in sorted(os.listdir(input_data)):

        try:
            if os.path.exists(os.path.join(input_data, ligand, 'refine.pdb')) and 'LIG' in open(os.path.join(input_data, ligand, 'refine.pdb'), 'r').read():
                print(ligand)
            
                separate_files(os.path.join(input_data, ligand, 'refine.pdb'))

                IFP = get_IFP()

            
                if list(IFP).count(0) < len(IFP):
                    ismiles.append(xtal_smiles[ligand])
                    ifrags.append(ligand)
                    ivecs.append(IFP)
                

        
        except:
            # raise
            pass
            # print(sys.exc_info()[1])
    
    return ismiles, ifrags, ivecs


def get_uniform_IFPs(ismiles, ifrags, ivecs):

    # IFPs are analysed to determine most common length (due to some models missing residues)
    # only IFPs of identical length are returned

    smiles = []
    frags = []
    vecs = []

    lengths = [len(i) for i in ivecs]
    length = statistics.mode(lengths)
    print(length)

    wrong = 0
    for i in range(len(ivecs)):
        if len(ivecs[i]) == length:
            vecs.append(ivecs[i])
            frags.append(ifrags[i])
            smiles.append(ismiles[i])
        else:
            wrong += 1

    return vecs, frags, smiles, wrong


def get_smiles_bits(vecs, smiles):

    # take IFP vectors and assign 'on' bits to the smiles strings responsible

    smiles_bits = {}

    for i in range(len(smiles)):
        if smiles[i] not in smiles_bits:
            smiles_bits[smiles[i]] = []

        for f in range(len(smiles)):
            mol1, mol2 = Chem.MolFromSmiles(smiles[i]), Chem.MolFromSmiles(smiles[f])
            mol_sim = DataStructs.DiceSimilarity(AllChem.GetMorganFingerprint(mol1,2), AllChem.GetMorganFingerprint(mol2, 2))
            if mol_sim == 1:
                smiles_bits[smiles[i]] += [b for b in range(len(vecs[f])) if vecs[f][b] != 0]
        
    for smiles in smiles_bits:
        smiles_bits[smiles] = list(set(smiles_bits[smiles]))
    
    return smiles_bits

parser = argparse.ArgumentParser()
parser.add_argument('-input', help='model_building directory full path')
parser.add_argument('-smiles', help='csv with format Mpro-x0001,SMILES_STRING')
args = vars(parser.parse_args())

INPUT_DATA = args['input']

smiles_strings = open(args['smiles'], 'r').readlines()
xtal_smiles = {}
for line in smiles_strings:
    xtal_smiles[line.split('.pdb')[0]] = line.split(',')[1].strip()


ifrags, ivecs, ismiles = get_IFP_vectors(INPUT_DATA)
vecs, frags, smiles, wrong = get_uniform_IFPs(ifrags, ivecs, ismiles)

print('IFPs:', len(vecs))
print('vectors of wrong length:', wrong)

smiles_bits = get_smiles_bits(vecs, smiles)

print('structures:', len(vecs), ', unique smiles:', len(smiles_bits))
json.dump(smiles_bits, open('Mpro_bits.json', 'w'))
