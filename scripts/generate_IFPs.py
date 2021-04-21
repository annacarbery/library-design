import os
from oddt.fingerprints_new import InteractionFingerprint, tanimoto
import oddt
import sys 
from pymol import cmd
import statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
import json


def separate_files(filepath):

    # takes in a protein-ligand complex and separates files into protein and ligand pdb files

    cmd.reinitialize()
    cmd.load(filepath, 'complex')
    cmd.select('lig', 'resn LIG')
    cmd.save('data/tmp/lig.pdb', 'lig')
    cmd.extract('hets', 'complex and HETATM')
    cmd.save('data/tmp/prot.pdb', 'complex')


def get_IFP():

    # uses the previously split protein and ligand files and caculates binary protein-ligand interaction fingerprint

    lig = next(oddt.toolkit.readfile('pdb', 'data/tmp/lig.pdb'))
    prot = next(oddt.toolkit.readfile('pdb', 'data/tmp/prot.pdb'))

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


def get_IFP_vectors(target):

    # all IFPs for a particular target are calculated

    ismiles = []
    ifrags = []
    ivecs = []

    for ligand in os.listdir(f'{DATA_DIR}/{target}'):

        try:
            if len(xtal_smiles[ligand]) > 1:
                xtal_smiles[ligand] = [min(xtal_smiles[ligand])]
            if Chem.MolToSmiles(Chem.MolFromSmiles(xtal_smiles[ligand][0])) in DSiP_smiles:

                separate_files(f'{DATA_DIR}/{target}/{ligand}')

                IFP = get_IFP()
            
                if list(IFP).count(0) < len(IFP):
                    ismiles.append(xtal_smiles[ligand][0])
                    ifrags.append(ligand)
                    ivecs.append(IFP)

        
        except:
            pass
    
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


DATA_DIR = '/dls/science/users/tyt15771/DPhil/Lib_activity/data'

target_data = {}
DSiP_smiles = get_DSiP_smiles()
print('DSiP total compounds:', len(set(DSiP_smiles)))

xtal_smiles = json.load(open('data/datafiles/xtal_smiles.json', 'r'))


for target in os.listdir(DATA_DIR):
# for target in ['70X']:

    try:
        print(target)

        ifrags, ivecs, ismiles = get_IFP_vectors(target)
        vecs, frags, smiles, wrong = get_uniform_IFPs(ifrags, ivecs, ismiles)

        print('IFPs:', len(vecs))
        print('vectors of wrong length:', wrong)

        smiles_bits = get_smiles_bits(vecs, smiles)

        print('structures:', len(vecs), ', unique smiles:', len(smiles_bits))
        target_data[target] = smiles_bits
        json.dump(target_data, open('data/datafiles/smiles_bits.json', 'w'))
    
    except:
        print(target, 'error')
        print(sys.exc_info()[1])

