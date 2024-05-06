#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import (Draw, PandasTools, AllChem, rdDepictor, rdMolDescriptors, 
                        rdFMCS, rdmolops, Descriptors)
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
from molvs import standardize_smiles
import pandas as pd
import io
from PIL import Image
import os
rdDepictor.SetPreferCoordGen(True)


# In[ ]:


sample_data = pd.read_csv('Raw.csv',encoding='utf_8_sig')
sample_data


# In[ ]:


def standardise(smiles):
    std_smiles = standardize_smiles(smiles)
    return standardize_smiles (smiles)

sample_data['Standardise_SMILES'] = sample_data['SMILES'].apply(standardise)
idx = sample_data.columns.get_loc('SMILES') + 1
# Remove the column and reinsert it to the desired location
standard_smiles_col = sample_data.pop('Standardise_SMILES')
sample_data.insert(idx, 'Standardise_SMILES', standard_smiles_col)
sample_data


# In[ ]:


standard_smiles_to_id = {}

def assign_id(standard_smiles):
    if standard_smiles not in standard_smiles_to_id:
        standard_smiles_to_id[standard_smiles] = 'MC-' + str(len(standard_smiles_to_id) + 4593).zfill(4) #Continue with MC-4592
    return standard_smiles_to_id[standard_smiles]

# Assuming sample_data is your DataFrame and standard_smiles_to_id is your dictionary
sample_data['ID'] = sample_data['Standardise_SMILES'].apply(assign_id)


# Move the ID col to the first row
cols = sample_data.columns.tolist()
cols = [cols[-1]] + cols[:-1]  ## WARNING: DO NOT RUN THIS PART OF CODE MORE THAN ONCE!!
sample_data = sample_data[cols]
sample_data = sample_data.sort_values(by='ID')


# In[ ]:


def get_macrocycle_ring_mol(smi, strip=False):
    mol = Chem.MolFromSmiles(smi)
    Chem.RemoveStereochemistry(mol)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    ri = mol.GetRingInfo()
    atoms_ring = max(ri.AtomRings(), key=len)
    if strip:
        macrocycle_smiles = AllChem.MolFragmentToSmiles(mol, atomsToUse=list(atoms_ring))
    else:
        query = Chem.MolFromSmarts('[!#1;!#7]~[R]*')
        matches = mol.GetSubstructMatches(query)
        atoms_to_use = list(atoms_ring) + [m[0] for m in matches]
        macrocycle_smiles = AllChem.MolFragmentToSmiles(mol, atomsToUse=atoms_to_use)
    macrocycle_mol = Chem.MolFromSmiles(macrocycle_smiles)
    mol_frags = rdmolops.GetMolFrags(macrocycle_mol, asMols=True)
    largest_mol = max(mol_frags, default=macrocycle_mol, key=lambda m: m.GetNumAtoms())
    return largest_mol


def get_macrocycle_ring_size(smi):
    macrocycle_mol = get_macrocycle_ring_mol(smi, strip = True)
    ring_size = macrocycle_mol.GetNumBonds()
    return ring_size

def get_free_amide_count(smi):
    macrocycle_mol = get_macrocycle_ring_mol(smi, strip=False)
    free_amide_smart = Chem.MolFromSmarts('[Nh][CX3]=[O]') 
    free_amide_count = len(macrocycle_mol.GetSubstructMatches(free_amide_smart))
    return free_amide_count

def get_sub_amide_count(smi):
    macrocycle_mol = get_macrocycle_ring_mol(smi, strip=False)
    free_amide_count = get_free_amide_count(smi)
    all_amide_smart = Chem.MolFromSmarts('[#6](=[#8])-[#7]-[#6]') 
    sub_amide_count = (len(macrocycle_mol.GetSubstructMatches(all_amide_smart))- free_amide_count)/2
    return sub_amide_count

def get_macrocycle_core_smiles(smi):
    macrocycle_mol = get_macrocycle_ring_mol(smi, strip = True)
    macrocycle_core_smiles = Chem.MolToSmiles(macrocycle_mol)
    return macrocycle_core_smiles

def get_macrocycle_peripheral_smiles(smi):
    macrocycle_mol = get_macrocycle_ring_mol(smi, strip = False)
    macrocycle_peripheral_smiles = Chem.MolToSmiles(macrocycle_mol)
    return macrocycle_peripheral_smiles


# In[ ]:


smiles = sample_data["Standardise_SMILES"].to_list() 

ring_size_list = []
free_amide_count_list = []
sub_amide_count_list = []
macrocycle_core_smiles_list = []
macrocycle_peripheral_smiles_list = []

for smile in smiles:
    ring_size = get_macrocycle_ring_size(smile)
    free_amide_count = get_free_amide_count(smile)
    sub_amide_count = int(get_sub_amide_count(smile))
    macrocycle_core_smiles = get_macrocycle_core_smiles(smile)
    macrocycle_peripheral_smiles = get_macrocycle_peripheral_smiles(smile)
    
    ring_size_list.append(ring_size)
    free_amide_count_list.append(free_amide_count)
    sub_amide_count_list.append(sub_amide_count)
    macrocycle_core_smiles_list.append(macrocycle_core_smiles)
    macrocycle_peripheral_smiles_list.append(macrocycle_peripheral_smiles)
    
sample_data['Macrocycle_Ring_Size'] = ring_size_list
sample_data['Macrocycle_Free_Amide_Count'] = free_amide_count_list
sample_data['Macrocycle_Substituted_Amide_Count'] = sub_amide_count_list
sample_data['Macrocycle_Overall_Amide_Count'] = sample_data['Macrocycle_Free_Amide_Count'] + sample_data['Macrocycle_Substituted_Amide_Count']
sample_data["Macrocycle_Core_smiles"] = macrocycle_core_smiles_list
sample_data["Macrocycle_Peripheral_smiles"] = macrocycle_peripheral_smiles_list
sample_data["Macrocycle_Free_Amide_Ratio"] = sample_data['Macrocycle_Free_Amide_Count']*3/sample_data['Macrocycle_Ring_Size']
sample_data['Macrocycle_Overall_Amide_Ratio'] = sample_data['Macrocycle_Overall_Amide_Count']*3/sample_data['Macrocycle_Ring_Size']


# In[ ]:


smiles_list = sample_data["Standardise_SMILES"].tolist()

# convert SMILES to RDKit molecule objects
mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]

# calculate descriptors and add them to the DataFrame
sample_data['Num_Rings'] = [mol.GetRingInfo().NumRings() for mol in mols]

# Get aromatic rings only
def count_aromatic_rings(mol):
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    return len(aromatic_rings)
sample_data['Num_Aromatic_Rings'] = [count_aromatic_rings(mol) for mol in mols]
sample_data['cLogP'] = [Descriptors.MolLogP(mol) for mol in mols]
sample_data['Molecular_Weight'] = [Descriptors.MolWt(mol) for mol in mols]
sample_data['Num_H_Acceptors'] = [Descriptors.NumHAcceptors(mol) for mol in mols]
sample_data['Num_H_Donors'] = [Descriptors.NumHDonors(mol) for mol in mols]
sample_data['Num_Heavy_Atoms'] = [Descriptors.HeavyAtomCount(mol) for mol in mols]

# Count Carbon atoms and Sp3 Carbon
def count_carbon_atoms(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
sample_data['Num_Carbon_Atoms'] = [count_carbon_atoms(mol) for mol in mols]

def fraction_sp3_carbons(mol):
    sp3_carbons = 0
    total_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':  # Check if the atom is a carbon atom
            total_carbons += 1
            if atom.GetHybridization() == Chem.HybridizationType.SP3:
                sp3_carbons += 1

    return sp3_carbons / total_carbons if total_carbons > 0 else 0


# Assuming 'mols' is a list of RDKit molecule objects
sample_data['Fraction_SP3_Carbons'] = [fraction_sp3_carbons(mol) for mol in mols]
sample_data['TPSA'] = [Descriptors.TPSA(mol) for mol in mols]
sample_data['Num_Rotatable_Bonds'] = [Descriptors.NumRotatableBonds(mol) for mol in mols]

# Number of charged atoms
def count_charged_atoms(mol):
    return sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0)
sample_data['Num_Charged_Atoms'] = [count_charged_atoms(mol) for mol in mols]
sample_data['Net_Charge'] = [rdmolops.GetFormalCharge(mol) for mol in mols]

sample_data["Kier_index"] = [(Descriptors.Kappa1(mol)*Descriptors.Kappa2(mol)/Descriptors.HeavyAtomCount(mol)) 
                             for mol in mols]
sample_data["InchiKey"] = [Chem.inchi.MolToInchiKey(mol) for mol in mols]


# In[ ]:


# Define directory name
dir_name = 'Desktop/Overall_images'

# Check if the directory exists, if not create it
if not os.path.exists(dir_name):
    os.makedirs(dir_name)

# Assuming 'sample_data' is defined
# Drop duplicates
sample_data.drop_duplicates(subset='Standardise_SMILES', keep="first", inplace=True)

# Create Images and save in target folder
for index, row in sample_data.iterrows():
    molecule = row['Standardise_SMILES']
    if molecule is not None:  
        drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
        drawer.drawOptions().clearBackground = False
        drawer.drawOptions().addStereoAnnotation = False
        drawer.DrawMolecule(Chem.MolFromSmiles(molecule))
        drawer.FinishDrawing()
        img_data = drawer.GetDrawingText()  # Get image data as bytes
        file_path = os.path.join(dir_name, f"{row['ID']}.png")
        with open(file_path, 'wb') as f:
            f.write(img_data)  # Write the bytes to a file
    else:
        print(f"SMILES string at index {index} could not be converted into a molecule.")


# In[ ]:


dir_name_2 = 'Desktop/Overall_sdf'

# Check if the directory exists, if not create it
if not os.path.exists(dir_name_2):
    os.makedirs(dir_name_2)

# Create Images and save in target folder
for index, row in sample_data.iterrows():
    molecule = row['Standardise_SMILES']
    if molecule:
        mol = Chem.MolFromSmiles(molecule)  # Assuming you need to convert SMILES to a molecule
        if mol is not None:
            file_path_2 = os.path.join(dir_name_2, f"{row['ID']}.sdf")
            Chem.MolToMolFile(mol, file_path_2)  # Save directly to file


# In[ ]:


sample_data.to_csv('New.csv',encoding='utf_8_sig',index=False)

