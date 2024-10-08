# -*- coding: utf-8 -*-
"""Molecular_Processing_Class.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/16lFVNqsVCKydbTNaZxLwnR4z4qFNissE
"""

from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import (Draw, PandasTools, AllChem, rdDepictor, rdMolDescriptors,
                        rdFMCS, rdmolops, Descriptors)
from rdkit.Chem.Draw import IPythonConsole, rdMolDraw2D
from molvs import standardize_smiles
import pandas as pd
from tqdm import tqdm
import io
from PIL import Image
import os
rdDepictor.SetPreferCoordGen(True)

sample_data = pd.read_csv('Data_2023July_2024August.csv',encoding='utf_8_sig')

def standardise(smiles):
    std_smiles = standardize_smiles(smiles)
    return standardize_smiles (smiles)

sample_data['Standardise_SMILES'] = sample_data['SMILES'].apply(standardise)
idx = sample_data.columns.get_loc('SMILES') + 1
# Remove the column and reinsert it to the desired location
standard_smiles_col = sample_data.pop('Standardise_SMILES')
sample_data.insert(idx, 'Standardise_SMILES', standard_smiles_col)
standard_smiles_to_id = {}

def assign_id(standard_smiles):
    if standard_smiles not in standard_smiles_to_id:
        standard_smiles_to_id[standard_smiles] = 'MC-' + str(len(standard_smiles_to_id) +4592).zfill(4)
    return standard_smiles_to_id[standard_smiles]

# Assuming sample_data is your DataFrame and standard_smiles_to_id is your dictionary
sample_data['ID'] = sample_data['Standardise_SMILES'].apply(assign_id)
smiles_list = sample_data["Standardise_SMILES"].tolist()

# convert SMILES to RDKit molecule objects



# Move the ID col to the first row
cols = sample_data.columns.tolist()
cols = [cols[-1]] + cols[:-1]  ## WARNING: DO NOT RUN THIS PART OF CODE MORE THAN ONCE!!
sample_data = sample_data[cols]
smiles_list = sample_data['Standardise_SMILES']
mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
def count_aromatic_rings(mol):
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    return len(aromatic_rings)
sample_data['Num_Aromatic_Rings'] = [count_aromatic_rings(mol) for mol in mols]

def _get_macrocycle_ring_mol(mol, strip=False):
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

class process:
    def __init__(self, list_of_smiles):
        self.smiles = list_of_smiles
        #self.standardized_smiles = [standardize_smiles(smi) for smi in list_of_smiles]
        self.mols = [Chem.MolFromSmiles(smi) for smi in self.smiles]

    def get_macrocycle_ring_size(self):
        ring_sizes = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=True)
            ring_size = macrocycle_mol.GetNumBonds()
            ring_sizes.append(ring_size)
        return ring_sizes

    def get_macrocycle_ring_size_for_mol (self,mol):
        macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=True)
        ring_size = macrocycle_mol.GetNumBonds()
        return ring_size


    def get_free_amide_count(self):
        free_amide_counts = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
            free_amide_smart = Chem.MolFromSmarts('[Nh][CX3]=[O]')
            free_amide_count = len(macrocycle_mol.GetSubstructMatches(free_amide_smart))
            free_amide_counts.append(free_amide_count)
        return free_amide_counts

    def get_free_amide_count_for_mol(self, mol):
        macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
        free_amide_smart = Chem.MolFromSmarts('[Nh][CX3]=[O]')
        return len(macrocycle_mol.GetSubstructMatches(free_amide_smart))

    def get_sub_amide_count(self):
        sub_amide_counts = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
            free_amide_count = self.get_free_amide_count_for_mol(mol)
            all_amide_smart = Chem.MolFromSmarts('[#6](=[#8])-[#7]-[#6]')
            all_amide_count = len(macrocycle_mol.GetSubstructMatches(all_amide_smart))
            sub_amide_count = (all_amide_count - free_amide_count) / 2
            sub_amide_counts.append(sub_amide_count)
        return sub_amide_counts

    def get_sub_amide_count_for_mol(self, mol):
        macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
        free_amide_count = self.get_free_amide_count_for_mol(mol)
        all_amide_smart = Chem.MolFromSmarts('[#6](=[#8])-[#7]-[#6]')
        all_amide_count = len(macrocycle_mol.GetSubstructMatches(all_amide_smart))
        sub_amide_count = (all_amide_count - free_amide_count) / 2
        return sub_amide_count

    def get_overall_amide_count(self):
        overall_amide_counts = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
            free_amide_count = self.get_free_amide_count_for_mol(mol)
            sub_amide_count = self.get_sub_amide_count_for_mol(mol)
            overall_amide_count = free_amide_count + sub_amide_count
            overall_amide_counts.append(overall_amide_count)
        return overall_amide_counts

    def get_macrocycle_core_smiles(self):
        core_smiles = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=True)
            macrocycle_core_smiles = Chem.MolToSmiles(macrocycle_mol)
            core_smiles.append(macrocycle_core_smiles)
        return core_smiles

    def get_macrocycle_peripheral_smiles(self):
        peripheral_smiles = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
            macrocycle_peripheral_smiles = Chem.MolToSmiles(macrocycle_mol)
            peripheral_smiles.append(macrocycle_peripheral_smiles)
        return peripheral_smiles

    def get_macrocycle_free_amide_ratio (self):
        macrocycle_free_amide_ratios = []
        for mol in tqdm(self.mols):
            free_amide_count = self.get_free_amide_count_for_mol(mol)
            ring_size = self.get_macrocycle_ring_size_for_mol(mol)
            free_amide_ratio = free_amide_count*3/ring_size
            macrocycle_free_amide_ratios.append(free_amide_ratio)
        return macrocycle_free_amide_ratios

    def get_macrocycle_amide_ratio (self):
        macrocycle_amide_ratios = []
        for mol in tqdm(self.mols):
            macrocycle_mol = _get_macrocycle_ring_mol(mol, strip=False)
            free_amide_count = self.get_free_amide_count_for_mol(mol)
            sub_amide_count = self.get_sub_amide_count_for_mol(mol)
            overall_amide_count = free_amide_count + sub_amide_count
            ring_size = self.get_macrocycle_ring_size_for_mol(mol)
            amide_ratio = overall_amide_count*3/ring_size
            macrocycle_amide_ratios.append(amide_ratio)
        return macrocycle_amide_ratios

    def get_num_of_rings (self):
        num_rings = []
        for mol in tqdm(self.mols):
            num_ring = mol.GetRingInfo().NumRings()
            num_rings.append(num_ring)
        return num_rings

    def get_cLogP(self):
        cLogP = []
        for mol in tqdm(self.mols):
            cLogP.append(Descriptors.MolLogP(mol))
        return cLogP

    def get_molecular_weight(self):
        molecular_weights = []
        for mol in tqdm(self.mols):
            molecular_weight = Descriptors.MolWt(mol)
            molecular_weights.append(molecular_weight)
        return molecular_weights

    def get_Num_H_Acceptors(self):
        Num_H_Acceptors = []
        for mol in tqdm(self.mols):
            Num_H_Acceptor = Descriptors.NumHAcceptors(mol)
            Num_H_Acceptors.append(Num_H_Acceptor)
        return Num_H_Acceptors

    def get_Num_H_donors(self):
        Num_H_donors = []
        for mol in tqdm(self.mols):
            Num_H_donor = Descriptors.NumHDonors(mol)
            Num_H_donors.append(Num_H_donor)
        return Num_H_donors

    def get_Num_Heavy_Atoms(self):
        Num_Heavy_Atoms = []
        for mol in tqdm(self.mols):
            Num_Heavy_Atom = Descriptors.HeavyAtomCount(mol)
            Num_Heavy_Atoms.append(Num_Heavy_Atom)
        return Num_Heavy_Atoms

    def get_Num_Carbon_Atoms(self):
        Num_Carbon_Atoms = []
        for mol in tqdm(self.mols):
            Num_Carbon_Atom = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
            Num_Carbon_Atoms.append(Num_Carbon_Atom)
        return Num_Carbon_Atoms

    def get_fraction_sp3_carbons(self):
        fraction_sp3_carbons = []
        sp3_carbons = 0
        total_carbons = 0
        for mol in tqdm(self.mols):
            for atom in mol.GetAtoms():
              if atom.GetSymbol() == 'C':  # Check if the atom is a carbon atom
                total_carbons += 1
              if atom.GetHybridization() == Chem.HybridizationType.SP3:
                sp3_carbons += 1
            if total_carbons == 0:
              fraction_sp3_carbon = 0
            else:
              fraction_sp3_carbon = sp3_carbons / total_carbons
              
            fraction_sp3_carbons.append(fraction_sp3_carbon)
        return fraction_sp3_carbons

    def get_TPSA(self):
        TPSA = []
        for mol in tqdm(self.mols):
            tpsa = Descriptors.TPSA(mol)
            TPSA.append(tpsa)
        return TPSA

    def get_Num_Rotatable_Bonds(self):
        Num_Rotatable_Bonds = []
        for mol in tqdm(self.mols):
            Num_Rotatable_Bond = Descriptors.NumRotatableBonds(mol)
            Num_Rotatable_Bonds.append(Num_Rotatable_Bond)

        return Num_Rotatable_Bonds

    def get_Num_Charged_Atoms(self):
        Num_Charged_Atoms = []
        for mol in tqdm(self.mols):
            Num_Charged_Atom = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() != 0)
            Num_Charged_Atoms.append(Num_Charged_Atom)
        return Num_Charged_Atoms
    
    def get_Net_Charge(self):
        Net_Charge = []
        for mol in tqdm(self.mols):
            net_charge = rdmolops.GetFormalCharge(mol)
            Net_Charge.append(net_charge)
        return Net_Charge

    def get_Kier_index(self):
        Kier_index = []
        for mol in tqdm(self.mols):
            kier_index = Descriptors.Kappa1(mol)*Descriptors.Kappa2(mol)/Descriptors.HeavyAtomCount(mol)
            Kier_index.append(kier_index)

        return Kier_index


    def get_InchiKey (self):
        InchiKey = []
        for mol in tqdm(self.mols):
            inchi_key = Chem.inchi.MolToInchiKey(mol)
            InchiKey.append(inchi_key)
        return InchiKey


    def result(self):
        Macrocycle_Ring_Sizes = self.get_macrocycle_ring_size()
        Macrocycle_Free_Amide_Counts = self.get_free_amide_count()
        Macrocycle_Substituted_Amide_Count = self.get_sub_amide_count()
        Macrocycle_Overall_Amide_Count = self.get_overall_amide_count()
        Macrocycle_Core_Smiles = self.get_macrocycle_core_smiles()
        Macrocycle_Peripheral_Smiles = self.get_macrocycle_peripheral_smiles()
        Macrocycle_free_amide_ratios = self.get_macrocycle_free_amide_ratio()
        Macrocycle_amide_ratios = self.get_macrocycle_amide_ratio()
        Num_Rings = self.get_num_of_rings()
        #Num_Aromatic_Rings = self.get_num_of_aromatic_rings()
        cLogP = self.get_cLogP()
        Molecular_Weight = self.get_molecular_weight()
        Num_H_Acceptors = self.get_Num_H_Acceptors()
        Num_H_donors = self.get_Num_H_donors()
        Num_Heavy_Atoms = self.get_Num_Heavy_Atoms()
        Num_Carbon_Atoms = self.get_Num_Carbon_Atoms()
        fraction_sp3_carbons = self.get_fraction_sp3_carbons()
        TPSA = self.get_TPSA()
        Num_Rotatable_Bonds = self.get_Num_Rotatable_Bonds()
        Num_Charged_Atoms = self.get_Num_Charged_Atoms()
        Net_Charge = self.get_Net_Charge()
        Kier_index = self.get_Kier_index()
        InchiKey = self.get_InchiKey()



        result_df = pd.DataFrame({
            "Macrocycle_Ring_Size": Macrocycle_Ring_Sizes,
            "Macrocycle_Free_Amide_Count": Macrocycle_Free_Amide_Counts,
            "Macrocycle_Substituted_Amide_Count": Macrocycle_Substituted_Amide_Count,
            "Macrocycle_Overall_Amide_Count": Macrocycle_Overall_Amide_Count,
            "Macrocycle_Core_Smiles": Macrocycle_Core_Smiles,
            "Macrocycle_Peripheral_Smiles": Macrocycle_Peripheral_Smiles,
            "Macrocycle_free_amide_ratios":Macrocycle_free_amide_ratios,
            "Macrocycle_amide_ratios":Macrocycle_amide_ratios,
            "Num_Rings":Num_Rings,
            #"Num_Aromatic_Rings":Num_Aromatic_Rings,
            "cLogP":cLogP,
            "Molecular_Weight":Molecular_Weight,
            "Num_H_Acceptors":Num_H_Acceptors,
            "Num_H_donors":Num_H_donors,
            "Num_Heavy_Atoms":Num_Heavy_Atoms,
            "Num_Carbon_Atoms":Num_Carbon_Atoms,
            "Fraction_SP3_Carbons":fraction_sp3_carbons,
            "TPSA":TPSA,
            "Num_Rotatable_Bonds":Num_Rotatable_Bonds,
            "Num_Charged_Atoms":Num_Charged_Atoms,
            "Net_Charge":Net_Charge,
            "Kier_index":Kier_index,
            "InchiKey":InchiKey})
        return result_df

rdkit_featurizer = process(smiles_list)
rdkit_features = rdkit_featurizer.result()
result_df = pd.concat([sample_data, rdkit_features], axis =1)
result_df.to_csv("Data_2023July_2024August_Result.csv",encoding='utf_8_sig')

""" # Define directory name
dir_name = '/Desktop/Overall_images'

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

dir_name_2 = '/Desktop/Overall_sdf'

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
            Chem.MolToMolFile(mol, file_path_2)  # Save directly to file """