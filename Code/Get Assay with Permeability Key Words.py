#!/usr/bin/env python
# coding: utf-8

# *Searching data on ChEMBL using its Python web resource client. **PAMPA** is used as an example.*
# 
# *It can be replaced by other assays (like Caco-2, MDCK, RRCK...)*

# ### Import packages and client

# In[ ]:


import csv
from chembl_webresource_client.new_client import new_client

# Create a new client
client = new_client


# ### Searching by keywords

# In[ ]:


# Search for PAMPA assays
assays = client.assay.filter(description__icontains='PAMPA')

# Open a new CSV file for writing
with open('pampa_data.csv', 'w', newline='') as csvfile:
    fieldnames = ['CHEMBL_ID', 'SMILES', 'PAMPA_value', 'PAMPA_units', 'PAMPA_doi', 'PAMPA_assay_description']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()

    # For each assay, get the activities (which includes PAMPA values)
    for assay in assays:
        activities = client.activity.filter(assay_chembl_id=assay['assay_chembl_id'])

        # For each activity, save the CHEMBL ID, SMILES string, value, units, doi, and assay description
        for activity in activities:
            molecule = client.molecule.get(activity['molecule_chembl_id'])
            smiles = None
            if molecule['molecule_structures'] is not None:
                smiles = molecule['molecule_structures']['canonical_smiles']
            try:
                document = client.document.get(activity['document_chembl_id'])
                doi = document['doi'] if document['doi'] else "Not Available"
            except Exception:
                doi = "Not Available"
            writer.writerow({
                'CHEMBL_ID': activity['molecule_chembl_id'],
                'SMILES': smiles,
                'PAMPA_value': activity['value'],
                'PAMPA_units': activity['units'],
                'PAMPA_doi': doi,
                'PAMPA_assay_description': assay['description'] if assay['description'] else "Not Available"
            })


# ### Fliter the macrocycle molecules

# In[ ]:


from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
from rdkit.Chem.Draw import IPythonConsole
import pandas as pd


# In[ ]:


sample_data = pd.read_csv('pampa_data.csv')
comp_ids = sample_data["CHEMBL_ID"].to_list()
smiles = sample_data["SMILES"].to_list() 

def get_max_ring_size(mol):
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    if len(atom_rings) == 0:  # if there are no rings, return 0
        return 0
    else:
        maxring = max(len(ra) for ra in atom_rings)
        return maxring
    
molecules_ring_size = []
for comp_id, smile in  zip(comp_ids, smiles):
    mol = Chem.MolFromSmiles(str(smile))
    if mol is not None:
        ring_size = get_max_ring_size(mol)
        molecules_ring_size.append(ring_size)
    else:
        print(f"Invalid SMILES: {smile}")
        molecules_ring_size.append(None)  


# In[ ]:


sample_data['Ring_Size'] = molecules_ring_size
macrocycles = sample_data[sample_data["Ring_Size"] >=12] #Attention! Should be ">="" 12
macrocycles.to_csv("macrocycle_pampa.csv")

