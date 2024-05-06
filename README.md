# NPMMP-DB

This folder contatins codes and data uesd for NPMMPDB database. 
The codes were provided in both **Jupyter Notebook(.ipynb)** and **Python Script(.py)**

## Content
```
NPMMP-DB
├── Code
│   ├── ChEMBL_Python_Client.ipynb
│   ├── ChEMBL_Python_Client.py
│   ├── Molecular_Processing.ipynb
│   ├── Molecular_Processing.py
│   ├── Split_and_Count.ipynb
│   └── Split_and_Count.py
├── Data
│   ├── csv
│   ├── img
│   └── sdf
└── README.md
```

## Descriptions

### Code

- **ChEMBL_Python_Client** : Python code we used for collecting macrocyclic molecular permeability data from ChEMBL. [(Davies et al., 2015)](https://academic.oup.com/nar/article/43/W1/W612/2467881)

  
- **Molecular_Processing**: We standardised the SMILES and gave an unique ID to each individual molecule. Then macrocycle, amide bond and serveral molecular feataures were calulated. We also create png images (saved in `/Data/img`) and sdf files(saved in `/Data/sdf`) for each molecules. Finally the overall data was saved in a csv file (saved in `/Data/csv`).
  
- **Split_and_Count**: The overall data were divided into different subsets accroding to their assays and endpoints. The subsets were saved in `/Data/csv`

### Data

- csv:
- img
- sdf

