# NPMMP-DB

This folder contatins codes and data uesd for Non-peptidic Macrocycle Membrane Permeability (NPMMP) database construc. 

The codes were provided in both **Jupyter Notebook(.ipynb)** and **Python Script(.py)**.

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
│   ├── Raw.csv
│   ├── csv
│   ├── img
│   └── sdf
└── README.md


```

## Working Environment
This projected was opreated on 
- macOS Sonoma (Version 14.4.1)
- Python (Version 3.7.16)
- RDKit (Version 2020.09.1)

## Descriptions

### Code

- **ChEMBL_Python_Client** : Python code we used for collecting macrocyclic molecular permeability data from ChEMBL. [(Davies et al., 2015)](https://academic.oup.com/nar/article/43/W1/W612/2467881)

  
- **Molecular_Processing**: We standardised the SMILES and gave an unique ID to each individual molecule from the raw data(`/Data/Raw.csv`). Then macrocycle, amide bond and serveral molecular feataures were calulated. We also create png images (saved in `/Data/img`) and sdf files(saved in `/Data/sdf`) for each molecules. Finally the overall data was saved in a csv file (saved in `/Data/csv`).
  
- **Split_and_Count**: The overall data were divided into different subsets accroding to their assays and endpoints. The subsets were saved in `/Data/csv`

### Data

- csv: Include overall data calculated by `/Code/Molecular_Processing.py` and subsets created by `/Code/Split_and_Count.py`
- img: png images of each molecules created by `/Code/Molecular_Processing.py`
- sdf: sdf files of each molecules created by `/Code/Molecular_Processing.py`


