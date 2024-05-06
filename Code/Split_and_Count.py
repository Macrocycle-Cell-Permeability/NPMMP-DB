#!/usr/bin/env python
# coding: utf-8

# *After Molecular processing, separate the overall dataset into subsets.*

# In[ ]:


import pandas as pd

# Load the dataset  
df = pd.read_csv('Overall.csv', encoding='utf_8_sig')


# #### Group by 'Assay' only and create subsets

# In[ ]:


grouped = df.groupby(['Assay'])
subsets = {group: data for group, data in grouped}
# Create a DataFrame to store the counts
counts_assay = pd.DataFrame(columns=['File', 'Number of Entries'])
# Optionally, to save these subsets as separate CSV files
for (assay), subset in subsets.items():
    file_name = f'{assay}.csv'
    subset.to_csv(file_name, index=False, encoding='utf_8_sig')
    # Add the count information to the DataFrame
    counts_assay = counts_assay.append({'File': file_name, 'Number of Entries': len(subset)}, ignore_index=True)
# Display the counts DataFrame
counts_assay_2 = counts_assay.sort_values(by = ['Number of Entries'],ascending=False)
counts_assay_2


# #### Group by Both 'Endpoint' and 'Assay' and create subsets

# In[ ]:


grouped = df.groupby(['Standardized_Endpoint', 'Assay'])
subsets = {group: data for group, data in grouped}
counts_sub = pd.DataFrame(columns=['File', 'Number of Entries'])
for (std_endpoint,assay), subset in subsets.items():
    file_name = f'{assay}  {std_endpoint}.csv'
    subset.to_csv(file_name, index=False, encoding='utf_8_sig')
    counts_sub = counts_sub.append({'File': file_name, 'Number of Entries': len(subset)}, ignore_index=True)
counts_sub_2 = counts_sub.sort_values(by = ['Number of Entries'],ascending=False)
counts_sub_2

