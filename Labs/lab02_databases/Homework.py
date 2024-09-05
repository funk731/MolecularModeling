'''
    Let's use what we have learned in Databases lab. For this homework,
    we will pull data for a target protein. I already chose a UniProt ID
    that we will be collecting information about. Let's start ...

    Remember you only need to edit the code when you see a TODO or ...
'''

# import necessary libraries
import pandas as pd
from chembl_webresource_client.new_client import new_client
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

# EGFR kinase 
uniprot_id = "P00533"

# define target and activity APIs
target_api = new_client.target
activity_api = new_client.activity

# Homework question 1: How many targets are pulled from ChEMBL
# you will ONLY need information on "organism", "pref_name", 
# "target_type", and "target_chembl_id"
targets = target_api.get(target_components__accession=uniprot_id).only(
    "target_chembl_id", ... # TODO: add more information to include
)
targets = pd.DataFrame.from_records(targets)
print(targets)
print(f"Number of Records: {len(targets)}")
print(f"Number of unique Target Types: {len(targets['target_type'].unique())}")
print(f"Unique target types: {targets['target_type'].unique()}")


# Homework question 2: Since there are multiple records for the target,
# we will choose the first one to proceed.
target = targets.iloc[0]
print(f"\nTarget to collect bioactivities:\n{target}", end='\n\n')
chembl_id = target.target_chembl_id
print(f"Selected target ChEMBL ID is: {chembl_id}")

# Now, we want to fetch the related bioactivity data.
# First, use the filter method to retrieve only the data
# you need. Also, obtain ONLY the following information:
# "activity_id", "assay_chembl_id", "assay_description",
# "assay_type", "molecule_chembl_id", "type", "relation",
# "target_chembl_id", "target_organism", "standard_units",
# and "standard_value"
bioactivities = activity_api.filter(
    target_chembl_id=chembl_id, 
    # TODO fill following parameters 
    type=...,
    relation=...,
    assay_type=..., 
).only(
    "activity_id",
    # TODO: add more information to the data
    ...
)
bioactivities = pd.DataFrame.from_dict(bioactivities)
print(bioactivities.head())
print(f"Number of activity records: {len(bioactivities)}")

# Next step would be curation and preprocessing of the data
# for visualizations and analysis, but we'll come to that later

