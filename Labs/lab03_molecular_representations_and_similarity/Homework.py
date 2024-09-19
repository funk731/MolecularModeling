'''
    Finally, this homework we will work on a small dataset
    of molecules. We will demonstrate how different molecular
    representations can result in different similarity scores.

    Let's start by loading the dataset and defining some helper
    functions.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import DistanceMetric
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator


# helper functions
def my_molecular_representation(mols):
    feature_vectors = []
    for mol in mols:
        mol_weight = Descriptors.MolWt(mol)
        log_p = Descriptors.MolLogP(mol)
        num_atoms = mol.GetNumAtoms()
        num_heavy_atoms = mol.GetNumHeavyAtoms()
        num_rings = mol.GetRingInfo().NumRings()
        num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        num_h_donors = Descriptors.NumHDonors(mol)
        num_h_acceptors = Descriptors.NumHAcceptors(mol)
        feature_vector = [mol_weight, log_p, num_atoms, num_heavy_atoms, 
                        num_rings, num_rotatable_bonds, num_h_donors, 
                        num_h_acceptors]
        feature_vectors.append(feature_vector)
    out = np.nan_to_num(feature_vectors)
    return out

def get_rdkit_descriptors(mols):
    feature_vectors = []
    for mol in mols:
        desc = []
        for desc_name, desc_func in Descriptors._descList:
            calculated_desc = desc_func(mol)
            desc.append(calculated_desc)
        feature_vectors.append(desc)
    out = np.nan_to_num(feature_vectors)
    return out

def get_morgan_fingerprints(mols):
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
    feature_vectors = []
    for mol in mols:
        fp = mfpgen.GetFingerprint(mol)
        feature_vectors.append(fp)
    return np.array(feature_vectors)

def get_rdkit_fingerprints(mols):
    rdkgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048, minPath=1, maxPath=7)
    feature_vectors = []
    for mol in mols:
        fp = rdkgen.GetFingerprint(mol)
        feature_vectors.append(fp)
    return np.array(feature_vectors)

def get_distance_matrix(feature_vectors, metric):
    dist = DistanceMetric.get_metric(metric)
    return dist.pairwise(feature_vectors)


# load the dataset
sample_smiles = pd.read_csv('data/homework_data.csv')
mols = [Chem.MolFromSmiles(smile) for smile in sample_smiles['SMILES']]
print(f"Number of molecules: {len(mols)}")

# calculate molecular representations
my_mol_rep = my_molecular_representation(mols)
rdkit_desc = get_rdkit_descriptors(mols)
morgan_fp = get_morgan_fingerprints(mols)
rdkit_fp = get_rdkit_fingerprints(mols)

# Question: What distance metric (e.g., euclidean, jaccard, etc.) 
# should we use for each molecular representation? Why?
# Hint: Think about the type of data each representation
# generates and how it can be compared.
# TODO: Assign the appropriate distance metric to each
my_mol_rep_dist = get_distance_matrix(my_mol_rep, "euclidean")
rdkit_desc_dist = get_distance_matrix(rdkit_desc, "minkowski")
morgan_fp_dist = get_distance_matrix(morgan_fp, "hamming")
rdkit_fp_dist = get_distance_matrix(rdkit_fp, "jaccard")


# plot and save the distance matrices
plt.figure(figsize=(12, 12))
sns.clustermap(my_mol_rep_dist, cmap='Reds_r', xticklabels=False, yticklabels=False)
plt.title("My Molecular Representation")
plt.savefig('data/my_mol_rep_dist.png')

plt.figure(figsize=(12, 12))
sns.clustermap(rdkit_desc_dist, cmap='Reds_r', xticklabels=False, yticklabels=False)
plt.title("RDKit Descriptors")
plt.savefig('data/rdkit_desc_dist.png')

plt.figure(figsize=(12, 12))
sns.clustermap(morgan_fp_dist, cmap='Reds_r', xticklabels=False, yticklabels=False)
plt.title("Morgan Fingerprints")
plt.savefig('data/morgan_fp_dist.png')

plt.figure(figsize=(12, 12))
sns.clustermap(rdkit_fp_dist, cmap='Reds_r', xticklabels=False, yticklabels=False)
plt.title("RDKit Fingerprints")
plt.savefig('data/rdkit_fp_dist.png')

print("Done!")

