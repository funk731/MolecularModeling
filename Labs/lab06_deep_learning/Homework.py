'''
    Similar to the previous homework, we will use the curated dataset
    to build a QSAR model. This time, we will use a deep learning model
    to build a regression model.

    In the lab notebook, we discussed the different types of deep learning
    models we can use for QSAR modeling. For this homework, choose one
    of the deep learning models we discussed in the lab and build a QSAR
    model for a regression task. Plot the predicted values against the
    true values to see how well your model performs. You will submit your
    evaluation metrics (R^2, MSE, MAE) and the plot of predicted vs true
    values in the homework submission.
'''

# Import necessary libraries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Load dataset
df = pd.read_csv('data/curated_solubility_dataset.csv')

# to make the code faster, let's use only a subset of the data
df = df.sample(frac=0.1, random_state=42)

# Convert Mols to molecular fingerprints
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
def mols_to_fingerprints(mols, fp_gen=mfpgen, counts=False):
    feature_vectors = []
    for mol in mols:
        if counts:
            fp = fp_gen.GetCountFingerprint(mol)
        else:
            fp = fp_gen.GetFingerprint(mol)
        feature_vectors.append(fp)
    return np.array(feature_vectors)

df['Fingerprint'] = mols_to_fingerprints([Chem.MolFromSmiles(smi) for smi in df['SMILES']]).tolist()

# Prepare features (X) and target (y) for regression and classification
X = np.array(df['Fingerprint'].tolist())
X_smiles = df['SMILES'].values
y_regression = df['LogS']

# Split into train and test sets
X_train_reg, X_test_reg, y_train_reg, y_test_reg = train_test_split(X, y_regression, test_size=0.2, random_state=42)
X_train_smiles, X_test_smiles, y_train_smiles, y_test_smiles = train_test_split(X_smiles, y_regression, test_size=0.2, random_state=42)


import keras
from keras import layers

# TODO: Build your neural network model for regression


# TODO: Evaluate the model using R^2, MSE, and MAE metrics 
#      and plot the predicted vs true values

