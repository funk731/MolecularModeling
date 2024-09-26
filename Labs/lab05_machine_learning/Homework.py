'''
    Finally, we get to do some machine learning! In this homework, 
    we will use the curated dataset you prepared in the last homework
    to build a QSAR model. 

    In the lab notebook, we discussed the different types of machine
    learning models we can use for QSAR modeling. For this homework,
    choose one model that was not discussed in the lab from sklearn
    and build a QSAR model for a regression task. Plot the predicted
    values against the true values to see how well your model performs.
    You will submit your evaluation metrics (R^2, MSE, MAE) and the plot
    of predicted vs true values in the homework submission.
'''

# Import necessary libraries
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

# Load dataset
df = pd.read_csv('data/curated_solubility_dataset.csv')

# Convert Mols to molecular fingerprints
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
def mols_to_fingerprints(mols, fp_gen=mfpgen):
    feature_vectors = []
    for mol in mols:
        fp = mfpgen.GetFingerprint(mol)
        feature_vectors.append(fp)
    return np.array(feature_vectors)

df['Fingerprint'] = mols_to_fingerprints([Chem.MolFromSmiles(smi) for smi in df['SMILES']]).tolist()

# Prepare features (X) and target (y) for regression
X = np.array(df['Fingerprint'].tolist())
y = df['LogS']

# Split into train and test sets
X_train, X_test, y_train_reg, y_test_reg = train_test_split(X, y, test_size=0.2, random_state=42)

# Standardize the features
scaler = MinMaxScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)


# TODO: Choose a model from sklearn and build a QSAR model for regression


# TODO: Evaluate the model using R^2, MSE, and MAE metrics 
#      and plot the predicted vs true values

