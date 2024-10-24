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
import matplotlib.pyplot as plt
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
from sklearn.metrics import r2_score
from sklearn import metrics

# TODO: Build your neural network model for regression
reg_model = keras.Sequential([
    layers.Dense(256, activation='relu', input_shape=(X_train_reg.shape[1],)),
    layers.Dense(64, activation='relu'),
    layers.Dense(8, activation='relu'),
    layers.Dense(1)  # Regression output
])

# Compiling the regression model
reg_model.compile(optimizer='adam', loss='mean_squared_error', metrics=['mean_absolute_error'])

# Training the regression model
history_reg = reg_model.fit(X_train_reg, y_train_reg, epochs=10, validation_split=0.2)

# TODO: Evaluate the model using R^2, MSE, and MAE metrics 
#      and plot the predicted vs true values

y_pred_reg = reg_model.predict(X_test_reg).flatten()

r2 = metrics.r2_score(y_test_reg, y_pred_reg)
mae = metrics.mean_absolute_error(y_test_reg, y_pred_reg)
mse = metrics.mean_squared_error(y_test_reg, y_pred_reg)
plt.figure(figsize=(6, 6))
plt.scatter(y_test_reg, y_pred_reg, alpha=0.5)
plt.annotate(f'R2 Score: {r2:.2f}\nMAE: {mae:.2f}\nMSE: {mse:.2f}' , xy=(0.05, 0.9), xycoords='axes fraction')
plt.xlabel('True LogS')
plt.ylabel('Predicted LogS')
plt.title('Regression Model: True vs. Predicted LogS')
plt.show()
