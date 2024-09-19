'''
    For this homework, we will do something slightly different.
    In the next lab, we will start doing QSAR modeling, and we 
    need to have a dataset to work with. So, we will prepare our
    dataset in this homework. 

    Pick a dataset from the data folder, and using the functions
    we have discussed in the notebook, curate your dataset. You 
    can decide on the curation steps you want to take, but keep
    in mind that to avoid errors and have a good final QSAR model,
    you need to have a clean dataset.
'''


import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdDistGeom


## curation functions
def convert_to_logS(log_sol, mw):
    # solubility is in log(ug/mL)
    sol = 10**log_sol           # ug/mL
    uM = sol * 1000 / mw        # uM
    logS = np.log10(uM*1e-6)         # log(uM)
    return logS

def sanitize(mol):
    flag = Chem.SanitizeMol(mol, catchErrors=True)
    if flag != Chem.SanitizeFlags.SANITIZE_NONE:
        return False
    return True

def neutralize_mol(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def embed_3d(mol):
    tmp_mol = Chem.AddHs(mol)
    flag_3d = False
    try:
        rdDistGeom.EmbedMolecule(tmp_mol)
        rdDistGeom.MMFFOptimizeMolecule(tmp_mol)
        flag_3d = True
    except:
        pass

    return flag_3d

def curate_mixture(mol, n_frags=2):
    num_frags = len(Chem.GetMolFrags(mol))
    return num_frags <= n_frags

def keep_largest_fragment(mol):
    frags = Chem.GetMolFrags(mol)
    if len(frags) == 1:
        return mol
    frag_sizes = [frag.GetNumAtoms() for frag in frags]
    largest_frag = frags[np.argmax(frag_sizes)]
    return largest_frag

def curate_inorganic(mol):
    valid_atoms = ["H", "C", "N", "O", "P", "S", "CL", "F", "I", "BR", "B"]
    flag_organic = True
    for atom in mol.GetAtoms():
        if atom.GetSymbol().toUpperCase() not in valid_atoms:
            flag_organic = False
            break

    return flag_organic

def curate_boron(mol):
    flag_boron = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol().toUpperCase() == "B":
            flag_boron = True
            break

    return not flag_boron


# TODO: load your dataset, pick only one
# Hint: you can copy paste from the notebook


# TODO: curation steps
# Hint: you can do the same as in the notebook


# save the curated dataset
# df.to_csv("data/curated_dataset.csv", index=False)

