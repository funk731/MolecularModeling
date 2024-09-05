'''
    Let's start with a very simple task. Here is a SMILES string of a molecule.
    Find 15 different chemical/physical properties of this molecule using RDKit
    and print them out with their names. You can use the RDKit documentation
    to find the functions that you need to use. RDKit has different modules for
    different types of properties. You're free to use any of them. 

    I'd suggest that you start by googling rdkit descriptors or you can take a 
    look at this blog :)
    https://greglandrum.github.io/rdkit-blog/posts/2022-12-23-descriptor-tutorial.html
'''

# helper function to print the properties as a table
def print_dict_as_table(d: dict):
    print(f"{'Property'.ljust(30)}| {'Value'}")
    print('-' * 45)
    for key, value in d.items():
        print(f'{key.ljust(30)}| {value:.2f}')


everyones_favorite = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'


# Import Descriptors module from rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors

# Write a function that takes a SMILES string as input 
# and returns a dictionary of 15 different chemical/physical properties
# of the molecule. Your function should return a dictionary with the
# property names as keys and the property values as values. Example:
# {'Molecular weight': 194.19, 'LogP': 0.5, ...}
def get_my_properties(smiles: str) -> dict:
    # TODO: convert the smiles into RDKit mol object
    mol = ... 

    # TODO: calculate 15 properties which should include MW, QED, 
    # maximum partial charge, TPSA, and number of H acceptors
    # to answer the question on canvas. Here is one from me :)
    properties = {}
    properties['MolWt'] = Descriptors.MolWt(mol)
    ...
    return properties


# Check if the function works
properties = get_my_properties(everyones_favorite)
assert len(properties) >= 15, "You should return at least 15 properties"
print_dict_as_table(properties)
