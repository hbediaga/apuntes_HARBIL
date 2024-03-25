import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np

smiles_list = ['C(C(=O)O)N', 'C1=CC=C(C=C1)C[C@@H](C(=O)O)N ', 'C1=C(NC=N1)C[C@@H](C(=O)O)N', 'C([C@@H](C(=O)O)N)S']

mol_list = []
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    mol_list.append(mol)


glycine = mol_list[0]
cysteine = mol_list[3]

bi = {}
fp = AllChem.GetMorganFingerprintAsBitVect(glycine, 2, nBits=1024, bitInfo=bi)
fp_arr =np.zeros((1,))
DataStructs.ConvertToNumpyArray(fp, fp_arr)
np.nonzero(fp_arr)

bi2 = {}
fp2 = AllChem.GetMorganFingerprintAsBitVect(cysteine, 2, nBits=1024, bitInfo=bi)
fp_arr2 =np.zeros((1,))
DataStructs.ConvertToNumpyArray(fp2, fp_arr2)
np.nonzero(fp_arr2)

# Draw molecule fingerprints
# prints = [(glycine, x, bi) for x in fp.GetOnBits()]
# Draw.DrawMorganBits(prints, molsPerRow=4, legends=[str(x) for x in fp.GetOnBits()])

tanimoto = DataStructs.TanimotoSimilarity(fp,fp2)