"""计算药物的 Tanimoto 相似度"""
import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs
import scipy.io as sio

path = 'C:/Users/1/Desktop/664_drugs/drug_info_664.csv'

smiles = pd.read_csv(path, header=0)['SMILES']
# print(smiles)
sim = np.zeros((len(smiles), len(smiles)))
for i, s1 in enumerate(smiles):
    if i % 10 == 0: print(i, end=' ')
    mol1 = Chem.MolFromSmiles(s1)
    fp1 = Chem.RDKFingerprint(mol1)
    for j, s2 in enumerate(smiles):
        mol2 = Chem.MolFromSmiles(s2)
        fp2 = Chem.RDKFingerprint(mol2)
        tani = DataStructs.TanimotoSimilarity(fp1, fp2)
        sim[i, j] = tani
print(sim)
print(np.sum(sim == 0))
# sio.savemat('Tani_sim_664.mat', {'tani': sim})

