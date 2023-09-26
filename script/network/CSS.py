import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit import DataStructs
from tqdm import tqdm


def get_molecualr_similarity(file, threshold):
    assert os.path.exists(file), "File does not existing, please read after cheaking!"
    a = pd.read_csv(file, sep='\t')
    mol_list, fp_list = [], []
    for smi in a.SMILES:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(smi)
            continue
        # mol_list = [Chem.MolFromSmiles(smi) for smi in a.SMILES]
        # fp_list = []

        fp_list.append(GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024))
    # fp_list = [GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024) for mol in mol_list]
    node1_smi, node2_smi = [], []
    node1_id, node2_id = [], []
    t_ls = []
    node1_type0, node1_type1, node1_type2 = [], [], []
    node2_type0, node2_type1, node1_type2 = [], [], []
    for i in tqdm(a.index):
        for j in range(i):
            tanimoto = DataStructs.TanimotoSimilarity(fp_list[i], fp_list[j])
            if tanimoto >= threshold:
                node1_smi.append(a.SMILES[i])
                node2_smi.append(a.SMILES[j])
                node1_type0.append(a.NAME[i])
                node2_type0.append(a.NAME[j])
                node1_type1.append(a.MIBIG[i])
                node2_type1.append(a.MIBIG[j])
                node1_id.append(i)
                node2_id.append(j)
                t_ls.append(tanimoto)

        '''else:
            node1_smi.append(a.SMILES[i])
            node2_smi.append(a.SMILES[i])
            node1_type0.append(a.CID[i])
            node2_type0.append(a.CID[i])
            node1_type1.append(a.CLASS[i])
            node2_type1.append(a.CLASS[i])
            node1_type2.append(a.SUBCLASS[i])
            node2_type2.append(a.SUBCLASS[i])
            node1_type3.append(a.PATHWAY[i])
            node2_type3.append(a.PATHWAY[i])
            node1_type4.append(a.MIBIG[i])
            node2_type4.append(a.MIBIG[i])
            node1_id.append(i)
            node2_id.append(i)
            t_ls.append(0)'''

    pd.DataFrame({
        'node1_id': node1_id,
        'node2_id': node2_id,
        'node1_smi': node1_smi,
        'node2_smi': node2_smi,
        'tanimoto': t_ls,
        'node1_type0': node1_type0,
        'node2_type0': node2_type0,
        'node1_type1': node1_type1,
        'node2_type1': node2_type1
    }).to_csv(f'Curvularin_t{threshold}.tsv', sep='\t', index=False)


if __name__ == '__main__':
    for i in [0.6]:
        get_molecualr_similarity("smiles_input.tsv", i)
