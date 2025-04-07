import os
import torch
import numpy as np

import argparse
import h5py
import pandas as pd
from tqdm.auto import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem

########### define your constants ###########
# atom types in zinc
atoms_zc = ['C', 'O', 'N', 'F', 'S', 'Cl', 'Br', 'I']
# atom types in TarDrug3D
#157->B, 185->Si, 240->Pt, 253->Si, Na (removed)
atoms_td = ['C', 'O', 'N', 'F', 'S', 'Cl', 'Br', 'I', 'P']
seed = 42  # refer to difflinker

'''
def trans_to_MolAnd3D(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('invalid smiles provided.')
    mol = Chem.AddHs(mol)
    rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 10:
        maxconfs = 150
    elif rot_bonds <= 20:
        maxconfs = 300
    elif rot_bonds <= 30:
        maxconfs = 350
    elif rot_bonds <= 40:
        maxconfs = 400
    else:
        maxconfs = 750
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=maxconfs, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, randomSeed=seed, numThreads=1)
    converged_res = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=16)
    cenergy = [i[1] for i in converged_res]
    sortedcids = sorted(cids, key=lambda cid: cenergy[cid])
    mol = Chem.RemoveHs(mol)
    if len(sortedcids) > 0:
        lowest_conf = sortedcids[0]
        mol.SetProp('_Energy', str(cenergy[lowest_conf]))
        pos = torch.tensor(mol.GetConformer(lowest_conf).GetPositions())
    return mol, pos

'''
# for test process_protac() and process_adc() function
def trans_to_MolAnd3D(smiles):  # for test fun
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.Compute2DCoords(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True)
    AllChem.UFFOptimizeMolecule(mol)
    mol = Chem.RemoveHs(mol)
    pos = mol.GetConformer().GetPositions()

    return mol, pos


def process_protac(raw_file_dir, save_path):
    dataset = []
    n_fail = 0  # record fail num
    all_data = pd.read_excel(raw_file_dir, header=0)
    for i in tqdm(range(len(all_data))):
        id, protac_smi, warhead_smi, linker_smi, e3lig_smi = (all_data.iloc[i]['PROTAC ID'], all_data.iloc[i]['PROTAC Smiles'], all_data.iloc[i]['Warhead Smiles'], all_data.iloc[i]['Linker Smiles'], all_data.iloc[i]['E3 Ligand Smiles'])
        if protac_smi is None or warhead_smi is None or linker_smi is None or e3lig_smi is None:
            print(f'fail i: {i}')
            continue
        # get 3d coords
        try:
            protac_mol, atoms_pos_in_protac = trans_to_MolAnd3D(protac_smi)
        except ValueError as e:
            print(f'protac_id: {id}, error: {e}')
            continue
        # get atom types
        atoms_num = protac_mol.GetNumAtoms()
        atomic_numbers = torch.zeros(atoms_num, dtype=torch.long)  
        for j in range(atoms_num):
            atom = protac_mol.GetAtomWithIdx(j)
            atomic_numbers[j] = atom.GetAtomicNum()
            
        # get bond_info
        full_bonds_info = []
        for bond in protac_mol.GetBonds():
            start = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType().name
            full_bonds_info.append((start, end, bond_type))

        # get fragment and linker indices
        warhead_idxs = protac_mol.GetSubstructMatches(Chem.MolFromSmiles(warhead_smi, sanitize=False))
        linker_idxs = protac_mol.GetSubstructMatches(Chem.MolFromSmiles(linker_smi, sanitize=False))
        e3lig_idxs = protac_mol.GetSubstructMatches(Chem.MolFromSmiles(e3lig_smi, sanitize=False))
        if len(warhead_idxs) == 0:
            print(str(i) + 'warhead_idxs error !')
        if len(linker_idxs) == 0:
            print(str(i) + 'linker_idxs error !')
        if len(e3lig_idxs) == 0:
            print(str(i) + 'e3lig_idxs error !')
        flag = False
        for warhead_idx in warhead_idxs:
            for linker_idx in linker_idxs:
                for e3lig_idx in e3lig_idxs:
                    warhead_idx, linker_idx, e3lig_idx = list(warhead_idx), list(linker_idx), list(e3lig_idx)
                    if ((len(set(warhead_idx).intersection(set(linker_idx)))==0) and (len(set(linker_idx).intersection(set(e3lig_idx)))==0) and (len(set(e3lig_idx).intersection(set(warhead_idx)))==0) and (len(warhead_idx)+len(linker_idx)+len(e3lig_idx)==protac_mol.GetNumAtoms())):
                        flag = True
                        break
                    if flag:
                        break
                if flag:
                    break
            if flag:
                break
        if not flag:
            print(str(i) + 'validity check fail!')
            continue
        # get anchors
        anchors = torch.zeros(atoms_num, dtype=int)
        frag_mask = torch.zeros(atoms_num, dtype=int)
        frag_mask[warhead_idx] = 1
        frag_mask[e3lig_idx] = 2
        for bond in protac_mol.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if frag_mask[start_idx] != frag_mask[end_idx]:
                anchors[start_idx] = 1
                anchors[end_idx] = 1
        anchors *= frag_mask
        # get submol
        warhead_mol = Chem.PathToSubmol(protac_mol, warhead_idx)
        linker_mol = Chem.PathToSubmol(protac_mol, linker_idx)
        e3lig_mol = Chem.PathToSubmol(protac_mol, e3lig_idx)

        dataset.append({
            'id': id,
            'smiles': [protac_smi, warhead_smi, linker_smi, e3lig_smi],
            'pos': atoms_pos_in_protac,
            'atoms_atomic_numbers': atomic_numbers,
            'frag1_idx': warhead_idx,
            'frag2_idx': e3lig_idx,
            'linker_idx': linker_idx,
            'linker_size': len(linker_idx),
            'bond': full_bonds_info,
            'anchors': anchors.bool(),
        })
    print('n fail: ', n_fail)
    print('saving data')
    with h5py.File(save_path, 'w') as f:
        grp = f.create_group('protac')
    
        for idx, data in enumerate(dataset):
            sample_grp = grp.create_group(f'{idx}')
            sample_grp.attrs['id'] = data['id']
            sample_grp.attrs['linker_size'] = data['linker_size']
            smiles = np.array(data['smiles'], dtype=h5py.string_dtype(encoding='utf-8'))
            sample_grp.create_dataset('smiles', data=smiles)
            for array_field in ['pos', 'atoms_atomic_numbers', 'frag1_idx', 'frag2_idx', 'linker_idx', 'anchors']:
                if array_field in data:
                    sample_grp.create_dataset(array_field, data=data[array_field])
            bonds = data['bond']
            bond_dtype = np.dtype([
                ('atom1', 'i4'),  
                ('atom2', 'i4'),  
                ('bond_type', h5py.string_dtype(encoding='utf-8'))  
            ])
    
            bond_array = np.array(
                [(b[0], b[1], b[2]) for b in bonds],  
            dtype=bond_dtype
            )
            sample_grp.create_dataset('bond', data=bond_array)

    print('length processed data: ', f'{len(dataset)}')


def process_adc(raw_file_dir, save_path):
    dataset = []
    n_fail = 0  # record fail num
    all_data = pd.read_excel(raw_file_dir, header=0)
    all_data = all_data[0:4]
    for i in tqdm(range(len(all_data))):
        id, clp_smi, conjugate_smi, linker_smi, payload_smi = (all_data.iloc[i]['ADC ID'], all_data.iloc[i]['Conjugate-Linker-Payload Isosmiles'], all_data.iloc[i]['Conjugate Smiles'], all_data.iloc[i]['Post-reaction Linker Isosmiles'], all_data.iloc[i]['Post-reaction Payload Isosmiles'])
        if clp_smi is None or conjugate_smi is None or linker_smi is None or payload_smi is None:
            print(f'fail i: {i}')
            continue
        # get 3d coords
        try:
            clp_mol, atoms_pos_in_clp = trans_to_MolAnd3D(clp_smi)
        except ValueError as e:
            print(e)
        # get atom types
        atoms_num = clp_mol.GetNumAtoms()
        atomic_numbers = torch.zeros(atoms_num, dtype=torch.long)  
        for j in range(atoms_num):
            atom = clp_mol.GetAtomWithIdx(j)
            atomic_numbers[j] = atom.GetAtomicNum()
        # get bond_info
        full_bonds_info = []
        for bond in clp_mol.GetBonds():
            start = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            bond_type = bond.GetBondType().name
            full_bonds_info.append((start, end, bond_type))

        # get fragment and linker indices
        conjugate_idxs = clp_mol.GetSubstructMatches(Chem.MolFromSmiles(conjugate_smi))
        linker_idxs = clp_mol.GetSubstructMatches(Chem.MolFromSmiles(linker_smi))
        payload_idxs = clp_mol.GetSubstructMatches(Chem.MolFromSmiles(payload_smi))
        if len(conjugate_idxs) == 0:
            print(str(i) + 'conjugate_idxs error !')
        if len(linker_idxs) == 0:
            print(str(i) + 'linker_idxs error !')
        if len(payload_idxs) == 0:
            print(str(i) + 'payload_idxs error !')
        flag = False
        for conjugate_idx in conjugate_idxs:
            for linker_idx in linker_idxs:
                for payload_idx in payload_idxs:
                    conjugate_idx, linker_idx, payload_idx = list(conjugate_idx), list(linker_idx), list(payload_idx)
                    if ((len(set(conjugate_idx).intersection(set(linker_idx)))==0) and (len(set(linker_idx).intersection(set(payload_idx)))==0) and (len(set(payload_idx).intersection(set(conjugate_idx)))==0) and (len(conjugate_idx)+len(linker_idx)+len(payload_idx)==clp_mol.GetNumAtoms())):
                        flag = True
                        break
                    if flag:
                        break
                if flag:
                    break
            if flag:
                break
        if not flag:
            print(str(i) + 'validity check fail!')
            continue
        # get anchors
        anchors = torch.zeros(atoms_num, dtype=int)
        frag_mask = torch.zeros(atoms_num, dtype=int)
        frag_mask[conjugate_idx] = 1
        frag_mask[payload_idx] = 2
        for bond in clp_mol.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if frag_mask[start_idx] != frag_mask[end_idx]:
                anchors[start_idx] = 1
                anchors[end_idx] = 1
        anchors *= frag_mask
        # get submol
        conjugate_mol = Chem.PathToSubmol(clp_mol, conjugate_idx)
        linker_mol = Chem.PathToSubmol(clp_mol, linker_idx)
        payload_mol = Chem.PathToSubmol(clp_mol, payload_idx)

        dataset.append({
            'id': id,
            'smiles': [clp_smi, conjugate_smi, linker_smi, payload_smi],
            'pos': atoms_pos_in_clp,
            'atoms_atomic_numbers': atomic_numbers,
            'frag1_idx': conjugate_idx,
            'frag2_idx': payload_idx,
            'linker_idx': linker_idx,
            'linker_size': len(linker_idx),
            'bond': full_bonds_info,
            'anchors': anchors.bool(),
        })
        #print('n fail: ', n_fail)
    print('saving data')
    with h5py.File(save_path, 'w') as f:
        grp = f.create_group('adc')
    
        for idx, data in enumerate(dataset):
            sample_grp = grp.create_group(f'{idx}')
            sample_grp.attrs['id'] = data['id']
            sample_grp.attrs['linker_size'] = data['linker_size']
            smiles = np.array(data['smiles'], dtype=h5py.string_dtype(encoding='utf-8'))
            sample_grp.create_dataset('smiles', data=smiles)
            for array_field in ['pos', 'atoms_atomic_numbers', 'frag1_idx', 'frag2_idx', 'linker_idx', 'anchors']:
                if array_field in data:
                    sample_grp.create_dataset(array_field, data=data[array_field])
            bonds = data['bond']
            bond_dtype = np.dtype([
                ('atom1', 'i4'),  
                ('atom2', 'i4'),  
                ('bond_type', h5py.string_dtype(encoding='utf-8'))  
            ])
    
            bond_array = np.array(
                [(b[0], b[1], b[2]) for b in bonds],  
            dtype=bond_dtype
            )
            sample_grp.create_dataset('bond', data=bond_array)
    print('length processed data: ', len(dataset))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_path', type=str, required=True)
    parser.add_argument('--dataset', type=str, default='protac', choices=['protac', 'adc'])
    parser.add_argument('--dest', type=str, required=True)
    parser.add_argument('--mode', type=str, default='full', choices=['tiny', 'full'])  # for zinc only
    args = parser.parse_args()

    # example for run: python process_rawdata.py --raw_path tarDrug/protac.xlsx --dest tarDrug/protac.h5 --dataset protac
    if args.dataset == 'protac':
        process_protac(args.raw_path, args.dest)
    elif args.dataset == 'adc':
        process_adc(args.raw_path, args.dest)
    else:
        raise NotImplementedError