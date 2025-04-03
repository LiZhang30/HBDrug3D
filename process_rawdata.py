import os
import torch
import numpy as np

import argparse
import pickle
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
'''

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
        try:
            atoms_num = protac_mol.GetNumAtoms()
            atom2idx = {atom: idx for idx, atom in enumerate(atoms_td)}
            one_hot = torch.zeros((atoms_num, len(atoms_td)), dtype=torch.float)
            for j in range(atoms_num):
                atom = protac_mol.GetAtomWithIdx(j)
                atom_symbol = atom.GetSymbol()
                if atom_symbol in atom2idx:
                    one_hot[j, atom2idx[atom_symbol]] = 1.0
                else:
                    print('find new atom type!')
                    raise ValueError('unknown atom type encountered!')
        except ValueError:
            continue
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
            'atom_one_hot': one_hot,
            'warhead_idx': warhead_idx,
            'e3lig_idx': e3lig_idx,
            'linker_idx': linker_idx,
            'linker_size': len(linker_idx),

            # reserved field
            'mol': protac_mol,
            'bond': full_bonds_info,
            'warhead_mol': warhead_mol,
            'e3lig_mol': e3lig_mol,
            'linker_mol': linker_mol,
            'anchors': anchors.bool(),
        })
    print('n fail: ', n_fail)
    print('saving data')
    with open(save_path, 'wb') as f:
        pickle.dump(dataset, f)
    print('length processed data: ', f'{len(dataset)}')

def process_adc(raw_file_dir, save_path):
    dataset = []
    n_fail = 0  # record fail num
    all_data = pd.read_excel(raw_file_dir, header=0)
    for i in tqdm(range(len(all_data))):
        id, clp_smi, conjugate_smi, linker_smi, payload_smi = (all_data.iloc[i]['ADC ID'], all_data.iloc[i]['Conjugate-Linker-Payload Smiles'], all_data.iloc[i]['Conjugate Smiles'], all_data.iloc[i]['Linker Smiles'], all_data.iloc[i]['Payload Smiles'])
        if clp_smi is None or conjugate_smi is None or linker_smi is None or payload_smi is None:
            print(f'fail i: {i}')
            continue
        # get 3d coords
        try:
            clp_mol, atoms_pos_in_clp = trans_to_MolAnd3D(clp_smi)
        except ValueError as e:
            print(e)
        # get atom types
        try:
            atoms_num = clp_mol.GetNumAtoms()
            atom2idx = {atom: idx for idx, atom in enumerate(atoms_td)}
            one_hot = torch.zeros((atoms_num, len(atoms_td)), dtype=torch.float)
            for j in range(atoms_num):
                atom = clp_mol.GetAtomWithIdx(j)
                atom_symbol = atom.GetSymbol()
                if atom_symbol in atom2idx:
                    one_hot[j, atom2idx[atom_symbol]] = 1.0
                else:
                    print(str(i) + 'find new atom type!')
                    raise ValueError('unknown atom type encountered!')
        except ValueError:
            continue
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
            'atom_one_hot': one_hot,
            'conjugate_idx': conjugate_idx,
            'playload_idx': payload_idx,
            'linker_idx': linker_idx,
            'linker_size': len(linker_idx),

            # reserved field
            'mol': clp_mol,
            'bond': full_bonds_info,
            'conjugate_mol': conjugate_mol,
            'payload_mol': payload_mol,
            'linker_mol': linker_mol,
            'anchors': anchors.bool(),
        })
        #print('n fail: ', n_fail)
    print('saving data')
    with open(save_path, 'wb') as f:
        pickle.dump(dataset, f)
    print('length processed data: ', len(dataset))


def process_zinc_from_difflinker(raw_file_dir, save_path, mode):
    datasets = {}
    if mode == 'full':
        all_subsets = ['train', 'val', 'test']
    else:
        all_subsets = ['val', 'test']
    for subset in all_subsets:
        data_list = []
        n_fail = 0

        # load data
        all_data = torch.load(os.path.join(raw_file_dir, f'zinc_final_{subset}.pt'), map_location='cpu')
        full_rdmols = Chem.SDMolSupplier(os.path.join(raw_file_dir, f'zinc_final_{subset}_mol.sdf'), sanitize=False)
        frag_rdmols = Chem.SDMolSupplier(os.path.join(raw_file_dir, f'zinc_final_{subset}_frag.sdf'), sanitize=False)
        link_rdmols = Chem.SDMolSupplier(os.path.join(raw_file_dir, f'zinc_final_{subset}_link.sdf'), sanitize=False)
        assert len(all_data) == len(full_rdmols) == len(frag_rdmols) == len(link_rdmols)
        for i in tqdm(range(len(all_data)), desc=subset):
            data = all_data[i]  # extract data
            mol, frag_mol, link_mol = full_rdmols[i], frag_rdmols[i], link_rdmols[i]  # extract link_mol
            if mol is None or frag_mol is None or link_mol is None:
                 print('Fail i: ', i)
                 n_fail += 1
                 continue
            pos = data['positions']
            fragment_mask, linker_mask = data['fragment_mask'].bool(), data['linker_mask'].bool()
            linker_idx = np.where(linker_mask)[0]  # extract linker_idx
            linker_size = np.sum(np.array(data['linker_mask']))  # extract linker_size

            # align full mol with positions, etc.
            mol_pos = mol.GetConformer().GetPositions()
            mapping = np.linalg.norm(mol_pos[None] - data['positions'].numpy()[:, None], axis=-1).argmin(axis=1)
            assert len(np.unique(mapping)) == len(mapping)
            new_mol = Chem.RenumberAtoms(mol, mapping.tolist())  # extract mol
            Chem.SanitizeMol(new_mol)
            # check frag mol and link mol are aligned
            assert np.allclose(frag_mol.GetConformer().GetPositions(), pos[fragment_mask].numpy())
            assert np.allclose(link_mol.GetConformer().GetPositions(), pos[linker_mask].numpy())
            #print(mapping)
            # Note: anchor atom index may be wrong!
            #print(data['anchors'].nonzero(as_tuple=True)[0].tolist())

            # construct new_fragment_mask
            frag_mols = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=False)
            assert len(frag_mols) == 2
            all_frag_atom_idx = set((fragment_mask == 1).nonzero()[:, 0].tolist())
            frag_atom_idx_list = []
            for m1 in new_mol.GetSubstructMatches(frag_mols[0]):
                for m2 in new_mol.GetSubstructMatches(frag_mols[1]):
                    if len(set(m1).intersection(set(m2))) == 0 and set(m1).union(set(m2)) == all_frag_atom_idx:
                        frag_atom_idx_list = [m1, m2]  # extract frag1_idx and frag2_idx
                        break
            try:
                assert len(frag_atom_idx_list) == 2 and all([x is not None and len(x) > 0 for x in frag_atom_idx_list])
            except:
                print('Fail i: ', i)
                n_fail += 1
                continue
            new_fragment_mask = torch.zeros_like(fragment_mask).long()
            new_fragment_mask[list(frag_atom_idx_list[0])] = 1
            new_fragment_mask[list(frag_atom_idx_list[1])] = 2

            # extract frag mol directly from new_mol, in case the Kekulize error
            bond_ids, full_bonds_info = [], []
            for bond_idx, bond in enumerate(new_mol.GetBonds()):
                start = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                bond_type = bond.GetBondType().name
                full_bonds_info.append((start, end, bond_type))  # extract bond
                if (new_fragment_mask[start] > 0) == (new_fragment_mask[end] == 0):
                    bond_ids.append(bond_idx)
            assert len(bond_ids) == 2
            break_mol = Chem.FragmentOnBonds(new_mol, bond_ids, addDummies=False)
            frags = [f for f in Chem.GetMolFrags(break_mol, asMols=True)
                     if f.GetNumAtoms() != link_mol.GetNumAtoms()
                     or not np.allclose(f.GetConformer().GetPositions(), link_mol.GetConformer().GetPositions())]
            assert len(frags) == 2
            new_frag_mol = Chem.CombineMols(*frags)  #extract frag_mol
            assert np.allclose(new_frag_mol.GetConformer().GetPositions(), frag_mol.GetConformer().GetPositions())

            data_list.append({
                'id': data['uuid'],
                'smiles': data['name'],
                'pos': data['positions'],
                'atom_one_hot': data['one_hot'],
                'frag1_idx': list(frag_atom_idx_list[0]),
                'frag2_idx': list(frag_atom_idx_list[1]),
                'linker_idx': linker_idx,
                'linker_size': linker_size,

                # reserved field
                'bond': full_bonds_info,
                'mol': new_mol,
                'frag_mol': new_frag_mol,
                'link_mol': link_mol,
                'anchors': data['anchors'].bool()
            })
        print('n fail: ', n_fail)
        datasets[subset] = data_list

    print('saving data')
    with open(save_path, 'wb') as f:
        pickle.dump(datasets, f)
    print('length processed data: ', [f'{x}: {len(datasets[x])}' for x in datasets])


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_path', type=str, required=True)
    parser.add_argument('--dataset', type=str, default='protac', choices=['protac', 'adc', 'zinc_difflinker' ])
    parser.add_argument('--dest', type=str, required=True)
    parser.add_argument('--mode', type=str, default='full', choices=['tiny', 'full'])  # for zinc only
    args = parser.parse_args()

    # example for run: python process_rawdata.py --raw_path zincFD/ --dest IndexFile/zinc.pkl
    if args.dataset == 'protac':
        process_protac(args.raw_path, args.dest)
    elif args.dataset == 'adc':
        process_adc(args.raw_path, args.dest)
    elif args.dataset == 'zinc_difflinker':
        process_zinc_from_difflinker(args.raw_path, args.dest, args.mode)
    else:
        raise NotImplementedError

    '''
    file_path = '/data1/private/zhangli/pj2/zenodo/protac.pkl'
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    print(data['train']['protac'][1524])
    for i in range(len(data['train']['protac'])):
        #print(i, '\n')
        num1 = len(data['train']['protac'][i]['frag1_idx'])+len(data['train']['protac'][i]['link_idx'])+len(data['train']['protac'][i]['frag2_idx'])
        num2 = len(data['train']['protac'][i]['atom_one_hot'])
        if num1 != num2:
            print('num1, num2: ', num1, num2 )
    '''
    '''
    file_path = '/data1/private/zhangli/pj2/zenodo/tarDrug.h5'
    with h5py.File(file_path, 'r') as h5_file:
        tarDrug_train,  tarDrug_val, tarDrug_test = h5_file['train'], h5_file['val'], h5_file['test']
        train_protac, val_protac,  test_protac = tarDrug_train['protac'], tarDrug_val['protac'] , tarDrug_test['protac']
        train_adc, val_adc,  test_adc = tarDrug_train['adc'], tarDrug_val['adc'] , tarDrug_test['adc']
        print("protac:", len(train_protac), len(val_protac), len(test_protac))
        print("adc:", len(train_adc), len(val_adc), len(test_adc))
        print(test_adc.keys())
        print(type(train_adc['11']))
        print(train_adc['11']['anchors'])
    '''