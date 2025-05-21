import os
import torch
import numpy as np

import argparse
import h5py
import pandas as pd
from tqdm.auto import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')  
rdBase.DisableLog('rdApp.warning')  

########### define your constants ###########
# atom types in zinc
atoms_zc = ['C', 'O', 'N', 'F', 'S', 'Cl', 'Br', 'I']
# atom types in HBDrug3D
#157->B, 185->Si, 240->Pt, 253->Si, Na (removed)
atoms_td = ['C', 'O', 'N', 'F', 'S', 'Cl', 'Br', 'I', 'P']
seed = 42  # refer to difflinker


def remove_all_hydrogens(mol):
    rwmol = Chem.RWMol(mol)
    for i in reversed(range(rwmol.GetNumAtoms())):
        if rwmol.GetAtomWithIdx(i).GetAtomicNum() == 1: 
            rwmol.RemoveAtom(i)
    return rwmol.GetMol()

def remove_all_hydrogens_with_charge(mol):
    rw_mol = Chem.RWMol(mol)
    for atom in sorted(mol.GetAtoms(), key=lambda a: a.GetIdx(), reverse=True):
        if atom.GetAtomicNum() == 1: 
            rw_mol.RemoveAtom(atom.GetIdx())  
    new_mol = rw_mol.GetMol()
    for new_atom, old_atom in zip(new_mol.GetAtoms(), mol.GetAtoms()):
        if old_atom.GetAtomicNum() != 1:  
            new_atom.SetFormalCharge(old_atom.GetFormalCharge())
    Chem.SanitizeMol(new_mol)
    return new_mol

def align_sdf_to_smiles_with_charge(sdf_path):

    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    aligned_mols = []
    
    for orig_mol in supplier:
        if not orig_mol:
            continue
            
        orig_mol_noH = remove_all_hydrogens_with_charge(orig_mol)

        for atom in orig_mol_noH.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        
        smiles = Chem.MolToSmiles(orig_mol_noH, 
                                canonical=True,
                                isomericSmiles=True,
                                allHsExplicit=False)
        
        new_mol = Chem.MolFromSmiles(smiles)
        
        # GetNumConformers() <- HasConformers()
        if orig_mol.GetNumConformers() > 0:
            orig_conf = orig_mol.GetConformer()
            new_conf = Chem.Conformer(new_mol.GetNumAtoms())
            
            # orig_noH map orig
            orig_to_origNoH = [i for i, atom in enumerate(orig_mol.GetAtoms()) 
                             if atom.GetAtomicNum() != 1]
            
            match = []
            for new_idx, atom in enumerate(new_mol.GetAtoms()):
                orig_noH_idx = atom.GetAtomMapNum() - 1
                match.append((orig_noH_idx, new_idx))
            
            for orig_noH_idx, new_idx in match:
                orig_idx = orig_to_origNoH[orig_noH_idx]
                pos = orig_conf.GetAtomPosition(orig_idx)
                new_conf.SetAtomPosition(new_idx, pos)
            
            new_mol.AddConformer(new_conf)
        
        for atom in new_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        
        aligned_mols.append(new_mol)
    
    return aligned_mols

def trans_to_MolAnd3D(sdf_path, protac_name):
    file_path = os.path.join(sdf_path, f"{protac_name}.sdf")
    if not os.path.isfile(file_path):
        print(f"Files not exist: {file_path}")
        return [] 
    supplier = Chem.SDMolSupplier(file_path)
    sdf_path_new = sdf_path + '_new'
    os.makedirs(sdf_path_new, exist_ok=True)
    output_path = os.path.join(sdf_path_new, f"{protac_name}.sdf")
    mols = [mol for mol in supplier if mol is not None]
    with Chem.SDWriter(output_path) as writer:
        for mol in supplier:
                writer.write(mol)
    return mols    


def process_protac(raw_file_dir, save_path, sdf_path):
    dataset = []
    n_fail = 0  # record fail num
    all_data = pd.read_excel(raw_file_dir, header=0)
    index = 0

    for i in tqdm(range(len(all_data))):
        id, protac_smi, warhead_smi, linker_smi, e3lig_smi = (all_data.iloc[i]['PROTAC ID'], all_data.iloc[i]['PROTAC Smiles'], all_data.iloc[i]['Warhead Smiles'], all_data.iloc[i]['Linker Smiles'], all_data.iloc[i]['E3 Ligand Smiles'])
        if protac_smi is None or warhead_smi is None or linker_smi is None or e3lig_smi is None:
            print(f'fail i: {i}')
            continue
        # get 3d coords
        try:
            protac_mols = trans_to_MolAnd3D(sdf_path, id)
        except ValueError as e:
            print(f'protac_id: {id}, error: {e}')
            continue
        try:
            linker_mol = Chem.MolFromSmiles(linker_smi)
            if linker_mol.GetNumAtoms() < 5:
                continue
        except:
            continue
        
        for protac_mol in protac_mols:
            atoms_pos_in_protac = protac_mol.GetConformer().GetPositions()
            
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
            #warhead_mol = Chem.PathToSubmol(protac_mol, warhead_idx)
            #linker_mol = Chem.PathToSubmol(protac_mol, linker_idx)
            #e3lig_mol = Chem.PathToSubmol(protac_mol, e3lig_idx)

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
            output_path = sdf_path + '_datasets/' + str(index) +'_.sdf'
            with Chem.SDWriter(output_path) as writer:
                writer.write(protac_mol)
            index += 1
        #print('n fail: ', n_fail)
        print(str(id) + 'succeed !')
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


def process_adc(raw_file_dir, save_path, sdf_path):
    dataset = []
    n_fail = 0  # record fail num
    all_data = pd.read_excel(raw_file_dir, header=0)
    os.mkdir(sdf_path + '_datasets/')
    index = 0
    
    for i in tqdm(range(len(all_data))):
        clp_name, clp_smi, conjugate_smi, linker_smi, payload_smi = (all_data.iloc[i]['Conjugate-Linker-Payload Name'], all_data.iloc[i]['Conjugate-Linker-Payload Isosmiles'], all_data.iloc[i]['Conjugate Smiles'], all_data.iloc[i]['Post-reaction Linker Isosmiles'], all_data.iloc[i]['Post-reaction Payload Isosmiles'])
        if clp_smi is None or conjugate_smi is None or linker_smi is None or payload_smi is None:
            print(f'fail i: {i}')
            continue
        linker_mol = Chem.MolFromSmiles(linker_smi)
        if linker_mol.GetNumAtoms() < 5:
            continue
        # get 3d coords
        try:
            clp_mols = trans_to_MolAnd3D(sdf_path, clp_name)
        except ValueError as e:
            print(e)
        if len(clp_mols) == 0:
            continue
        
        for clp_mol in clp_mols:
            atoms_pos_in_clp = clp_mol.GetConformer().GetPositions()
            
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
                'clp_name': clp_name,
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
    
            output_path = sdf_path + '_datasets/' + str(index) +'_.sdf'
            with Chem.SDWriter(output_path) as writer:
                writer.write(clp_mol)
            index += 1
            print(str(i) + "succeed !")
            #print('n fail: ', n_fail)
    
    with h5py.File(save_path, 'w') as f:
        grp = f.create_group('adc')
    
        for idx, data in enumerate(dataset):
            sample_grp = grp.create_group(f'{idx}')
            sample_grp.attrs['clp_name'] = data['clp_name']
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


def process_pdc(raw_file_dir, save_path, sdf_path):
    dataset = []
    data_record = []
    n_fail = 0  # record fail num
    all_data = pd.read_excel(raw_file_dir, header=0)
    index = 0
    os.mkdir(sdf_path + '_datasets/')
    
    for i in tqdm(range(len(all_data))):
        pdc_name, pdc_smi, peptide_smi, linker_smi, payload_smi = (all_data.iloc[i]['PDC ID'], all_data.iloc[i]['PDC Smiles'], all_data.iloc[i]['Peptide Smiles'], all_data.iloc[i]['Linker Smiles'], all_data.iloc[i]['Payload Smiles'])
        if pdc_smi is None or peptide_smi is None or linker_smi is None or payload_smi is None:
            print(f'fail i: {i}')
            continue
        linker_mol = Chem.MolFromSmiles(linker_smi)
        if linker_mol.GetNumAtoms() < 5:
            continue
        # get 3d coords
        try:
            pdc_mols = trans_to_MolAnd3D(sdf_path, pdc_name)
            #pdc_mols = trans_to_MolAnd3D(pdc_smi)
        except ValueError as e:
            print(e)
        if len(pdc_mols) == 0:
            continue
        
        for pdc_mol in pdc_mols:
            atoms_pos_in_pdc = pdc_mol.GetConformer().GetPositions()
            
            # get atom types
            atoms_num = pdc_mol.GetNumAtoms()
            atomic_numbers = torch.zeros(atoms_num, dtype=torch.long)  
            for j in range(atoms_num):
                atom = pdc_mol.GetAtomWithIdx(j)
                atomic_numbers[j] = atom.GetAtomicNum()
            
            # get bond_info
            full_bonds_info = []
            for bond in pdc_mol.GetBonds():
                start = bond.GetBeginAtomIdx()
                end = bond.GetEndAtomIdx()
                bond_type = bond.GetBondType().name
                full_bonds_info.append((start, end, bond_type))

            # get fragment and linker indices
            #peptide_idxs = pdc_mol.GetSubstructMatches(Chem.MolFromSmiles(peptide_smi))
            try:
                linker_idxs = pdc_mol.GetSubstructMatches(Chem.MolFromSmiles(linker_smi))
                payload_idxs = pdc_mol.GetSubstructMatches(Chem.MolFromSmiles(payload_smi))
            except:
                continue
            if len(linker_idxs) == 0:
                print(str(i) + 'linker_idxs error !')
            if len(payload_idxs) == 0:
                print(str(i) + 'payload_idxs error !')
            flag = False
            pdc_atoms_idx = set(range(pdc_mol.GetNumAtoms()))
            for linker_idx in linker_idxs:
                for payload_idx in payload_idxs:
                    if (len(set(linker_idx).intersection(set(payload_idx)))!=0):
                        continue
                    else:
                        peptide_idx = pdc_atoms_idx - set(linker_idx + payload_idx)
                        peptide_mol = Chem.RWMol(pdc_mol)
                        for atom in reversed(range(peptide_mol.GetNumAtoms())):
                            if atom not in peptide_idx:
                                peptide_mol.RemoveAtom(atom)
                        peptide_mol = peptide_mol.GetMol()
                        try:
                            frags = Chem.GetMolFrags(peptide_mol, asMols=True)
                        except:
                            continue
                        if len(frags) == 1:
                            flag = True
                    if flag:
                        break
                if flag:
                    break
            if not flag:
                print(str(i) + 'validity check fail!')
                continue

            # get anchors
            peptide_idx, linker_idx, payload_idx = list(peptide_idx), list(linker_idx), list(payload_idx)
            anchors = torch.zeros(atoms_num, dtype=int)
            frag_mask = torch.zeros(atoms_num, dtype=int)
            frag_mask[peptide_idx] = 1
            frag_mask[payload_idx] = 2
            for bond in pdc_mol.GetBonds():
                start_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if frag_mask[start_idx] != frag_mask[end_idx]:
                    anchors[start_idx] = 1
                    anchors[end_idx] = 1
            anchors *= frag_mask
            
            # get submol
            #peptide_mol = Chem.PathToSubmol(pdc_mol, peptide_idx)
            #linker_mol = Chem.PathToSubmol(pdc_mol, linker_idx)
            #payload_mol = Chem.PathToSubmol(pdc_mol, payload_idx)
            #linker_mol = Chem.MolFromSmiles(linker_smi)
            #payload_mol = Chem.MolFromSmiles(payload_smi)
                
                dataset.append({
                    'pdc_name': pdc_name,
                    'smiles': [pdc_smi, peptide_smi, linker_smi, payload_smi],
                    'pos': atoms_pos_in_pdc,
                    'atoms_atomic_numbers': atomic_numbers,
                    'frag1_idx': peptide_idx,
                    'frag2_idx': payload_idx,
                    'linker_idx': linker_idx,
                    'linker_size': len(linker_idx),
                    'bond': full_bonds_info,
                    'anchors': anchors.bool(),
                })
                
                output_path = sdf_path + '_datasets/' + str(index) +'_.sdf'
                with Chem.SDWriter(output_path) as writer:
                    writer.write(pdc_mol)
                index += 1
                print(str(i) + "succeed !")
            #print('n fail: ', n_fail)
            except:
                continue
    
    with h5py.File(save_path, 'w') as f:
        grp = f.create_group('pdc')
    
        for idx, data in enumerate(dataset):
            sample_grp = grp.create_group(f'{idx}')
            sample_grp.attrs['pdc_name'] = data['pdc_name']
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
    parser.add_argument('--dataset', type=str, default='protac', choices=['protac', 'adc', 'pdc'])
    parser.add_argument('--dest', type=str, required=True)
    parser.add_argument('--sdf_path', type=str)
    parser.add_argument('--mode', type=str, default='full', choices=['tiny', 'full'])  # for zinc only
    args = parser.parse_args()

    # example for run: python process_rawdata.py --raw_path /protac.xlsx --dest /protac.h5 --dataset protac --sdf_path /protac_sdf
    if args.dataset == 'protac':
        process_protac(args.raw_path, args.dest, args.sdf_path)
    elif args.dataset == 'adc':
        process_adc(args.raw_path, args.dest, args.sdf_path)
    elif args.dataset == 'pdc':
        process_pdc(args.raw_path, args.dest, args.sdf_path)
    else:
        raise NotImplementedError