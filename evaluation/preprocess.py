import h5py
import numpy as np
import torch
from tqdm import tqdm
import os
import argparse
import pickle
import json
from rdkit import Chem
from itertools import product
from rdkit.Chem import AllChem, rdmolops

import sys
sys.path.append('..')
from utils.constants import *


# bond mapping
bond_dict = {'SINGLE': 0, 'DOUBLE': 1, 'TRIPLE': 2, 'AROMATIC': 3}
number_to_bond= {0: Chem.rdchem.BondType.SINGLE, 1: Chem.rdchem.BondType.DOUBLE,
                 2: Chem.rdchem.BondType.TRIPLE, 3: Chem.rdchem.BondType.AROMATIC}


def get_one_hot(charge):
    one_hot = np.zeros(DRUG_NUMBER_OF_ATOM_TYPES)
    one_hot[DRUG_ATOM2IDX[DRUG_CHARGES2ATOM[charge]]] = 1
    return one_hot


def get_mask(num_atoms, *mask_index):
    mask = np.zeros(num_atoms, dtype=float)
    for index in mask_index:
        mask[index] = 1
    return mask


def extract_mol_from_smiles(frag_idx, mol):

    frag_mol = Chem.RWMol()
    for atom_idx in frag_idx:
        atom = mol.GetAtomWithIdx(atom_idx)
        frag_mol.AddAtom(atom)

    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in frag_idx and bond.GetEndAtomIdx() in frag_idx:
            begin_idx = frag_idx.index(bond.GetBeginAtomIdx())
            end_idx = frag_idx.index(bond.GetEndAtomIdx())
            frag_mol.AddBond(begin_idx, end_idx, bond.GetBondType())

    frag_mol = frag_mol.GetMol()
    return frag_mol


def align_mol_to_frags(smi_molecule, smi_linker, smi_frags):
    # Load SMILES as molecules
    mol = Chem.MolFromSmiles(smi_molecule)
    if not (frags := Chem.MolFromSmiles(smi_frags)):
        frags = Chem.MolFromSmiles(smi_frags, sanitize=False)
    if not (linker := Chem.MolFromSmiles(smi_linker)):
        linker = Chem.MolFromSmiles(smi_linker, sanitize=False)

    # Renumber molecule based on frags (incl. dummy atoms)
    aligned_mols = []

    sub_idx = []
    # Get matches to fragments and linker
    frags_matches = list(mol.GetSubstructMatches(frags))
    linker_matches = list(mol.GetSubstructMatches(linker))

    # Loop over matches
    for frag_match, linker_match in product(frags_matches, linker_matches):
        if len(set(frag_match + linker_match)) == len(frag_match) + len(linker_match):
            break
    # Add frag indices
    sub_idx += frag_match
    # Add linker indices to end
    sub_idx += [idx for num, idx in enumerate(linker_match) if
                linker.GetAtomWithIdx(num).GetAtomicNum() != 0 and idx not in sub_idx]

    aligned_mols.append(Chem.rdmolops.RenumberAtoms(mol, sub_idx))
    aligned_mols.append(frags)
    nodes_to_keep = [i for i in range(len(frag_match))]

    return (aligned_mols[0], aligned_mols[1]), (list(frag_match), list(linker_match)), nodes_to_keep, sub_idx


def need_kekulize(mol):
    for bond in mol.GetBonds():
        if bond_dict[str(bond.GetBondType())] >= 3:
            return True
    return False


def to_graph_mol(mol):
    if mol is None:
        return [], []
    # Kekulize it
    if need_kekulize(mol):
        rdmolops.Kekulize(mol)
        if mol is None:
            return None, None
    # remove stereo information, such as inward and outward edges
    Chem.RemoveStereochemistry(mol)

    edges = []
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        begin_idx, end_idx = min(begin_idx, end_idx), max(begin_idx, end_idx)
        if mol.GetAtomWithIdx(begin_idx).GetAtomicNum() == 0 or mol.GetAtomWithIdx(end_idx).GetAtomicNum() == 0:
            continue
        else:
            edges.append((begin_idx, bond_dict[str(bond.GetBondType())], end_idx))
            assert bond_dict[str(bond.GetBondType())] != 3

    return edges


def get_submol_by_indices(mol, atom_indices):
    rw_mol = Chem.RWMol()
    
    idx_map = {}
    
    for orig_idx in atom_indices:
        atom = mol.GetAtomWithIdx(orig_idx)
        new_idx = rw_mol.AddAtom(atom)
        idx_map[orig_idx] = new_idx
    
    for bond in mol.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        if begin in idx_map and end in idx_map:
            rw_mol.AddBond(
                idx_map[begin],
                idx_map[end],
                bond.GetBondType()
            )
    
    if mol.GetNumConformers() > 0:
        conf = mol.GetConformer()
        new_conf = Chem.Conformer(rw_mol.GetNumAtoms())
        for orig_idx, new_idx in idx_map.items():
            pos = conf.GetAtomPosition(orig_idx)
            new_conf.SetAtomPosition(new_idx, pos)
        rw_mol.AddConformer(new_conf)
    
    return rw_mol.GetMol()


def mol_to_sdf(data_idx_list, linker_idx_list, frags_idx_list, sdf_path, save_dir, split):
    mol_fp = open(os.path.join(save_dir, f'tarDrug_{split}_molecules.sdf'), 'w')
    frag_fp = open(os.path.join(save_dir, f'tarDrug_{split}_frag.sdf'), 'w')
    link_fp = open(os.path.join(save_dir, f'tarDrug_{split}_link.sdf'), 'w')
    
    for i in range(len(data_idx_list)):
        mol = Chem.SDMolSupplier(os.path.join(sdf_path, str(data_idx_list[i])+'_.sdf'))
        mol = mol[0]
        frag_mol = get_submol_by_indices(mol, frags_idx_list[i])
        linker_mol = get_submol_by_indices(mol, linker_idx_list[i])
        writer = Chem.SDWriter(mol_fp)
        writer.SetKekulize(False)
        writer.write(mol)
        writer = Chem.SDWriter(frag_fp)
        writer.SetKekulize(False)
        writer.write(frag_mol)   
        writer = Chem.SDWriter(link_fp)
        writer.SetKekulize(False)
        writer.write(linker_mol)
    

def read_indices_file(ids_filepath):
    with open(ids_filepath, "r", encoding="utf-8") as f:
        indices = [int(line.strip()) for line in f if line.strip()]
    return indices


def preprocess_to_pt(raw_data_path, ids_txt_path, save_dir, sdf_path, split="train"):
    data = []
    linker_smiles, mol_smiles, frag_smiles, mol_pos, mol_charges = [], [], [], [], []

    indices = read_indices_file(ids_txt_path)
    data_idx_list = []
    linker_idx_list = []
    frags_idx_list = []
    with h5py.File(raw_data_path, 'r') as f:
        
        all_items = f["pdc"] # revise according to data type
        uuid = 0
        for idx in tqdm(indices, desc=f"Processing {split} split"):            
            key = str(idx)
            try:
                item = all_items[key]
            except Exception as e:
                print(f"{idx} reading data failure: {e}")
                continue
            name = item["smiles"][:][0]
            frag1_idx, frag2_idx, linker_idx = item["frag1_idx"][:], item["frag2_idx"][:], item["linker_idx"][:]
            anchors = item["anchors"][:]
            pos = item["pos"][:]
            charges = item["atoms_atomic_numbers"][:]
            num_atoms = len(charges)

            # transfer atomic numbers to one-hot representation
            one_hot = np.array([get_one_hot(charge) for charge in charges], dtype=float)

            fragment_mask = get_mask(num_atoms, frag1_idx, frag2_idx)
            linker_mask = get_mask(num_atoms, linker_idx)
            frag1_mask = get_mask(num_atoms, frag1_idx)
            frag2_mask = get_mask(num_atoms, frag2_idx)

            # CoM
            fragment_idx = np.where(fragment_mask > 0)[0]
            pos = pos - pos[fragment_idx].mean(axis=0)

            data.append({
                "uuid": uuid,
                "name": name,
                "positions": torch.tensor(pos, dtype=torch.float),
                "one_hot": torch.tensor(one_hot, dtype=torch.float),
                "charges": torch.tensor(charges, dtype=torch.float),
                "anchors": torch.tensor(anchors, dtype=torch.float),
                "fragment_mask": torch.tensor(fragment_mask, dtype=torch.float),
                "linker_mask": torch.tensor(linker_mask, dtype=torch.float),
                #"frag1_mask": torch.tensor(frag1_mask, dtype=torch.float),
                #"frag2_mask": torch.tensor(frag2_mask, dtype=torch.float),
                "num_atoms": num_atoms
            })

            mol_smi, frag1_smi, linker_smi, frag2_smi = item["smiles"][:]
            mol_smiles.append(mol_smi.decode('utf-8'))
            linker_smiles.append(linker_smi.decode('utf-8'))
            frag_smiles.append(f'{frag1_smi.decode("utf-8")}.{frag2_smi.decode("utf-8")}')
            mol_pos.append(pos)
            mol_charges.append(charges)

            uuid += 1
            data_idx_list.append(idx)
            linker_idx = list(linker_idx)
            linker_idx = [int(idx) for idx in linker_idx]
            linker_idx_list.append(linker_idx)
            frags_idx = list(frag1_idx)+list(frag2_idx)
            frags_idx = [int(idx) for idx in frags_idx]
            frags_idx_list.append(frags_idx)
            print(uuid)

    torch.save(data, os.path.join(save_dir, f'{split}.pt'))

    # write molecular structures to .sdf file
    mol_to_sdf(data_idx_list, linker_idx_list, frags_idx_list, sdf_path, save_dir, split)
    

    with open(os.path.join(save_dir, f'molecules_tarDrug_{split}_final.smi'), 'w') as f:
        # Format: full_mol (SMILES), linker (SMILES), fragments (SMILES), distance (Angstrom), angle (Radians)
        for mol_smi, linker_smi, frag_smi in zip(mol_smiles, linker_smiles, frag_smiles):
            # linker_smi_with_anchor, frags_smi_with_anchor = mark_anchors_by_indices(mol, linker, frag, anchor_idx)
            f.write("%s %s %s %s %s\n" % (mol_smi, linker_smi, frag_smi, '0.0', '0.0'))
    
    processed_data = []
    n_exceed, max_atom_num = 0, 200 # revise according to data type
    for data_idx, (mol_smi, linker_smi, frag_smi, d) in enumerate(zip(mol_smiles, linker_smiles, frag_smiles, data)):
        (mol_out, mol_in), (frag_idx, link_idx), nodes_to_keep, sub_idx = align_mol_to_frags(mol_smi, linker_smi, frag_smi)
        assert len(sub_idx) == len(d['positions']), 'molecule to graph failed'
        if len(sub_idx) >= max_atom_num:
            print(f"split {split}, index {data_idx} out of max atom num: {len(sub_idx)} > {max_atom_num}")
            n_exceed += 1
            continue
        exit_points = np.where(d['anchors'] == 1)[0].tolist()
        if len(exit_points) != 2:
            continue
        # re-indexing
        exit_points = [sub_idx.index(i) for i in exit_points]
        pos_out = d['positions'].numpy()[sub_idx]
        nodes_out = d['one_hot'].numpy().astype(np.int64)[sub_idx]
        pos_in = pos_out[:len(frag_idx)]
        nodes_in = nodes_out[:len(frag_idx)]
        edges_out = to_graph_mol(mol_out)
        # edges_in = to_graph_mol(mol_in)
        edges_in = []
        for edge in edges_out:
            begin_idx, _, end_idx = edge
            if begin_idx in frag_idx and end_idx in frag_idx:
                edges_in.append([frag_idx.index(begin_idx), edge[1], frag_idx.index(end_idx)])
        processed_data.append({
            'graph_in': edges_in,
            'graph_out': edges_out,
            'node_features_in': nodes_in.tolist(),
            'node_features_out': nodes_out.tolist(),
            'smiles_out': mol_smi,
            'smiles_in': frag_smi,
            'v_to_keep': nodes_to_keep,
            'exit_points': exit_points,
            'abs_dist': ["0.00", "0.00"],
            'positions_out': pos_out.tolist(),
            'positions_in': pos_in.tolist()
        })
    print(f"number of exceed mol: {n_exceed}")
    with open(os.path.join(save_dir, f'molecules_tarDrug_{split}_final.json'), 'w') as f:
        json.dump(processed_data, f)
    if split == 'train':
        with open(os.path.join(save_dir, 'tarDrug_train_linkers.smi'), 'w') as f:
            for linker_smi in linker_smiles:
                f.write(f"{linker_smi}\n")
    elif split == 'test':
        with open(os.path.join(save_dir, 'tarDrug_test_smiles.smi'), 'w') as f:
            for mol_smi, frag_smi in zip(mol_smiles, frag_smiles):
                f.write(f"{mol_smi} {frag_smi}\n")


def extract_atoms_with_mask(mol, mask):
    assert len(mask) == mol.GetNumAtoms(), \
        f"Mask length {len(mask)} does not match the number of atoms {mol.GetNumAtoms()}"

    editable_mol = Chem.EditableMol(Chem.Mol())
    atom_map = {}

    for i, atom in enumerate(mol.GetAtoms()):
        if mask[i] == 1:
            new_idx = editable_mol.AddAtom(atom)
            atom_map[i] = new_idx

    for bond in mol.GetBonds():
        start_idx, end_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if start_idx in atom_map and end_idx in atom_map:
            editable_mol.AddBond(
                atom_map[start_idx],
                atom_map[end_idx],
                bond.GetBondType()
            )

    new_mol = editable_mol.GetMol()

    for conf in mol.GetConformers():                                                                
        new_conf = Chem.Conformer(len(atom_map))                                                    
        for old_idx, new_idx in atom_map.items():                                                   
            pos = conf.GetAtomPosition(old_idx)                                                     
            new_conf.SetAtomPosition(new_idx, pos)                                                  
        new_mol.AddConformer(new_conf, assignId=True)
    
    return new_mol


def preprocess_to_pkl(raw_data_dir, save_path):
    datasets = {}
    all_subsets = ["train", "val", "test"]
    for subset in all_subsets:
        data_list = []
        all_data = torch.load(os.path.join(raw_data_dir, f'{subset}.pt'), map_location='cpu')
        full_rdmols = Chem.SDMolSupplier(os.path.join(raw_data_dir, f'tarDrug_{subset}_molecules.sdf'), sanitize=False)
        assert len(all_data) == len(full_rdmols)
        for i in tqdm(range(len(all_data)), desc=subset):
            data = all_data[i]
            mol = full_rdmols[i]
            pos = data['positions']
            # check positions are aligned
            # assert np.allclose(mol.GetConformer().GetPositions(), pos.numpy(), rtol=1e-3)

            fragment_mask, linker_mask = data['fragment_mask'].bool(), data['linker_mask'].bool()
            frag1_mask, frag2_mask = data['frag1_mask'].bool(), data['frag2_mask'].bool()

            new_fragment_mask = torch.zeros_like(fragment_mask).long()
            new_fragment_mask[frag1_mask] = 1
            new_fragment_mask[frag2_mask] = 2

            frag_mol = extract_atoms_with_mask(mol, fragment_mask.numpy().astype(np.int64))
            #frag2_mol = extract_atoms_with_mask(mol, frag2_mask.numpy().astype(np.int64))
            #frag_mol = Chem.CombineMols(*[frag1_mol, frag2_mol])
            link_mol = extract_atoms_with_mask(mol, linker_mask.numpy().astype(np.int64))

            data_list.append({
                'id': data['uuid'],
                'smiles': data['name'],
                'mol': mol,
                'frag_mol': frag_mol,
                'link_mol': link_mol,
                'fragment_mask': new_fragment_mask,
                'atom_indices_f1': torch.where(new_fragment_mask == 1)[0].numpy().tolist(),
                'atom_indices_f2': torch.where(new_fragment_mask == 2)[0].numpy().tolist(),
                'linker_mask': linker_mask
            })
        datasets[subset] = data_list

    print('Saving data')
    with open(save_path, 'wb') as f:
        pickle.dump(datasets, f)
    print('Length processed data: ', [f'{x}: {len(datasets[x])}' for x in datasets])


def parse():
    p = argparse.ArgumentParser(description='process for train/val/test')
    p.add_argument('--raw_data_name', type=str, required=True, help='dataset name')
    p.add_argument('--raw_data_path', type=str, required=True, help='path to the raw dataset file (.h5)')
    p.add_argument('--sdf_path', type=str, required=True, help='path to the sdf file (.h5)')
    p.add_argument('--save_dir', type=str, required=True, help='path to save the processed dataset')
    return p.parse_args()


if __name__ == "__main__":
    args = parse()

    base_dir = os.path.dirname(args.raw_data_path)
    for split in ["train", "val", "test"]:
        ids_txt_path = os.path.join(base_dir, f"{args.raw_data_name}_{split}_ids.txt")
        preprocess_to_pt(args.raw_data_path, ids_txt_path, args.save_dir, args.sdf_path, split)

    raw_data_dir = os.path.split('/evaluation')[0]
    save_path = os.path.join('/evaluation', 'index_full.pkl')
    preprocess_to_pkl('/evaluation', save_path)