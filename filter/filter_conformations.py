import os
import torch
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem


def init_gpu():
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA GPU not available")
    device = torch.device('cuda')
    print(f"Using GPU: {torch.cuda.get_device_name(0)}")
    return device


def calculate_rmsd(mol1, mol2, device):
    aligner = AllChem.AlignMol(mol1, mol2, includeChirality=False)
    if not aligner:
        return float('inf')
    
    mol1_aligned = Chem.Mol(mol1)
    mol2_aligned = Chem.Mol(mol2)
    AllChem.TransformMol(mol1_aligned, aligner[1])

    conf1 = mol1_aligned.GetConformer()
    conf2 = mol2_aligned.GetConformer()
    conf1 = mol1.GetConformer()
    conf2 = mol2.GetConformer()
    atoms1 = torch.tensor([conf1.GetAtomPosition(i) for i in range(mol1.GetNumAtoms())],
                        dtype=torch.float32, device=device)
    atoms2 = torch.tensor([conf2.GetAtomPosition(i) for i in range(mol2.GetNumAtoms())],
                        dtype=torch.float32, device=device)

    delta = atoms1 - atoms2
    sq_diff = delta ** 2
    rmsd = torch.sqrt(sq_diff.sum(dim=0)).mean().item()
    return rmsd


def process_single_file(input_path, output_path, threshold, device):
    supplier = Chem.SDMolSupplier(input_path)
    mols = [mol for mol in supplier if mol is not None]
    
    if not mols:
        print(f"Skip Empty Files: {input_path}")
        return None
    if len(mols) == 1:
        print(f"No 3D Structures: {input_path}")
        return None 
    mols = mols[1:]
    keep_indices = []
    reference_mols = []

    for i, ref_mol in enumerate(tqdm(mols, desc="Processing conformers")):
        is_duplicate = False
        
        for j, kept_mol in enumerate(reference_mols):
            rmsd = calculate_rmsd(ref_mol, kept_mol, device)
            if rmsd < threshold:
                is_duplicate = True
                break
        
        if not is_duplicate:
            keep_indices.append(i)
            reference_mols.append(ref_mol)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    writer = Chem.SDWriter(output_path)
    for idx in keep_indices:
        writer.write(mols[idx])
    writer.close()

    print(f"Process Finished: {input_path} -> Reserve {len(keep_indices)} conformations")
    return keep_indices


def process_directory(input_dir, output_dir, threshold=0.5, use_gpu=True):
    device = init_gpu() if use_gpu else torch.device('cpu')
    
    results = []
    for filename in tqdm(os.listdir(input_dir), desc="Processing files"):
        if filename.endswith(".sdf"):
            input_path = os.path.join(input_dir, filename)
            output_path = os.path.join(output_dir, f"{filename}")
            
            try:
                result = process_single_file(input_path, output_path, threshold, device)
                if result is not None:
                    results.append(result)
            except Exception as e:
                print(f"Processing {filename} Error: {str(e)}")
    
    return results


if __name__ == "__main__":

    process_directory(
        input_dir="./protac_conformations_sdf",
        output_dir="./protac_conformations_sdf_filtered",
        threshold=0.5,
        use_gpu=True
    )