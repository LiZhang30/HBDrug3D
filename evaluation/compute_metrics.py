#!/usr/bin/env python
import csv
import numpy as np
import pandas as pd
import sys
import argparse
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl

from tqdm import tqdm
from rdkit import Chem
from pdb import set_trace
from networkx.algorithms import isomorphism
from utils.delinker_utils import calc_SC_RDKit, frag_utils, sascorer
from rdkit.Chem import MolStandardize, QED, rdMolAlign, rdMolDescriptors


def disable_rdkit_logging():

    logger = rkl.logger()
    logger.setLevel(rkl.ERROR)
    rkrb.DisableLog('rdApp.error')

disable_rdkit_logging()


parser = argparse.ArgumentParser()
parser.add_argument('--data_set', action='store', type=str, required=True, default='ZINC', choices=['ZINC', 'HBDrug3D'])
parser.add_argument('--data_set_path', action='store', type=str, required=True)
parser.add_argument('--gen_smi_file', action='store', type=str, required=True)
parser.add_argument('--train_set_path', action='store', type=str, required=True)
parser.add_argument('--n_cores', action='store', type=int, default=5, required=True)
parser.add_argument('--verbose', action='store_true', default=False)
parser.add_argument('--restrict', action='store', type=str, default='None')
parser.add_argument('--pains_smarts_loc', action='store', type=str, default='./resources/wehi_pains.csv')
parser.add_argument('--method', action='store', type=str, required=True, default='unilinker',
                    choices=['unilinker', 'difflinker', '3dlinker', 'delinker', 'linkernet'])

args = parser.parse_args()

data_set = args.data_set  # Options: ZINC, HBDrug3D
data_set_path = args.data_set_path
gen_smi_file = args.gen_smi_file  # Path to generated molecules
train_set_path = args.train_set_path  # Path to training set
n_cores = args.n_cores  # Number of cores to use
verbose = args.verbose  # Output results
if args.restrict == "None":
    restrict = None
else:
    restrict = int(args.resitrct)  # Set to None if don't want to restrict
pains_smarts_loc = args.pains_smarts_loc  # Path to PAINS SMARTS
method = args.method

if verbose:
    print("##### Start Settings #####")
    print("Data set:", data_set)
    print("Generated smiles file:", gen_smi_file)
    print("Training set:", train_set_path)
    print("Number of cores:", n_cores)
    print("Verbose:", verbose)
    print("Restrict data:", restrict)
    print("PAINS SMARTS location:", pains_smarts_loc)
    print("#####  End Settings  #####")


# Load molecules
# FORMAT: (Starting fragments (SMILES), Original molecule (SMILES), Generated molecule (SMILES), Generated linker)
data = []
with open(gen_smi_file, 'r') as f:
    for line in tqdm(f.readlines()):
        parts = line.strip().split(' ')
        data.append({
            'fragments': parts[0],
            'true_molecule': parts[1],
            'pred_molecule': parts[2],
            'pred_linker': parts[3] if len(parts) > 3 else '',
        })

if restrict is not None:
    data = data[:restrict]


summary = {}


# -------------- Validity -------------- #

def is_valid(pred_mol_smiles, frag_smiles):
    pred_mol = Chem.MolFromSmiles(pred_mol_smiles)
    frag = Chem.MolFromSmiles(frag_smiles)
    if frag is None:
        return False
    if pred_mol is None:
        return False
    try:
        Chem.SanitizeMol(pred_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:
        return False
    if len(pred_mol.GetSubstructMatch(frag)) != frag.GetNumAtoms():
        return False
    return True


valid_cnt = 0
total_cnt = 0
for obj in tqdm(data):
    valid = is_valid(obj['pred_molecule'], obj['fragments'])
    obj['valid'] = valid
    valid_cnt += valid
    total_cnt += 1

validity = valid_cnt / total_cnt * 100
print(f'Validity: {validity:.3f}%')
summary['validity'] = validity


# ----------------- QED ------------------ #

qed_values = []
for obj in tqdm(data):
    if not obj['valid']:
        obj['qed'] = None
        continue

    qed = QED.qed(Chem.MolFromSmiles(obj['pred_molecule']))
    obj['qed'] = qed
    qed_values.append(qed)

print(f'Mean QED: {np.mean(qed_values):.3f}')
summary['qed'] = np.mean(qed_values)


# ----------------- SA ------------------ #

sa_values = []
for obj in tqdm(data):
    if not obj['valid']:
        obj['sa'] = None
        continue

    sa = sascorer.calculateScore(Chem.MolFromSmiles(obj['pred_molecule']))
    obj['sa'] = sa
    sa_values.append(sa)

print(f'Mean SA: {np.mean(sa_values):.3f}')
summary['sa'] = np.mean(sa_values)


# ----------------- Number of Rings ------------------ #

rings_n_values = []
for obj in tqdm(data):
    if not obj['valid']:
        obj['rings_n'] = None
        continue

    try:
        rings_n = rdMolDescriptors.CalcNumRings(Chem.MolFromSmiles(obj['pred_linker']))
    except:
        continue

    obj['rings_n'] = rings_n
    rings_n_values.append(rings_n)

print(f'Mean Number of Rings: {np.mean(rings_n_values):.3f}')
summary['rings_n'] = np.mean(rings_n_values)


# -------------- Uniqueness -------------- #

true2samples = dict()
for obj in tqdm(data):
    if not obj['valid']:
        continue

    true_mol = obj['true_molecule']
    true_frags = obj['fragments']
    key = f'{true_mol}_{true_frags}'
    true2samples.setdefault(key, []).append(obj['pred_molecule'])

unique_cnt = 0
total_cnt = 0
for samples in tqdm(true2samples.values()):
    unique_cnt += len(set(samples))
    total_cnt += len(samples)

uniqueness = unique_cnt / total_cnt * 100
print(f'Uniqueness: {uniqueness:.3f}%')
summary['uniqueness'] = uniqueness


# ----------------- Novelty ---------------- #

linkers_train = set()
with open(train_set_path, 'r') as f:
    for line in f:
        linkers_train.add(line.strip())

novel_cnt = 0
total_cnt = 0
for obj in tqdm(data):
    if not obj['valid']:
        obj['pred_linker_clean'] = None
        obj['novel'] = False
        continue

    try:
        linker = Chem.RemoveStereochemistry(obj['pred_linker'])
        linker = MolStandardize.canonicalize_tautomer_smiles(Chem.MolToSmiles(linker))
    except Exception:
        linker = obj['pred_linker']

    novel = linker not in linkers_train
    obj['pred_linker_clean'] = linker
    obj['novel'] = novel
    novel_cnt += novel
    total_cnt += 1

novelty = novel_cnt / total_cnt * 100
print(f'Novelty: {novelty:.3f}%')
summary['novelty'] = novelty


# ----------------- Recovery ---------------- #

recovered_inputs = set()
all_inputs = set()
for obj in tqdm(data):
    if not obj['valid']:
        obj['recovered'] = False
        continue

    key = obj['true_molecule'] + '_' + obj['fragments']

    try:
        true_mol = Chem.MolFromSmiles(obj['true_molecule'])
        Chem.RemoveStereochemistry(true_mol)
        true_mol_smi = Chem.MolToSmiles(Chem.RemoveHs(true_mol))
    except:
        true_mol = Chem.MolFromSmiles(obj['true_molecule'], sanitize=False)
        Chem.RemoveStereochemistry(true_mol)
        true_mol_smi = Chem.MolToSmiles(Chem.RemoveHs(true_mol, sanitize=False))

    pred_mol = Chem.MolFromSmiles(obj['pred_molecule'])
    Chem.RemoveStereochemistry(pred_mol)
    pred_mol_smi = Chem.MolToSmiles(Chem.RemoveHs(pred_mol))

    recovered = true_mol_smi == pred_mol_smi
    obj['recovered'] = recovered
    if recovered:
        recovered_inputs.add(key)
    all_inputs.add(key)

recovery = len(recovered_inputs) / len(all_inputs) * 100
print(f'Recovery: {recovery:.3f}%')
summary['recovery'] = recovery


# ----------------- PAINS Filter ---------------- #

def check_pains(mol, pains):
    for pain in pains:
        if mol.HasSubstructMatch(pain):
            return False
    return True


with open(pains_smarts_loc, 'r') as f:
    pains_smarts = [Chem.MolFromSmarts(line[0], mergeHs=True) for line in csv.reader(f)]
    pains_smarts = set(pains_smarts)


passed_pains_cnt = 0
total_cnt = 0
for obj in tqdm(data):
    if not obj['valid']:
        obj['passed_pains'] = False
        continue

    pred_mol = Chem.MolFromSmiles(obj['pred_molecule'])
    passed_pains = check_pains(pred_mol, pains_smarts)
    obj['passed_pains'] = passed_pains
    passed_pains_cnt += passed_pains
    total_cnt += 1

pains_score = passed_pains_cnt / total_cnt * 100
print(f'Passed PAINS: {pains_score:.3f}%')
summary['pains'] = pains_score


# ----------------- RA Filter ---------------- #

def check_ring_filter(linker):
    check = True
    ssr = Chem.GetSymmSSSR(linker)
    for ring in ssr:
        for atom_idx in ring:
            for bond in linker.GetAtomWithIdx(atom_idx).GetBonds():
                if bond.GetBondType() == 2 and bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                    check = False
    return check


passed_ring_filter_cnt = 0
total_cnt = 0
for obj in tqdm(data):
    if not obj['valid']:
        obj['passed_ring_filter'] = False
        continue

    pred_linker = Chem.MolFromSmiles(obj['pred_linker'], sanitize=False)
    try:
        passed_ring_filter = check_ring_filter(pred_linker)
    except:
        obj['passed_ring_filter'] = False
        continue
    obj['passed_ring_filter'] = passed_ring_filter
    passed_ring_filter_cnt += passed_ring_filter
    total_cnt += 1

ra_score = passed_ring_filter_cnt / total_cnt * 100
print(f'Passed Ring Filter: {ra_score:.3f}%')
summary['ra'] = ra_score


# ---------------------------- Saving -------------------------------- #

out_path = gen_smi_file[:-3] + 'csv'
table = pd.DataFrame(data)
table.to_csv(out_path, index=False)

summary_path = gen_smi_file[:-4] + '_summary.csv'
summary_table = pd.DataFrame([summary])
summary_table.to_csv(summary_path, index=False)


# ----------------- 2D Filter ---------------- #

passed_2d_filter_cnt = 0
total_cnt = 0
for obj in tqdm(data):
    if not obj['valid']:
        obj['passed_2d_filter'] = False
        continue

    pred_linker = Chem.MolFromSmiles(obj['pred_linker'], sanitize=False)
    try:
        passed_ring_filter = check_ring_filter(pred_linker)
    except:
        obj['passed_ring_filter'] = False
        continue
    
    pred_mol = Chem.MolFromSmiles(obj['pred_molecule'])
    passed_pains = check_pains(pred_mol, pains_smarts)
    
    passed_2d_filter = (passed_ring_filter&passed_pains)
    obj['passed_2d_filter'] = passed_2d_filter
    passed_2d_filter_cnt += passed_2d_filter
    total_cnt += 1

passed_2d_filter_score = passed_2d_filter_cnt / total_cnt * 100
print(f'Passed 2D Filter: {passed_2d_filter_score:.3f}%')
summary['2d'] = passed_2d_filter_score


# ---------------------------- Saving -------------------------------- #

out_path = gen_smi_file[:-3] + 'csv'
table = pd.DataFrame(data)
table.to_csv(out_path, index=False)

summary_path = gen_smi_file[:-4] + '_summary.csv'
summary_table = pd.DataFrame([summary])
summary_table.to_csv(summary_path, index=False)


# ----------------------- RMSD --------------------- #

sdf_path = gen_smi_file[:-3] + 'sdf'
pred_mol_3d = Chem.SDMolSupplier(sdf_path)

if data_set == 'ZINC':
    # Use SMILES of test set generated for molecules processed by OpenBabel
    # (for consistency with other evaluation metrics)
    # Because SMILES produced by our model are also based on OpenBabel
    true_smi_path = f'{data_set_path}/zinc_final_test_smiles.smi'
    true_mol_path = f'{data_set_path}/zinc_final_test_molecules.sdf'
    true_smi = pd.read_csv(true_smi_path, sep=' ', names=['mol', 'frag']).mol.values
    true_mol_3d = Chem.SDMolSupplier(true_mol_path)
    true_smi2mol3d = dict(zip(true_smi, true_mol_3d))
elif data_set == 'HBDrug3D':
    true_smi_path = f'{data_set_path}/HBDrug3D_test_smiles.smi'
    true_mol_path = f'{data_set_path}/HBDrug3D_test_molecules.sdf'
    true_smi = pd.read_csv(true_smi_path, sep=' ', names=['mol', 'frag']).mol.values
    true_mol_3d = Chem.SDMolSupplier(true_mol_path)
    true_smi2mol3d = dict(zip(true_smi, true_mol_3d))
else:
    raise NotImplementedError


def find_exit(mol, num_frag):
    neighbors = []
    for atom_idx in range(num_frag, mol.GetNumAtoms()):
        N = mol.GetAtoms()[atom_idx].GetNeighbors()
        for n in N:
            if n.GetIdx() < num_frag:
                neighbors.append(n.GetIdx())
    return neighbors


def kabsch_numpy(P: np.ndarray, Q: np.ndarray):
    P = P.astype(np.float64)
    Q = Q.astype(np.float64)

    PC = np.mean(P, axis=0)
    QC = np.mean(Q, axis=0)

    UP = P - PC
    UQ = Q - QC

    C = UP.T @ UQ
    V, S, W = np.linalg.svd(C)

    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        V[:, -1] = -V[:, -1]

    R: np.ndarray = V @ W

    t = QC - PC @ R  # (3,)

    return (UP @ R + QC).astype(np.float32), R.astype(np.float32), t.astype(np.float32)


def compute_crmsd(X, Y, aligned=False):
    if not aligned:
        X_aligned, _, _ = kabsch_numpy(X, Y)
    else:
        X_aligned = X
    dist = np.sum((X_aligned - Y) ** 2, axis=-1)
    crmsd = np.sqrt(np.mean(dist))
    return float(crmsd)


# ----------------------------- SC-RDKit -------------------------- #

def calc_sc_rdkit_full_mol(gen_mol, ref_mol):
    try:
        _ = rdMolAlign.GetO3A(gen_mol, ref_mol).Align()
        sc_score = calc_SC_RDKit.calc_SC_RDKit_score(gen_mol, ref_mol)
        return sc_score
    except:
        return -0.5


sc_rdkit_list = []
for i, (obj, pred) in tqdm(enumerate(zip(data, pred_mol_3d)), total=len(data)):
    obj['sc_rdkit'] = None
    if not obj['valid']:
        continue

    true = true_smi2mol3d[obj['true_molecule']]
    score = calc_sc_rdkit_full_mol(pred, true)
    sc_rdkit_list.append(score)
    obj['sc_rdkit'] = score

sc_rdkit_list = np.array(sc_rdkit_list)

sc_rdkit_1 = (sc_rdkit_list > 0.1).sum() / len(sc_rdkit_list) * 100
sc_rdkit_2 = (sc_rdkit_list > 0.2).sum() / len(sc_rdkit_list) * 100
sc_rdkit_3 = (sc_rdkit_list > 0.3).sum() / len(sc_rdkit_list) * 100
sc_rdkit_4 = (sc_rdkit_list > 0.4).sum() / len(sc_rdkit_list) * 100
sc_rdkit_5 = (sc_rdkit_list > 0.5).sum() / len(sc_rdkit_list) * 100
sc_rdkit_6 = (sc_rdkit_list > 0.6).sum() / len(sc_rdkit_list) * 100
sc_rdkit_7 = (sc_rdkit_list > 0.7).sum() / len(sc_rdkit_list) * 100
sc_rdkit_8 = (sc_rdkit_list > 0.8).sum() / len(sc_rdkit_list) * 100
sc_rdkit_9 = (sc_rdkit_list > 0.9).sum() / len(sc_rdkit_list) * 100
sc_rdkit_mean = np.mean(sc_rdkit_list)

print(f'SC_RDKit > 0.1: {sc_rdkit_1:3f}%')
print(f'SC_RDKit > 0.2: {sc_rdkit_2:3f}%')
print(f'SC_RDKit > 0.3: {sc_rdkit_3:3f}%')
print(f'SC_RDKit > 0.4: {sc_rdkit_4:3f}%')
print(f'SC_RDKit > 0.5: {sc_rdkit_5:3f}%')
print(f'SC_RDKit > 0.6: {sc_rdkit_6:3f}%')
print(f'SC_RDKit > 0.7: {sc_rdkit_7:3f}%')
print(f'SC_RDKit > 0.8: {sc_rdkit_8:3f}%')
print(f'SC_RDKit > 0.9: {sc_rdkit_9:3f}%')
print(f'Mean SC_RDKit: {sc_rdkit_mean}')

summary['sc_rdkit_1'] = sc_rdkit_1
summary['sc_rdkit_2'] = sc_rdkit_2
summary['sc_rdkit_3'] = sc_rdkit_3
summary['sc_rdkit_4'] = sc_rdkit_4
summary['sc_rdkit_5'] = sc_rdkit_5
summary['sc_rdkit_6'] = sc_rdkit_6
summary['sc_rdkit_7'] = sc_rdkit_7
summary['sc_rdkit_8'] = sc_rdkit_8
summary['sc_rdkit_9'] = sc_rdkit_9
summary['sc_rdkit_mean'] = sc_rdkit_mean


# ---------------------------- Saving -------------------------------- #

out_path = gen_smi_file[:-3] + 'csv'
table = pd.DataFrame(data)
table.to_csv(out_path, index=False)

summary_path = gen_smi_file[:-4] + '_summary.csv'
summary_table = pd.DataFrame([summary])
summary_table.to_csv(summary_path, index=False)