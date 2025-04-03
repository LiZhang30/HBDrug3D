import argparse
import pandas as pd
from tqdm.auto import tqdm
from difflib import SequenceMatcher

from rdkit import Chem
from rdkit.Chem import rdmolops, DataStructs, RDKFingerprint

########### define your constants ###########

def split_protac(raw_file_dir, save_path):
    unimatch_num = 0  # record unique match
    null_match_idxs, mul_matches_idxs = [], []  # record not unique match
    data = pd.read_excel(raw_file_dir, header=0)  # read the raw file
    result_to_save = []  # buffer to save data
    final_records_num = 0  # the number of saved samples

    for i in tqdm(range(len(data))):
        protac_id = data.iloc[i]['PROTAC ID']
        protac_smi = data.iloc[i]['PROTAC Smiles']
        warhead_id = data.iloc[i]['Warhead ID']
        warhead_smi = data.iloc[i]['Warhead Smiles']
        linker_id = data.iloc[i]['Linker ID']
        linker_smi = data.iloc[i]['Linker Smiles']
        e3lig_id = data.iloc[i]['E3 Ligand ID']
        e3lig_smi = data.iloc[i]['E3 Ligand Smiles']
        protac_mol = Chem.MolFromSmiles(protac_smi)
        warhead_mol = Chem.MolFromSmiles(warhead_smi)
        linker_mol = Chem.MolFromSmiles(linker_smi)
        e3lig_mol = Chem.MolFromSmiles(e3lig_smi)

        # match linker from protac
        match = protac_mol.GetSubstructMatches(linker_mol)
        if not match:
            null_match_idxs.append(i)  # record null match
            continue
        elif len(match) == 1:
            unimatch_num += 1  # record unique match
        else:
            mul_matches_idxs.append(i)  # record multiple matches

        frag_warhead_fpsmi_buff, frag_e3lig_fpsmi_buff = [], []  # record fingerprint similarity
        warhead_save_buff, linker_save_buff, e3lig_save_buff = [], [], []  # record smiles
        # parse the protac with matched linker/linkers
        for one_match in match:
            one_match_atom_idxs = set(one_match)
            # extract bonds between linker and other two frags
            bonds_to_remove = [
                bond.GetIdx()
                for bond in protac_mol.GetBonds()
                if (bond.GetBeginAtomIdx() in one_match_atom_idxs) != (bond.GetEndAtomIdx() in one_match_atom_idxs)
            ]
            #assert len(bonds_to_remove) == 2, f"Expected 2 bonds to remove, but got {len(bonds_to_remove)}"
            if len(bonds_to_remove) != 2:
                continue

            # split protac into three frags by linker (warhead, linker and e3 ligand)
            RWprotac_mol = Chem.RWMol(protac_mol)
            for bond_idx in sorted(bonds_to_remove, reverse=True):
                bond = RWprotac_mol.GetBondWithIdx(bond_idx)
                RWprotac_mol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            frags = rdmolops.GetMolFrags(RWprotac_mol, asMols=True, sanitizeFrags=False)
            #assert len(frags) == 3, f"Expected 3 frags, but got {len(frags)}"
            if len(frags) != 3:
                continue
            frags_smi = [Chem.MolToSmiles(frag, canonical=True) for frag in frags]

            # tag unknown frags as warhead or e3 ligand
            warhead_fp = RDKFingerprint(warhead_mol)
            linker_fp = RDKFingerprint(linker_mol)
            e3lig_fp = RDKFingerprint(e3lig_mol)
            get_warhead, get_linker, get_e3lig = None, None, None  # empty the buffer
            get_frag_warhead_fpsmi, get_frag_e3lig_fpsmi = None, None  # empty the buffer
            for idx, frag in enumerate(frags):
                frag_fp = RDKFingerprint(frag)
                frag_linker_fpsmi = DataStructs.FingerprintSimilarity(frag_fp, linker_fp)
                frag_warhead_fpsmi = DataStructs.FingerprintSimilarity(frag_fp, warhead_fp)
                frag_e3lig_fpsmi = DataStructs.FingerprintSimilarity(frag_fp, e3lig_fp)
                if frag_linker_fpsmi > 0.99:  # 0.99 is the threshold
                    get_linker = frags_smi[idx]
                elif (frag_warhead_fpsmi>=0.6) and (frag_warhead_fpsmi>frag_e3lig_fpsmi):  # 0.6 is the threshold
                    get_warhead = frags_smi[idx]
                    get_frag_warhead_fpsmi = frag_warhead_fpsmi
                elif (frag_e3lig_fpsmi>=0.6) and (frag_e3lig_fpsmi>frag_warhead_fpsmi): # 0.6 is the threshold
                    get_e3lig = frags_smi[idx]
                    get_frag_e3lig_fpsmi = frag_e3lig_fpsmi
                '''else:
                    print(f'frag_warhead_fpsmi:', frag_warhead_fpsmi, '   frag_e3lig_fpsmi:', frag_e3lig_fpsmi)'''
            if get_warhead and get_linker and get_e3lig:
                warhead_save_buff.append(get_warhead)  # save smiles
                linker_save_buff.append(get_linker)  # save smiles
                e3lig_save_buff.append(get_e3lig)  # save smiles
                frag_warhead_fpsmi_buff.append(get_frag_warhead_fpsmi)  # save similarity
                frag_e3lig_fpsmi_buff.append(get_frag_e3lig_fpsmi)  # save similarity

        # skip cases where no matches meet the conditions
        if (len(warhead_save_buff))*(len(e3lig_save_buff)) != 0:
            #assert frag_warhead_fpsmi.index(max(frag_warhead_fpsmi))==frag_e3lig_fpsmi.index(max(frag_e3lig_fpsmi)), f"It should be equal, but it isn't."
            if frag_warhead_fpsmi_buff.index(max(frag_warhead_fpsmi_buff)) == frag_e3lig_fpsmi_buff.index(max(frag_e3lig_fpsmi_buff)):
                final_records_num += 1
                warhead_fsave, linker_fsave, e3lig_fsave = None, None, None  # empty the buffer
                '''print(f'frag_warhead_fpsmi_buff:', len(frag_warhead_fpsmi_buff), f'frag_e3lig_fpsmi_buff:', len(frag_e3lig_fpsmi_buff))
                print(frag_warhead_fpsmi_buff.index(max(frag_warhead_fpsmi_buff)), frag_e3lig_fpsmi_buff.index(max(frag_e3lig_fpsmi_buff)))'''
                idx_from_max_value = frag_warhead_fpsmi_buff.index(max(frag_warhead_fpsmi_buff))
                warhead_fsave = warhead_save_buff[idx_from_max_value]
                linker_fsave = linker_save_buff[idx_from_max_value]
                e3lig_fsave = e3lig_save_buff[idx_from_max_value]

                result_to_save.append({
                    'PROTAC ID': protac_id,
                    'PROTAC Smiles': protac_smi,
                    'Warhead ID': warhead_id,
                    'Warhead Smiles': warhead_fsave,
                    'Linker ID': linker_id,
                    'Linker Smiles': linker_fsave,
                    'E3 Ligand ID': e3lig_id,
                    'E3 Ligand Smiles': e3lig_fsave,
                })

    print(f'null_match:', len(null_match_idxs), f'uni_match:', unimatch_num, f'mul_matches:', len(mul_matches_idxs), f'final_records:', final_records_num)
    result_to_save_df = pd.DataFrame(result_to_save)
    result_to_save_df.to_excel(save_path, index=False)


def split_adc(raw_file_dir, save_path):
    data = pd.read_excel(raw_file_dir, header=0)
    result_to_save = []  # buffer to save data
    n_fail_f, n_fail_l = 0, 0  # record error

    for i in tqdm(range(len(data))):
         adc_id = data.iloc[i]['ADC ID']
         adcdb_id = data.iloc[i]['ADCdb ID']
         clp_name = data.iloc[i]['Conjugate-Linker-Payload Name']
         clp_smi = data.iloc[i]['Conjugate-Linker-Payload Smiles']
         conjugate_type = data.iloc[i]['Conjugate Type']
         conjugate_smi = data.iloc[i]['Conjugate Smiles']
         linker_name = data.iloc[i]['Linker Name']
         #linker_smi = data.iloc[i]['Linker Smiles']
         payload_name = data.iloc[i]['Payload Name']
         payload_smi = data.iloc[i]['Payload Smiles']

         clp_mol = Chem.MolFromSmiles(clp_smi)
         clp_atoms_idxs = set(range(clp_mol.GetNumAtoms()))
         conjugate_mol = Chem.MolFromSmiles(conjugate_smi)
         payload_mol = Chem.MolFromSmiles(payload_smi)

         # the process of filtering multi matches has been considered.
         conjugate_matches = clp_mol.GetSubstructMatches(conjugate_mol)
         payload_matches = clp_mol.GetSubstructMatches(payload_mol)
         fail_f = False
         if len(conjugate_matches)==0 or len(payload_matches)==0:
            print('data ' + str(adc_id) + ' fragment error.')
            n_fail_f += 1
            fail_f = True

         linker_smi_to_save = ''  # init the buff
         found_valid_linker = False  # the flag to jump for
         for conjugate_idxs in conjugate_matches:
             for payload_idxs in payload_matches:
                 frags_idxs = set(conjugate_idxs + payload_idxs)
                 linker_idxs = clp_atoms_idxs - frags_idxs

                 linker_mol_to_save = Chem.RWMol()
                 atom_map = {}
                 for j, atom_idx in enumerate(linker_idxs):
                     atom = clp_mol.GetAtomWithIdx(atom_idx)
                     atom_map[atom_idx] = linker_mol_to_save.AddAtom(atom)

                 for bond in clp_mol.GetBonds():
                    start_idx, end_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    if {start_idx, end_idx} <= set(linker_idxs):
                         linker_mol_to_save.AddBond(atom_map[start_idx], atom_map[end_idx], bond.GetBondType())

                 linker_smi_to_save = Chem.MolToSmiles(linker_mol_to_save)
                 if linker_smi_to_save is None or '.' in linker_smi_to_save:
                     linker_smi_to_save = ''
                 else:
                     found_valid_linker = True
                     break
             if found_valid_linker:
                 break
         if (linker_smi_to_save=='') & (~fail_f):
            print('data ' + str(adc_id) + ' linker error.')
            n_fail_l += 1

         result_to_save.append({
            'ADC ID': adc_id,
            'ADCdb ID': adcdb_id,
            'Conjugate-Linker-Payload Name': clp_name,
            'Conjugate-Linker-Payload Smiles': clp_smi,
            'Conjugate Type': conjugate_type,
            'Conjugate Smiles': conjugate_smi,
            'Linker Name': linker_name,
            'Linker Smiles': linker_smi_to_save,
            'Payload Name': payload_name,
            'Payload Smiles': payload_smi,
         })
    print('Num of fragment failed: ' + str(n_fail_f))
    print('Num of linker failed: ' + str(n_fail_l))
    result_to_save_df = pd.DataFrame(result_to_save)
    result_to_save_df.to_excel(save_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw_path', type=str, required=True)
    parser.add_argument('--dataset', type=str, default='adc', choices=['protac', 'adc'])
    parser.add_argument('--dest', type=str, required=True)
    args = parser.parse_args()

    # example for run: python filter_rawdata.py --raw_path TarDrug3D/PROTAC-DB-1.xlsx --dest IndexFile/protac.xlsx
    # example for run: python filter_rawdata.py --raw_path TarDrug3D/ADCdb-1.xlsx --dest IndexFile/adc.xlsx
    if args.dataset == 'protac':
        split_protac(args.raw_path, args.dest)
    elif args.dataset == 'adc':
        split_adc(args.raw_path, args.dest)
    else:
        raise NotImplementedError