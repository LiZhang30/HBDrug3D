# HBDrug3D
Descriptions:

## Molecular type, number, and atom type:

PROTAC: 5,607    

ADC: 254

PDC: 

This dataset includes the following atom types:
C, O, N, F, S, Cl, Br, I, and P.

## Data splits:
Training: 5,447 Protac examples + 214 ADC examples
Validation: 100 examples (80 Protac + 20 ADC)
Test: 100 examples (80 Protac + 20 ADC)

File Structure: The tarDrug dataset is stored in HDF5 format, organized with the following hierarchy:

Top-level groups for each data split: 'train', 'val', and 'test'.
Within each split, subgroups for each component type:
'train/protac', 'train/adc', and similarly for the validation and test splits, allowing for organized access to Protac and ADC data separately within each split.

---

## Data structure (each molecule entry)

Each entry in the HDF5 file represents a unique molecule and contains the following datasets:

id: A unique identifier for the molecule entry.

smiles: A list of SMILES strings for different parts of the molecule, formatted as [whole_smi, frag1_smi, linker_smi, frag2_smi].

pos: Atom positions for the molecule, stored as a [N, 3] array, where N is the number of atoms. Each entry provides the 3D coordinates of an atom.

atoms_atomic_numbers: Atom types stored as atomic numbers with shape [N], where each value represents the atomic number of the corresponding atom.

frag1_idx: An array containing the indices of atoms that belong to frag1 within the molecule.

frag2_idx: An array containing the indices of atoms that belong to frag2.

linker_idx: An array of indices for atoms that belong to the linker part of the molecule.

linker_size: An integer value indicating the total number of atoms in the linker segment.

bond: Bond information for the molecule, stored as an array of tuples with each tuple in the format (start_idx, end_idx, bond_type), where bond_type could be 'SINGLE', 'DOUBLE', 'TRIPLE', or 'AROMATIC'.

anchors: A boolean array of length N, where each entry indicates whether the atom acts as an anchor point connecting the linker to one or both fragments (True if the atom is an anchor, False otherwise).

---

## Implementation baselines:
