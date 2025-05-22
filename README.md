<p align="center">
    <img src="https://raw.githubusercontent.com/LiZhang30/HBDrug3D/main/images/logo.svg" alt="log" width="60%" />
    <br/>
</p>

<a href="https://pytorch.org/get-started/locally/"><img alt="PyTorch" src="https://img.shields.io/badge/PyTorch-ee4c2c?logo=pytorch&logoColor=white"></a>
</p>

_HBDrug3D_, the first benchmark dataset for accerating heterobifunctional drug linker design.

## Introduction

This is the official code repository for the paper 'HBDrug3D: A Dataset and Benchmark for AI-Driven Heterobifunctional Molecule Design', which aims to accelerate heterobifunctional drug discovery. The dataset covers three types of heterobifunctional molecules, as listed below:

| Molecule | Paper link | Data platform |
|------------|---------------------------------------------|-------------------------------------------|
| PROTAC | https://doi.org/10.1093/nar/gkae768 | http://cadd.zju.edu.cn/protacdb/ |
| ADC | https://doi.org/10.1093/nar/gkad831 | https://adcdb.idrblab.net/ |
| PDC | https://doi.org/10.1093/nar/gkae859 | https://pdcdb.idrblab.net/ |

It contains `6,279 heterobifunctional molecules` and `58,824 conformations`, spanning three druggable mechanisms，as shown below.
<p align="center">
    <img src="https://raw.githubusercontent.com/LiZhang30/HBDrug3D/main/images/Figure 1.png" alt="Cover" width="80%" />
    <br/>
</p>

Four state-of-the-art linker design models are benchmarked on the HBDrug3D dataset, as listed below.
| Model | Paper link | Github |
|------------|---------------------------------------------|-------------------------------------------|
| DeLinker | https://pubs.acs.org/doi/10.1021/acs.jcim.9b01120 | https://github.com/oxpig/DeLinker |
| 3DLinker | https://arxiv.org/abs/2205.07309 | https://github.com/YinanHuang/3DLinker |
| DiffLinker | https://doi.org/10.1038/s42256-024-00815-9 | https://github.com/igashov/DiffLinker |
| LinkerNet | https://openreview.net/forum?id=6EaLIw3W7c | https://github.com/guanjq/LinkerNet |

## Installation

### Create environment with basic packages

```
conda env create -f environment.yml
conda activate cbgbench
```

### Install pytorch and torch_geometric

```
conda install pytorch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2 pytorch-cuda=11.8 -c pytorch -c nvidia
conda install pyg pytorch-scatter pytorch-cluster -c pyg
```

### Install tools for chemistry

```
# install rdkit, efgs, obabel, etc.
pip install --use-pep517 EFGs
pip install biopython
pip install lxml
conda install rdkit openbabel tensorboard tqdm pyyaml easydict python-lmdb -c conda-forge

# install plip
mkdir tools
cd tools
git clone https://github.com/pharmai/plip.git
cd plip
python setup.py install
alias plip='python plip/plip/plipcmd.py'
cd ..
```

## Prepare Dataset

***Step1***  ->  Directly download the raw data from the data platform or use scripts in the collection directory to obtain the raw data.
<br/>
***Step2***  ->  Use filter_rawdata.py in the filter directory to preprocess the raw data.
<br/>
***Step3***  ->  First, use [Schrodinger ConfGen](https://www.schrodinger.com/platform/products/confgen/) to generate conformations, then run filter_conformations.py in the filter directory to filter out low-quality conformations.
<br/>
***Step4***  ->  Finally, use process_rawdata.py to construct samples and save them in H5 format.
<br/>

## Dataset
[HBDrug3D.h5](https://drive.google.com/drive/folders/1XdgJPCcVfQfMFQN8YnqXITgvpxl15A2v++)


### Dataset:

File Structure: The HBDrug3D dataset is stored in HDF5 format, organized with three subsets:
form protac.h5, adc.h5, pdc.h5

Molecule: _PROTAC_: 5,607; _ADC_: 254; _PDC_: 426
Atoms: This dataset includes these atom types: C, O, N, F, S, Cl, Br, I, and P.
Splits: PROTACs (28,170/400/400), ADCs(2,981/100/100), and PDCs (25,873/400/400).
from .txt files to get indxs in protac.h5 


### Data structure (each molecule entry)

Each entry in the HDF5 file represents a unique molecule and contains the following datasets:
```
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
```

## Benchmark
### Evaluation metrics:

### Implementation baselines:
