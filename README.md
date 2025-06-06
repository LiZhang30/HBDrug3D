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

### Create environment:

You can quickly create a new environment using the method below, or start from scratch, to perform model evaluation.
```
conda env create -f env.yml
conda activate HBDrug3D
```

## Prepare Dataset

***Step1***  ->  Directly download the raw data from the data platform or use scripts in the collection directory to obtain the raw data.
<br/>
***Step2***  ->  Use filter_rawdata.py in the filter directory to preprocess the raw data.
<br/>
***Step3***  ->  First, use [Schrodinger ConfGen](https://www.schrodinger.com/platform/products/confgen/) to generate conformations, then run filter_conformations.py in the filter directory to filter out low-quality conformations.
<br/>
***Step4***  ->  Finally, use process_rawdata.py to construct samples and save them in HDF5 format.
<br/>

## Dataset

HBDrug3D consists of three subsets: PROTAC, ADC, and PDC. The HDF5 file for each subset, along with the corresponding conformations and splits, can be downloaded from [HBDrug3D](https://drive.google.com/drive/folders/1XdgJPCcVfQfMFQN8YnqXITgvpxl15A2v).

| Subset | PROTAC | ADC | PDC |
|------------|-----------------------------------|-----------------------------------|-----------------------------------|
| File | protac.h5 | adc.h5 | pdc.h5 |
| Number of molecules | 5,607 | 246 | 426 |
| Number of entries | 28,970 | 3,181 | 26,673 |
| Splits (train/val/test) | 28,170/400/400 | 2,981/100/100 | 25,873/400/400 |

***Note***  ->  This dataset includes these atom types: C, O, N, F, S, Cl, Br, I, and P.

### Data structure:

Each entry in the HDF5 file represents a unique sample, and the details of each field are described as follows:

->  ***id***: A unique identifier for the entry.
<br/>
->  ***smiles***: A list of SMILES strings for different parts of the molecule, formatted as [whole_smi, frag1_smi, linker_smi, frag2_smi].
<br/>
->  ***pos***: Atom positions for the molecule, stored as a [N, 3] array, where N is the number of atoms. Each entry provides the 3D coordinates of an atom.
<br/>
->  ***atoms_atomic_numbers***: Atom types stored as atomic numbers with shape [N], where each value represents the atomic number of the corresponding atom.
<br/>
->  ***frag1_idx***: An array containing the indices of atoms that belongs to frag1 within the molecule.
<br/>
->  ***frag2_idx***: An array containing the indices of atoms that belongs to frag2 within the molecule.
<br/>
->  ***linker_idx***: An array of indices for atoms that belongs to the linker part of the molecule.
<br/>
->  ***linker_size***: An integer value indicating the total number of atoms in the linker segment.
<br/>
->  ***bond***: Bond information for the molecule, stored as an array of tuples with each tuple in the format (start_idx, end_idx, bond_type), where bond_type could be 'SINGLE', 'DOUBLE', 'TRIPLE', or 'AROMATIC'.
<br/>
->  ***anchors***: A boolean array of length N, where each entry indicates whether the atom acts as an anchor point connecting the linker to one or both fragments (True if the atom is an anchor, False otherwise).

## Benchmark

### Evaluation metrics:

->  ***Validity*** (chemical soundness)
<br/>
->  ***Novelty*** (non-overlap with training samples)
<br/>
->  ***Recovery Rate*** (exact matches with reference structures)
<br/>
->  ***Standardized 2D Filtering***
<br/>
->  ***SC_RDKit*** (rigorous 3D similarity evaluation through both geometric alignment and chemical feature matching with reference compounds)

### Evaluation steps:

***Step1 (Preprocessing)***  ->  First, use preprocess.py in the evaluation directory to generate the files required for step4.
<br/>
***Step2 (Sampling)***  ->  Next, use the trained model to sample linkers for the given fragments, following the sampling logic in sample.py in the evaluation directory.
<br/>
***Step3 (Reformatting)***  ->  Then, use reformat_obabel.py in the evaluation directory to reformat the sampling results.
<br/>
***Step4 (Calculating metrics)***  ->  Finally, use compute_metrics.py in the evaluation directory to obtain various metrics.

### Implementation baselines:

***DeLinker***  ->  We utilized the original source code ([DeLinker](https://github.com/oxpig/DeLinker)) under TensorFlow 1.10 framework, maintaining the same framework and hyperparameter configurations as the published work.
<br/>
***3DLinker***  ->  We utilized the original source code ([3DLinker](https://github.com/YinanHuang/3DLinker)) under PyTorch 1.11.0 framework, maintaining the same framework and hyperparameter configurations as the published work.
<br/>
***DiffLinker***  ->  We implemented DiffLinker using the official source code ([DiffLinker](https://github.com/igashov/DiffLinker)) with PyTorch 2.0, strictly maintaining the original model architecture and hyperparameters. The only modification involved adjusting the number of training epochs (from 500 to 600) to ensure proper convergence while keeping all other training parameters unchanged (batch size=32, learning rate=1e-4, 1000 diffusion steps). The model was trained on NVIDIA A100 GPUs with mixed-precision acceleration, following the same data preprocessing pipeline as described in the original work.
<br/>
***LinkerNet***  ->  We faithfully reproduced LinkerNet using the original source code ([LinkerNet](https://github.com/guanjq/LinkerNet)) implemented in PyTorch 2.5.1, maintaining all published model architectures without modification.
