# VANE
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)

![VANE-2](https://user-images.githubusercontent.com/75185329/227608114-a2cb2ba8-6245-4a26-a64e-7866ed66fa2e.png)

## Welcome to VANE!

VANE, is a python package that process cellular (and tissue) level epigenomic data, and creates celullar-specific protein-protein interaction networks to model omnigenics pathways of disease. 

VANE is able to process either a list of Tagging SNPs, or a list of fine-mapped variants. 

# Installation

### 1. Install dependencies either using PIP or Conda.

#### Installation using PIP:

```
$ git clone https://github.com/BaranziniLab/VANE.git
$ cd VANE/envs/
$ pip install -r requirements_pip.txt
```

#### Installation using Conda: 

```
$ git clone https://github.com/BaranziniLab/VANE.git
$ cd VANE/envs/
$ conda env create -f environment-vane.yml
$ conda activate vane
```

### 2. Download base files

VANE needs multiple cleaned and processed files to work. You can download the base file folder from: (link)

#### Linux:

```
$ cd VANE
$ wget (link)
$ unzip file
```
#### Mac:

```
$ cd VANE
$ curl (link)
$ unzip file
```


# Quick Start

### Complete example of creating cellular-level network for a subset of autoimmune-related tagging variants


