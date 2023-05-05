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
$ wget https://ucsf.box.com/shared/static/b5hbzpl84v73lz0d2tt72gsi2vvftpz3.zip
$ unzip b5hbzpl84v73lz0d2tt72gsi2vvftpz3.zip
```
#### Mac:

```
$ cd VANE
$ curl https://ucsf.box.com/shared/static/b5hbzpl84v73lz0d2tt72gsi2vvftpz3.zip
$ unzip b5hbzpl84v73lz0d2tt72gsi2vvftpz3.zip
```


> **A:** Your VANE folder and files, should look exactly like *this*.

    VANE
    ├── ...
    ├── base_files                   
    │   ├── cell_networks          # Folders with all the cell networks
    │   ├── clean_snp2gene.db         # SNP2 Gene SQL db from OpenTargets
    │   └── gene_protein_ensembl_map_table.tsv                # Gene-protein-ENSG_id mapping tables
    ├── VANE.py
    ├── snp_parser.py
    ├── regulome_db_parsing.py
    └── ld_parser.py

# Quick Start

### Complete example of creating cellular-level network for a subset of autoimmune-related tagging variants

#### Epigenomic Analysis
```python
from VANE import VANE

ms_variants = ['rs35486093','rs11256593', 'rs112344141', 'rs1323292', 'rs72928038', 'rs78727559', 'rs1800693', 'rs10801908', 'rs6670198',
'rs62420820', 'rs1738074', 'rs4939490', 'rs9843355', 'rs11809700', 'rs35540610', 'rs1026916', 'rs1014486', 'rs6589706', 'rs11079784',
'rs11749040', 'rs631204', 'rs4808760', 'rs12478539', 'rs7977720', 'rs3809627', 'rs12365699', 'rs6032662', 'rs60600003', 'rs2546890', 'rs1465697']

MS_variant = VANE(list_of_variants = ms_variants)
df_with_positions = MS_variant.process_tag_variants(ld=0.7, population = "1000GENOMES:phase_3:CEU", verbose=False)

snp_2_gene = MS_variant.query_snp_opentargets(df_with_positions, "base_files/clean_snp2gene.db")
gene_to_rsid = MS_variant.merge_and_clean(snp_2_gene, df_with_positions)
MS_variant.regulomedb_scoring(save=True, verbose=True)


cell_list = ['B', 'T', 'mono', 'astrocyte', 'neutrophil', 'keratinocyte', 'melanocyte']
organs_list = ['brain','spleen','thymus', 'lymphoid tissue', 'lymph node', 'skin']

MS_variant.process_cells(cell_list)
MS_variant.plot_cells()

MS_variant.process_organs(organs_list)
MS_variant.plot_organs()

```
#### Network Analysis

Calculate P value and CI of the cellular-network of the epigenomic-relevant target genes:

```python
from VANE import VANE_network

T_cell_network = VANE_network(network_folder = 'base_files/cell_networks/', df_for_network=MS_variant.df_for_network)
T_cell_network.create_cellular_network(cell = 'T', epigenome_threshold = 0)
CI, pval = T_cell_network.calculate_p_value()
T_cell_network.plot_bootstrap()

```




