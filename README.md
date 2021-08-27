[![Build Status](https://travis-ci.com/arzwa/wgd.svg?branch=dev)](https://travis-ci.com/arzwa/wgd)

VIB/UGent center for plant systems biology 
[Bioinformatics & evolutionary genomics group](https://www.vandepeerlab.org/)

# `wgd` - simple command line tools for the analysis of ancient whole-genome duplications

`wgd` is a python package designed for the inference of ancient whole-genome duplication (WGD) events from genome and transcriptome assembly. This tutorial elucidates the principle and usage of `wgd`, in accordance with a chapter of a soon-to-be published book. Please refer to the book chapter for more detailed discussion and elucidation.

## Introduction

Whole genome duplication (WGD), a long-known important evolutionary force for eukaryotes, leaves its unique footprints on genome in various aspects, including the explosion of massive gene duplicates, the preservation of syntenic blocks, the burst of reconciled duplication events on certain branches in the concerned phylogenetic tree, the increase of chromosome numbers and et al. `wgd` is the program enabling users to recover such imprinting by means of deciphering the Ks distribution and synteny. We provide an exemplified workflow practice of how to detect such “WGD signature” from a given genomic dataset.We divide the practice into four steps: Step1 Installation, Step2 Data Collection and Preparation, Step3 Construction of Ks Distribution and Step4 Syntenic Analysis.   

## Step1 Installation

To install `wgd` in a virtual environment, the following command lines could be used.

```
git clone <wgd repo>
cd wgd
virtualenv -p=python3 ENV
source ENV/bin/activate
pip install requirements.txt
pip install .
```

When met with permission problem in installation, please try the following command line.

```
pip install -e .
```

If multiply versions of `wgd` were installed in the system, please add the right path of interested version into the environment variables, for example

```
export PATH="$PATH:~/.local/bin/wgd"
```

## Step2 Data Collection and Preparation 

Since the synonymous distance Ks (the number of synonymous substitutions per synonymous site) is a feasible proxy for the age of gene duplicates (see detailed discussions in the book chapter),the protein-coding genes are exactly what we need for the construction of age distribution of whole paranome in Step3. The gene positional information is needed for the profiling of the synteny relationship in Step4. Thus, we need the CDS (protein-coding sequences) and GFF (General Feature Format) files of interested species, here as *Vitis vinifera and *Amborella trichopoda, both of which are downloaded from [PLAZA](https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v4_5_dicots/download/). A cleaning process for the CDS file is recommended, including: 1) only the longest transcripts are retained if alternatives are available; 2) genes with identical sequences or IDs are removed to eliminate redundancy; 3) the sequence length of CDS should be dividable by three and only contain ACGT characters while not contain any stop codons in the sequence (only at the end of sequence is allowed). The cleaned files are in the [data](https://github.com/heche-psb/wgd/tree/dev/data) directory.

## Citation
 
Please cite us at https://doi.org/10.1093/bioinformatics/bty915

```
Zwaenepoel, A., and Van de Peer, Y. 
wgd - simple command line tools for theanalysis of ancient whole genome
duplications. Bioinformatics., bty915,
https://doi.org/10.1093/bioinformatics/bty915
```

For citation of the tools used in wgd, please consult the documentation at
https://wgd.readthedocs.io/en/latest/index.html#citation.

