
# VCboost: reducing false positives in long-read variant calling for SNP and indel detection in challenging genomic regions
 
Email: holyterror@163.com  

----

## Introduction

VCboost effectively filters out a substantial number of false positive sites, leading to a significant improvement in accuracy and F1 score with minimal loss of true positive sites.

----

## Contents

* [Introduction](#introduction)
* [Usage](#usage)
* [Pre-trained Models](#pre-trained-models)
  * [Guppy3,4 Model](#pre-trained-models)
* [Installation](#installation)
  + [Option 1. Build an anaconda virtual environment](#option-1-build-an-anaconda-virtual-environment)
* [Quick Demo](#quick-demo)

* [Training Data](docs/training_data.md)
* [Model Training](docs/pileup_training.md)

----



---

## Pre-trained Models

In a bioconda installation, models are in `{CONDA_PREFIX}/bin/models/`.

|       Model name       |  Platform   |                       Training samples                       |   Basecaller   |
|:----------------------:| :---------: | :----------------------------------------------------------: |:--------------:|
|     r941_g422_1245     |     ONT r9.4.1    |                         HG001,2,4,5                          | Guppy4.2.2 |
|     r941_g422_1235     |     ONT r9.4.1    |                         HG001,2,3,5                          | Guppy4.2.2 |





## Installation

### Build an anaconda virtual environment


```bash
# Clone the Repository
git clone https://github.com/oranges7/VCboost.git
# Navigate to the Project Directory
cd ${repository}
# You can use Conda, pip, or other tools for installation.
conda create --name vcboost --file requirements.txt
# Activate the virtual environment
conda activate vcboost
# Display help information for the vcboost
sh vcboost.sh -h
```

## Usage

### General Usage


```bash
sh vcboost.sh \
  -o ${OUTPUT_PATH} \
  -b ${BAM_FILE} \ 
  -v ${ORIGINAL_VCF_FILE} \
  -m ${MODEL_PREFIX} \ 
  -r ${REFERENCE}

## VCboost final output file: ${OUTPUT_PATH}/vc_boost.vcf
```

### Options

**Required parameters:**

```bash
  Options:
  -o, -out_path        Output path.
  -b, -bam_file        BAM file path.
  -v, -vcf             VCF file path.
  -m, -model prefix    Model path.
  -r, -ref_path        Reference file path.

```

**Other parameters:**

```bash
  -t, -threads     Number of threads.The default is 32.
  -c, -contig      Contig to process.The default is chr1-22.
  -p, -phase       Disable phase.
  -h, -help        Display this help message.
```

**Train parameters:**

```bash
  -train   Train mode.
  -w, -work_path   Working directory path.
  -a, -aim_vcf     Aim VCF file path.
  -b, -bam_file    BAM file path.
  -r, -ref         Reference file path.
  -q, -vcf         Benchmark VCF file path.
  -m, -mode        Mode of operation.You can choose snp or indel. The default is both.
  -j, -object      Object to process.
  -p, -phase       Enable phase.
```
----

## Training Data

The training and testing datasets were derived from the public Genome in a Bottle (GIAB) and Human Pangenome Reference Consortium (HPRC) from individuals HG001, HG002, HG003, HG004, and HG005. The original nanopore sequencing data generated using R9.4.1 flow cells and basecalled with Guppy v4.2.2 were aligned to the reference GRCh38 to map sequencing reads.  Currently, two models have been trained and tested: one based on training with HG001, HG002, HG003, HG005, and their downsampled counterparts; and the other based on training with HG001, HG002, HG004, HG005, and their downsampled counterparts.

----
