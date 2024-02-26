
# VCboost - Improving the Precision of Variant Calling for Long-Read Sequencing by Filtering
 
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

### Option 1.  Build an anaconda virtual environment




### Option 2.  Bioconda


```bash
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "clair3"
# replace clair3 by clair3-illumina for using illumina data
conda create -n clair3 -c bioconda clair3 python=3.9.0 -y
conda activate clair3

# run clair3 like this afterward
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

run_clair3.sh \
  --bam_fn=input.bam \                 ## change your bam file name here
  --ref_fn=ref.fa \                    ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \ 
  --output=${OUTPUT_DIR}               ## output path prefix 
```

Check [Usage](#Usage) for more options. 



## Usage

### General Usage


```bash
sh gen_onehot_data.sh \
  -o ${OUTPUT_PATH} \
  -b ${BAM_FILE} \
  -v ${ORIGINAL_VCF_FILE} \
  -m ${MODEL_PATH}/${MODEL_PREFIX} \
  -r ${REFERENCE}

## VCboost final output file: ${OUTPUT_PATH}/vcboost.vcf
```

### Options

**Required parameters:**

```bash
  Options:
  -o, --out_path    Output path.
  -b, --bam_file    BAM file path.
  -v, --vcf         VCF file path.
  -m, --model       Model path.
  -r, --ref_path    Reference file path.

```

**Other parameters:**

```bash
  -t, --threads     Number of threads.The default is 32.
  -c, --contig      Contig to process.The default is chr1-22.
  -p, --phase       Disable phase.
  -h, --help        Display this help message.
```

**Train parameters:**

```bash
  -w, --work_path   Working directory path.
  -a, --aim_vcf     Aim VCF file path.
  -b, --bam_file    BAM file path.
  -r, --ref         Reference file path.
  -q, --vcf         Benchmark VCF file path.
  -m, --mode        Mode of operation.You can choose snp or indel. The default is both.
  -j, --object      Object to process.
  -p, --phase       Enable phase.
```
----

## Training Data

The training and testing datasets were derived from the public Genome in a Bottle (GIAB) and Human Pangenome Reference Consortium (HPRC) from individuals HG001, HG002, HG003, HG004, and HG005. The original nanopore sequencing data generated using R9.4.1 flow cells and basecalled with Guppy v4.2.2 were aligned to the reference GRCh38 to map sequencing reads.  Currently, two models have been trained and tested: one based on training with HG001, HG002, HG003, HG005, and their downsampled counterparts; and the other based on training with HG001, HG002, HG004, HG005, and their downsampled counterparts.

----
