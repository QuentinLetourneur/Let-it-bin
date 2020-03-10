# Let-it-bin

## Contents

- [Introduction](#introduction)
- [Prerequisites](#prerequisites)
- [Pipeline inputs and options](#Pipeline_inputs_and_options)

## Introduction

Let-it-bin allows to perform the binning of  metagenomic short paired-end reads from multiple samples into species.
The pipeline take raw paired-end reads from multiple samples as primary input (the pairs names MUST finish with _{1,2}.fq/fastq). It comprise 4 major steps, reads preprocessing, assembly, binning and evaluation.
The pipeline can be started from the second or third step by adding arguments to the command line, provided that you have the needed inputs.
You have to select the binning softwares that will run in the following list :
  binsanity, canopy, concoct, cocacola, maxbin, metabat, metabat2 and metagen
You just have to prefix the name of the wanted programms with '--' (Ex : --concoct).
If you want to use them all just use --all T

When path are needed please give **full path**.

If you run this pipeline on a cluster (what I recommend) you can use the given .config file to specify allocated memory per task and other cluster options. Memory values have been placed based on experience but can be changed.  
Be it locally or on a cluster **be sure to add the full path to let-it-bin.simg** (more details in the next section) in the config file


The output directory have the following layout :
```
[out]
|
|__assembly         Assembly and contigs annotation
|__Binnings         Folder of each chosen binning software
|  |__Metabat
|     |__checkm_res
|__cleaned_reads
|__khmer_res
|__mapping
```
## Prerequisites

To run this pipeline you will need Nextflow and Singularity (tested with version 19.10.0.5170 and 3.5.0 respectively).  
Here are the links to the installation instruction for [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html) and [Singularity](https://github.com/sylabs/singularity/blob/master/INSTALL.md)

The singularity image can be downloaded here (warning: the file is heavy ~ 1.9Go):
```
wget ftp://shiny01.hosting.pasteur.fr/pub/let-it-bin.simg
```
A recipe file named Singularity is also given.  
To build the image on an unix system move to let-it-bin repository and lauch  
`sudo singularity build let-it-bin.simg Singularity`  
This will take at least an hour.

## Pipeline inputs and options
```
    GENERAL ARGUMENTS :

  --reads [PATH] Directory containing unzipped paired reads files.
Only needed if you start from raw reads or reads from which contaminant have been removed
  --nb_samples [INT] Number of samples you are using
  --out [PATH] Directory were will be stored the results of the pipeline
  --sim_data [CHAR] Can be either F (Default) or T. Will change the execution of the pipeline depending on the analysed data (simulated or not).
  --cpus [INT] Number of cpus used for task with multiprocessing (Default 4)
  --min_contigs_length [INT] Minimum contigs length in base to be passed to binning programms (Default 1000)
  --nb_ref [INT] If you use simulated data specify the total number of different genomes present in the samples
  --dastool [CHAR] Can be either T (Default) or F. If you use multiple binning softwares you can use dastool to combine these results and try to extract best bins corresponding to the same microorganism.
  --local_scratch [CHAR] Can be either T (Default) or F. If you are on TARS or on a cluster with a /local/scratch space on nodes. You can use this option to speed up the execution of post processing of binning result for Canopy, Concoct, Cocacola and Metagen.
  --tmp_checkm [PATH] Directory were will be stored CheckM temporary files. The path length souhldn't exeed 65 chars. (Default [out]/tmp_checkm)
  --help Print the help message

    READS PREPROCESSING :

  --contaminant [PATH] Path and prefix of the bowtie2 index files of contaminant sequences (HAVE TO be computed before lauching the pipeline)
  --minlength [INT] Minimum length for trimed contigs (Default 45)
  --alienseq [PATH] Fasta file containing adaptaters sequence for AlienTrimmer
  --cleaned_readsDir [PATH] Folder were will be stored reads that have been filtered to eliminate contaminant and trimmed (Default [out]/cleaned_reads)
If there are already fastq files in the folder it will take it as input for khmer and skip the cleaning step
  --filt_readsDir [PATH] Directory containing Khmer results.
IF SPECIFIED the workflow will start at the assembly by taking as input the filtered fastq files in the directory.

    ASSEMBLY :

  --cpuassembly [INT] Number of cpus used for reads assembly (Default 10)
  --memassembly [INT] Quantity of RAM in Mb used for reads assembly. Default 160000 Mb with 2 retries. If the wanted memory is <= 160000, the allocated memory will grow according to the following formula : number of retry * [memmapping]
Else no retry will be done
  --qos_assembly [STRING] If you run the pipeline on a cluster with SLURM you can specify the queue of submision (qos) to use : fast, normal (Default) or long (on TARS you can only use 5 cpus in long)
  --mode [STRING] Name of the assembler to be used. Can be either spades (Default), clc, megahit or ray
  --multi_assembly [CHAR] By default a co-assembly of the samples is done. If you want to change that behaviour and do an assembly by sample set this parameter to T (Default F). The generated contigs will then be pulled in one file and filtered to lessen the redundancy but eliminating it is hard so there migth still be redundancy in the filtered contigs
  --contigs [PATH] If the assembly has already been done you can specify the path to the fasta file containing contigs. If provided the assembly steps will be skipped
  --refs_info [PATH] If you use a simulated dataset and specify the --contigs option. Give the path to the sum_contigs_length_per_annotation.tsv file contained in [out]/assembly/

    CONTIGS ANNOTATION :

  --blast_db [PATH] If you use simulated data. Path to the BLAST database containing reference sequences (HAVE TO be computed before running the pipeline)
  --coverage [INT] Coverage threshold used to filter alignments of contigs on reference genomes or a public database (Default 90)
  --identity [INT] Identity threshold used to filter alignments of contigs on reference genomes or a public database (Default 95)
  --mismatch [INT] Number of mismatch allowed in the seed aligment of BLAST (Default 1)
  --evalue [INT] E-value used for BLAST (Default 10)
  --hit [INT] Maximum number of hits for each querry in the BLAST output (Default 10)
  --link_ref_id_species [PATH] For simulated dataset tab-delimited file containing contigs IDs of reference sequence and the species to which they belong. Used to identify the target in the BLAST of contigs against the references
  --contigs_annotation [PATH] If you use simulated data and have specified the --contigs option specify the path to the (bh_blast_contigs_on_refs.tsv)(to be updated) file in [out]/assembly/

    MAPPING :

  --cpumapping [INT] Number of cpus used for mapping reads on the assembly (Default 4)
  One mapping per sample
  --memmapping [INT] Quantity of RAM in Mb for the mapping of each sample (Default 5000 )
  --bowtie2_indexDir [PATH] Directory were will be stored bowtie2 index of assembled sequences (Default [out]/bowtie2_index)
  --index_prefix [STRING] Prefix for the index files generated by bowtie2-build
  --bamDir [PATH] Directory were will be stored sorted BAM files generated by the mapping and if you use MetaGen indexed BAM files (Default [out]/mapping/bam). If there are sorted BAM files in the folder the mapping step will be skipped and they will be taken as input for the binning programs
  --count_matrix [PATH] If you use Canopy and have specified the --bamDir option. Provide the path to the count_matrix file contained by default in the [out]/mapping/comptage folder

    BINNING :

  --cpubinning [INT] Number of cpus for binning programs (Default same value as --cpus option)
  --nb_cluster [INT] Needed if you use Concoct. Rougth estimation of the expected number of species in the samples. This value will be a starting point for the program that will then refine it
  --bic_step [INT] MetaGen parameter corresponding to the step for the search of the optimal number of clusters (Default 5)
  --auto_method [INT] MetaGen parameter can be either 1 or 2. Recommended to be 2 for large datasets and is the default here

    BINNING EVALUATION :

  --min_bin_size [INT] Minimum size of bins in base (sum of the length of the sequences it contains) to be places in plots (Default 500000)
  --conta_threshold [FLOAT] Maximum contamination percent for good quality bins [0-1] (Default 0.1)
  --comp_threshold [FLOAT] Minimum completeness percent for good quality bins [0-1] (Default 0.6)
```
## Examples
For real data starting from raw reads or reads filtered from contaminant and co-assembly with Spades.
```
  ~/nextflow -c ~/let-it-bin/nextflow_slurm_singularity_common.config run -w [directory to store temporary files] let-it-bin.nf --reads ~/data/reads --out ~/results --cpus 4 --metabat --canopy --maxbin --index_prefix spades_contigs --tmp_checkm tmp
```
For simulated data starting from raw reads and co-assembly with megahit
```
  ~/nextflow -c ~/let-it-bin/nextflow_slurm_singularity_common.config run -w [directory to store temporary files] let-it-bin.nf --reads ~/data/reads --out ~/results --sim_data T --cpus 4
  --nb_ref 50 --metabat2 --cocacola --metagen --mode megahit
  --blast_db ~/blast_db/refs_seq --link_ref_id_species ~/link_id_species.tsv
  --index_prefix spades_contigs
```
