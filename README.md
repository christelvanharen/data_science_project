# Genomic polymorphisms in ribo-interactome may underlie variation in gene expression in people with cancer

The aim of this study is to demonstrate by means of an expression
quantitative trait loci (eQTL) analysis an association between SNPs in ribosomal and ribosome-associated proteins, and a difference in gene expression between cancer and non-cancer patients.
During the eQTL analysis, gene expression is compared by regression to SNPs in the ribo-interactome of individuals. This allows the effect of these SNPs on gene expression to be determined. For this research, a pipeline was created with which the analysis was performed for test data of the human chromosome 10.

## Clone the project
First, start by cloning this project using the following command:
git clone http://bio-inf.han.nl:3000/Bio-Inf_Jaar3_2122/BIN-3f.git

## Installation
Make sure all the packages from the requirements.txt have been
installed. STAR, GATK and Snakemake are all tools which needs to be installed on the server/PC as well. 

## SnakeMake
There are 2 ways to use the Snakemake. Firstly, if you are running SnakeMake on the server provided by the HAN and used by BIN-3f, then you only have to navigate to the folder where the Snakefile is located. In this case that will be /home/administrator/pipeline_final. Then you can just use the following command: Snakemake -F, or if you want a specific file to be generated: Snakemake -F "outputfile.txt". Note: The "outputfile.txt" is an example of a random outputfile, in the Snakefile there is no such file as "outputfile.txt" and you'll have to look in the Snakefile for an outputfile.

The second way of running the SnakeMake is by downloading Snakemake yourself and running it on Linux. The command stays the same. However one of the very important things to do is when you're running this Snakefile on a different server/PC then you'll need to change the directories of the input and outputfiles. This needs to be done in the Snakefile provided in this GIT. 

## STAR
To be able to run STAR it is highly necessary to download the 12 FastQ files of the 6 patients provided by STAR itself. 

