# BRF-sRNA-target
This repository contains source code and accessory files for a Random Forest model for the prediciton of non-coding RNA targets in bacteria, described in Mitra and Kushner, 2015 "Supervised Prediction of Regulatory Non-coding RNA Targets in Bacteria Using Alignment-independent Sequence Information". 

The R code in this repo utilizes two packages: package “randomForest” for training or using the models, and package “seqinr” for parsing FASTA input files. To install the packages, run the following commands from the R console:

> install.packages(“randomForest”,dependencies=TRUE)

> install.packages(“seqinr”,dependencies=TRUE)

####Training the BRF model
The file “BRF_modelTrain.r” contains the code for training the Balanced Random Forest model used in the study. To allow for the sampling procedures used in the algorithm, the training set is divided into the positive and negative sets in the files “posTrain.txt” and “negTrain.txt”. The script outputs the OOB results for the model, and saves the model in the current working directory as “BRF_model.Rdata”.

####Making new predictions
To make new predictions using the tuned BRF model, run the R script in “BRF_predict.r”. This program uses the files containing the kmer patterns and their names for sRNA (“srnaPatterns.txt”) and mRNA (“mrnaPatterns.txt”). The tuned model used for the predictions is stored in “BRF_model.Rdata”. The user needs to ensure that these files are present in the working directory.

The user input is in the form of two FASTA formatted files, one for the sRNA sequences (“srna_fasta.txt”) and the other for their corresponding mRNA sequences (“mrna_fasta.txt”). The order of the pairs of sequences must be maintained for interaction to be predicted. 

The program outputs probability values (random forest vote ratios) for each instance in the input.
