# NRFSE
* We propose a neighborhood-regularization method (NRFSE) that leverages multi-view data on drugs and side effects to predict the frequency of side effects.<br>
* First, we adopt a class-weighted non-negative matrix factorization to decompose the drug-side effect frequency matrix, in which Gaussian likelihood is used to model unknown drug-side effect pairs.<br>
* Second, we design a multi-view neighborhood regularization to integrate three drug attributes and two side effect attributes, respectively, which makes most similar drugs and most similar side effects have similar latent signatures. The regularization can adaptively determine the weights of different attributes.
* More details can be found in our paper (A neighborhood-regularization method leveraging multi-view data for predicting the frequency of drug side effects, published in *Bioinformatics*).


# Environment requirements
Matlab R2019b or later version

# Files
## Source files
* main.m--
* warm_test.m -- to perform 10-fold cross-validation under warm-start scenario.
* DecompositionAlgorithm_wei.m
* KNN_fun.m -- 
* auc.m -- to compute AUC and AUPR values.
## Data sets


# Getting started
1.Unarchive the data.rar and put the dataset and source code in the same directory.<br>
2.Run the main.m to perform 10-fold cross-validation under warm-start scenario.
