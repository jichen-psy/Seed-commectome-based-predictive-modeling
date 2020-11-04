# Seed-commectome-based-predictive-modeling

This by incorporating a supervised feature extraction procedure only within the training samples, and then take the rsFC of clusters showing significant associations with the target variable (continuous scores) as the features to train a model to predict individual scores in the test set.


To SPM, and REST to run the feature selection procedure which can be downloaded at https://www.fil.ion.ucl.ac.uk/spm/software/spm12/ and http://www.restfmri.net/forum/REST_V1.8.


Our predictive modelingnode is based on a relevance vector machine (RVM) (Tipping, 2001) as implemented in the SparseBayes package (http://www.miketipping.com/index.htm). The RVM refers to a specialization of the general Bayesian framework which has an identical functional form to the support vector machine. To implement our predictive modeling with different validation schemes (i.e., 10-fold cross-validation, leave-one-site-out cross-validation, as well as training models within-sample and then validate in independent samples), please add the downloaded SparseBayes package to your Matlab working directory.
