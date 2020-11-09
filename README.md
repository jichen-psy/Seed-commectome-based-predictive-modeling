# Seed-commectome-based-predictive-modeling

This project aims to realize seed commectome-based predictive modeling, i.e., predicting continuous variables using seed-to-whole-brain resting-state functional connectivity (rsFC) patterns. Ten-fold cross-validation is used to assess the out-of-sample prediction performance. A supervised feature extraction procedure is incorporated within cross-validation by performing general linear models (GLMs) on (only) the training set to detect clusters whose rsFC with the seed regions are significantly associated with the target variable to be predicted. Significant clusters are then used as masks to extract the rsFC values from the seed-to-whole-brain rsFC maps for both the training and the test sets. The extracted rsFC values are employed as features for training the predictive models within the traning set via a sparse Bayesian learning approach of relevance vector machine (RVM). Finally, the resulting RVM weights are applied to the rsFC values of test set to obtain the out-of-sample predicted individual scores. 

When the dataset includs multiple sites, leave-one-site-out cross-validations are encouraged to assess the generalizability of predictions across sites. As here a supervised feature selection is implemented, a validation sample that is independent of this feature selection procedure is required in order to avoid the potential issue of cicularity as discussed in Kriegeskorte et al., 2009. 

Supervised feature selection is a long-standing method described in machine learning literature (Guyon and Elisseeff, 2003; Kohavi and John, 1997), and has been proposed to effectively improve predictive performance while discarding irrelevant variables (Finn et al., 2015; Shen et al., 2017; Beaty et al., 2018). Our predictive modeling framework modifies the original “connectome-based predictive modeling” (CPM) introduced in Shen et al., 2017, and incorporates cluster-level inference corrected for multiple comparisons to derive significant clusters that are significantly associated with the target variable. 

##
The statistical parametric mapping (SPM12; https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) is used to perform the GLMs, and the REST package is additionally required for the feature selection procedure which can be freely downloaded at http://www.restfmri.net/forum/REST_V1.8.

Our predictive modeling is based on a relevance vector machine (RVM) (Tipping, 2001) as implemented in the SparseBayes package (http://www.miketipping.com/index.htm). The RVM refers to a specialization of the general Bayesian framework which has an identical functional form to the support vector machine. To implement our predictive modeling with different validation schemes (i.e., 10-fold cross-validation, leave-one-site-out cross-validation, as well as training models within-sample and then validate in independent samples), please add the downloaded SparseBayes package to your Matlab working directory.

## References
Guyon I, Elisseeff A. An introduction to variable and feature selection. Journal of machine learning research 2003; 3: 1157-1182. 

R. Kohavi and G. John. Wrappers for feature selection. Artificial Intelligence 1997: 97(1-2): 273–324. 

Shen X, Finn ES, Scheinost D, et al. Using connectome-based predictive modeling to predict individual behavior from brain connectivity. Nature protocol 2017; 12: 506-518.
 
Tipping, ME. Sparse Bayesian learning and the relevance vector machine. Journal of machine learning research 2001; 1: 211-244. 

Kriegeskorte N, Simmons WK, Bellgowan PSF, et al. Circular analysis in systems neuroscience: the dangers of double dipping. Nature neuroscience 2009; 12(5): 535.

Finn ES, Shen X, Scheinost D, et al. Functional connectome fingerprinting: identifying individuals using patterns of brain connectivity. Nature neuroscience 2015; 18: 1664-1671.
 
Beaty RE, Kenett YN, Christensen AP, et al. Robust prediction of individual creative ability from brain functional connectivity. Proc Natl Acad Sci USA 2018; 115: 1087-1092.
