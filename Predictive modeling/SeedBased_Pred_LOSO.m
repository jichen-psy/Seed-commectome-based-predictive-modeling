
function results = SeedBased_Pred_LOSO(X, y, Allsub, CovCatIdx, CovSiteIdx, CovName, datadir, writedir, SPMdir, ROIName, ROIMaskdir)


% This script may be used if you want to test the generalizability of your
% model across sites, as the leave-one-site-out (LOSO) cross-validation, when your 
% data have multiple sites, i.e., training the model on all sites but one and testing on the left-out site

%%INPUT
% X          : A matrix with n rows (observations/subjects) and d columns including both the dependent (to be predicted) and confounding variables (to be adjusted). 
%              Both the dependent variable and the predictors (i.e., the extracted seed-to-whole-brain resting-state functional connectivity [rsFC]) are adjusted 
%              according to the current recommendations in Pervaiz et al.(2020).
%              Please ensure that the last column is the dependent variable. 
% y          : A vector with the dependent (target) variable to be predicted (the same as the last column in X)
% Allsub     : A cell from the lookup file containging the names for all subjects in a column
% CovName    : A cell with one row indicating the name for each variable in X
% CovCatIdx  : A vector indicates which columns are categorical variables. Example: CovCatIdx=[1,2]
% CovSiteIdx : Numeric, indicating which column encodes the site (Example: CovSiteIdx=[2]). Numbering the site column, if there is only one site, then no need to run LOSO cross-validation. 
% datadir    : A cell indicates the directory containing each subject's seed-to-whole-brain rsFC maps for each seed separately, 
%              the number of datadir equals to the number of seeds e.g., three seed: a 1x3 cell.
% writedir   : A cell indicates the directory to store the output images, for each seed separately; e.g., three seed: a 1x3 cell.
% SPMdir     : a 1x1 cell indicates the working director for SPM
% ROIMaskdir : A cell indicates the mask files (i.e., the masks within which the seed-to-whole-brain rsFC were calculated). 
%              The mask for each seed can be created by substracting the seed region (binarized) from a whole-brain grey matter binary mask.
% ROIName    : A cell contains the names you have defined for the seed-to-whole-brain rsFC
%              maps. Example: assuming the name of the rsFC map for one subject and a
%              particular seed is sub001_STG.nii, then you please define the
%              ROIName{1}='_STG'.

% OUTPUT (a structure with results including performance measures)
% results.
%        R: correlation between the target variable and their predicted scores 
%        nRMSE: the normalized root-mean-square-error between the target variable and their predicted scores 
%        YPredAll: out-of-site predicted scores
%        YTrueAll: confound-adjusted true scores

% ---References
% Pervaiz U, Vidaurre D, Woolrich MW, Smith SM (2020): Optimising network 
% modelling methods for fMRI. NeuroImage 211, 116604.

% Tipping ME (2001): Sparse Bayesian learning and the relevance vector machine. 
% J Mach Learn Res 1: 211-244.

%   Last edited by Ji Chen on 17-Oct-2019 

%%
CVindices=X(:,CovSiteIdx);
NumSites=length(unique(CVindices));

  yPred_all = nan(size(y,1),1);
  yTrue_all= nan(size(y,1),1);

    
   for ith_site=1:NumSites 
    
TestIdx=CVindices(:,1)==ith_site;
TrainIdx=CVindices(:,1)~=ith_site;
    
    Trainy=y(find(TrainIdx==1),1);
    Testy = y(find(TestIdx==1),1);
    
    [FCAll, InddAll]= Extract_Feature_LOSO(TrainIdx,Allsub,X,CovName, CovSiteIdx, NumSites, ROIName, ROIMaskdir, SPMdir, datadir, writedir);
   if isempty(FCAll)
    continue;
   end
   model = [];
   model.hyperparameters = [];
   NumCov=1:(size(X,2)-1);
   NumCov(:,CovSiteIdx)=[];
   CovCate=intersect(NumCov,CovCatIdx);
   for Cate=1:length(CovCate)
   if CovCate(1,Cate)>CovSiteIdx
       CovCate(1,Cate)=CovCate(1,Cate)-1;
   end
   end
   DesignMatrix=[]; DesignMatrixTr=[];
   DesignMatrixTr=x2fx(X(TrainIdx,1:end-1),'linear',CovCatIdx); 
   DesignMatrix=x2fx(X(:,NumCov),'linear',CovCate);
 
    [~, reg_y] = regress_confounds(Trainy, DesignMatrix(TrainIdx, :)); 
    Trainy = regress_confounds(Trainy, DesignMatrixTr);

    [~, reg_x] = regress_confounds_x(FCAll, DesignMatrix(TrainIdx, :));
    FCAll = regress_confounds_x(FCAll, DesignMatrixTr);


   [n,d] = size(FCAll);
   
   % RVM
   [model.rvm, model.hyperparams, model.diagnostics] = SparseBayes('Gaussian', [FCAll, ones(n,1)], Trainy);

   % get the weights
   model.weights = zeros(d,1);

   model.b = 0;
   if model.rvm.Relevant(end)==(d+1) 
     model.b = model.rvm.Value(end);
     model.weights(model.rvm.Relevant(1:end-1)) = model.rvm.Value(1:end-1);
   else
     model.weights(model.rvm.Relevant) = model.rvm.Value;
   end 
    
   TestX=Apply_Masks(TestIdx,Allsub, ROIName, datadir, InddAll);  
    
   %Confound adjustment
   TestX = regress_confounds(TestX, DesignMatrix(TestIdx, :),reg_x);
   Testy = regress_confounds(Testy, DesignMatrix(TestIdx, :),reg_y);

   ybar.y = TestX*model.weights + model.b;   
   
  
   yPred_all(TestIdx) = ybar.y;
   yTrue_all(TestIdx)=Testy;

   end

  
%% create a single results structure

results = [];

% validation performance
results.yPredicted = yPred_all;
results.yTrue= yTrue_all;
rmse_all = sqrt(sum(((yTrue_all - yPred_all).^2))./(length(yPred_all)-1));
results.nRMSE = rmse_all./sqrt(sum(((yTrue_all-mean(yTrue_all)).^2))./(length(yPred_all)-1));
results.R=corr(yPred_all,yTrue_all,'type', 'Pearson', 'Rows', 'complete');  

end


    
    
    
    