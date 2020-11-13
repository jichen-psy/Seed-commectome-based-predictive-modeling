function results = SeedBased_Pred_kFold(X, y, k, NumRep, Allsub, CovCatIdx, CovSiteIdx, CovName, datadir, writedir, SPMdir, ROIName, ROIMaskdir)

%% INPUT 
% X          : A matrix with n rows (observations/subjects) and d columns including both the dependent (to be predicted) and confounding variables (to be adjusted). 
%              Both the dependent variable and the predictors (i.e., the extracted seed-to-whole-brain resting-state functional connectivity [rsFC]) are adjusted 
%              according to the current recommendations in Pervaiz et al.(2020).
%              Please ensure that the last column is the dependent variable. 
% y          : A vector with the dependent (target) variable to be predicted (the same as the last column in X)
% Allsub     : A cell from the lookup file containging the names for all subjects in a column
% k          : Number of folds
% NumRep     : Number of repeats 
% CovName    : A cell with one row indicating the name for each variable in X
% CovCatIdx  : A vector indicates which columns are categorical variables. Example: CovCatIdx=[1,2]
% CovSiteIdx : Numeric, indicating which column encodes the site (Example: CovSiteIdx=[2]). Numbering the site column, if there is only one site, just put ones 
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

%% OUTPUT (a structure with results including performance measures)
% results.
%        R: correlation between the target variable and their predicted
%        scores for each fold and repeat
%        nRMSE: normalized root-mean-square-error between the target variable 
%        and their predicted scores for each fold and repeat
%        YPredAll: out-of-sample predicted scores
%        YTrueAll: confound-adjusted true scores
%        R_all:correlation between the target variable and their predicted
%        scores per repeat
%        nRMSE_all: normalized root-mean-square-error between the target
%        variable and their predicted scores per repeat

% ---References
% Pervaiz U, Vidaurre D, Woolrich MW, Smith SM (2020): Optimising network 
% modelling methods for fMRI. NeuroImage 211, 116604.

% Tipping ME (2001): Sparse Bayesian learning and the relevance vector machine. 
% J Mach Learn Res 1: 211-244.

%   Last edited by Ji Chen on 17-Oct-2019 
%%
    
%Generate stratified folds
  sites = unique(X(:,CovSiteIdx));
  indices = {};
  for i=1:length(sites)
      indices{i} = find(X(:,CovSiteIdx)==sites(i));
  end

  for ith_repeat = 1:NumRep
     for i=1:length(sites)
         ind = indices{i};
         ith_partition = cvpartition(length(ind),'KFold',k);
         for ith_fold=1:k
           j = ith_partition.test(ith_fold);
           CVindices(ind(j),ith_repeat) = ith_fold;
         end
     end
  end


for ith_repeat=1:NumRep
  
  pea = nan(k,1);
  ypred_all = nan(size(X,1),1);
  yres_All= nan(size(X,1),1);
 
 for ith_fold=1:k 
    
    TestIdx=CVindices(:,ith_repeat)==ith_fold;
    TrainIdx=CVindices(:,ith_repeat)~=ith_fold; 
    Trainy=y(find(TrainIdx==1),1);   
    Testy = y(find(TestIdx==1),1);
    
    %feature extraction (only within the training samples)
    [FCAll, InddAll]= Extract_Feature_kFold(TrainIdx,Allsub,X,CovName, CovSiteIdx, ROIName, ROIMaskdir, SPMdir, datadir, writedir);
    if isempty(FCAll)
    continue;
    end
        
   model = [];
   model.hyperparameters = [];
   DesignMatrix=[];
   DesignMatrix=x2fx(X(:,1:end-1),'linear', CovCatIdx);

   %Confound adjustment
   [Trainy,reg_y] = regress_confounds(Trainy, DesignMatrix(TrainIdx, :));
   [FCAll, reg_x] = regress_confounds(FCAll, DesignMatrix(TrainIdx, :));

   [n,d] = size(FCAll);
   
   % RVM
   [model.rvm, model.hyperparams, model.diagnostics] = SparseBayes('Gaussian', [FCAll, ones(n,1)], Trainy);

   % get the RVM weights
   model.weights = zeros(d,1);

   model.b = 0;
   if model.rvm.Relevant(end)==(d+1) 
     model.b = model.rvm.Value(end);
     model.weights(model.rvm.Relevant(1:end-1)) = model.rvm.Value(1:end-1);
   else
     model.weights(model.rvm.Relevant) = model.rvm.Value;
   end 
   
   %Extract the features for the test set using the masks determined on the training set
   TestX=Apply_Masks(TestIdx,Allsub, ROIName, datadir, InddAll);
   
   %Confound adjustment
   TestX = regress_confounds(TestX, DesignMatrix(TestIdx, :), reg_x);
   Testy = regress_confounds(Testy, DesignMatrix(TestIdx, :), reg_y);
   
   %Get the predicted scores
   ybar.y = TestX*model.weights + model.b;   
   
   ypred_all(TestIdx) = ybar.y;
   yres_All(TestIdx)=Testy;
   
   % validation set performance
   pea(ith_fold) = corr(Testy, ybar.y,'type', 'Pearson', 'Rows', 'complete');
   rmse = sqrt(sum(((Testy - ybar.y).^2))./(length(ybar.y)-1));
   nrmse = rmse./sqrt(sum(((Testy-mean(Testy)).^2))./(length(ybar.y)-1));
   nrmse(ith_fold)=nrmse;
  end
 
  nRMSE(:,ith_repeat) =nrmse;
  Pearson(:,ith_repeat)  = pea;
  ypredicted(ith_repeat,:) = ypred_all;
  yRes(ith_repeat,:) = yres_All;

  % overall performance
  Pearson_all(ith_repeat) = corr(yres_All, ypred_all,'type', 'Pearson', 'Rows', 'complete');
  rmse = sqrt(sum(((yres_All - ypred_all).^2))./(length(ypred_all)-1));
  nrmse = rmse./sqrt(sum(((yres_All-mean(yres_All)).^2))./(length(ypred_all)-1));
  nRMSE_all(ith_repeat)=nrmse;

 end 


%% create a single results structure
results= [];

% validation
results.YPredAll = ypredicted;
results.YTrueAll=yRes;
results.R  = Pearson;
results.nRMSE =nRMSE;
results.R_all  = Pearson_all;
results.nRMSE_all = nRMSE_all;
end


    
    
    
    