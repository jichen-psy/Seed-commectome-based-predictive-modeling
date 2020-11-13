
function results = SeedBased_Pred_Validation(y, Confound, Allsub, Clust, GroupIdx, CovCatIdx, CovSiteIdx,ROIName,datadir) 

% For validation analysis: assuming you have two independent datasets,
% then you train a RVM model within one dataset and then apply it to the
% other dataset for validation. The training and the validation datasets
% used for this analysis should have the same number of confounds 

% y          :A vector with the dependent (continuous) variable to be predicted for both the datasets
% Allsub     :A cell contains the names as one column for all subjects 
% Clust      :A cell indicates the directory to the identified clusters (in NIFTI format)
%             per seed for extracting the resting-state functional connectivity (rsFC) values (i.e., the rsFC patterns you would like to validate).
%             For example, if you have three seeds, and each seed has revealed 2 significant/robust clusters,
%             then it should be a 1x3 cell array, and each cell contains two cells with each cell denotes one cluster 
% datadir    :A cell structure indicates the directory where the seed-to-whole-brain rsFC maps for each seed are stored  
% GroupIdx   :A numeric column array indexing which subjects are the training
%             sample and which are from the validation dataset
%             Example: GroupIdx=[1;1;1;1;1;2;2;2;2;2]; 
% Confound   :A matrix with confounding variables that would like to be adjusted, here both the
%             dependent variable and the predictors are adjusted according to the current
%             recommendations in Pervaiz et al. (2020)
% CovCateIdx :A vector indicates which columns in your confound matrix are categorical variables, 
%             Example: CovCatIdx=[1,2]
% CovSiteIdx :Numeric,indicating which column encodes the site in your confound matrix
%             Example: CovSiteIdx=[2]

% OUTPUT (a structure with results including performance measures)
% results.
%        R & P: correlation coefficient between target variable and their
%        predicted scores as well as the p value
%        YPred: predicted scores in the validation sample
%        YTrue: confound-adjusted true scores

% ---References
% Pervaiz U, Vidaurre D, Woolrich MW, Smith SM (2020): Optimising network 
% modelling methods for fMRI. NeuroImage 211, 116604.

% Tipping ME (2001): Sparse Bayesian learning and the relevance vector machine. 
% J Mach Learn Res 1: 211-244.

%----Last edited by Ji Chen on 17-Oct-2019 

prompt = {'Enter a number you encode the training sample in "GroupIdx"'};
dlgtitle = 'GroupIndex';
TrainNum = inputdlg(prompt,dlgtitle);

prompt = {'Enter a number you encode the independent validation sample in "GroupIdx"'};
dlgtitle = 'GroupIndex';
ValNum = inputdlg(prompt,dlgtitle);
TrainIdx=GroupIdx==str2num(TrainNum{1});
ValIdx=GroupIdx==str2num(ValNum{1});
TrainSub=Allsub(TrainIdx,:);
ValSub=Allsub(ValIdx,:);

yTrain=y(TrainIdx,:);
yVal=y(ValIdx,:);
CovTr=Confound(TrainIdx,:);
CovVal=Confound(ValIdx,:);

AllSeedFCTrain=[]; AllSeedFCVal=[];
%Extract features from the training dataset and build the RVM model
for Seed=1:length(Clust)
SeedFCTrain=[]; SeedFCVal=[];
for Cluster=1:length(Clust{Seed})
MASK=spm_vol(Clust{Seed}{Cluster});
MASKV=spm_read_vols(MASK);
ind=find(MASKV==1);
[maskXYZ(:,1), maskXYZ(:,2), maskXYZ(:,3)] = ind2sub(MASK.dim,ind);
for s= 1:numel(TrainSub) 
Cfil = dir(fullfile(datadir{Seed},strcat(TrainSub{s},ROIName{Seed})));
FCgm = spm_vol(fullfile(datadir{Seed},Cfil.name));
FC=spm_sample_vol(FCgm,maskXYZ(:,1),maskXYZ(:,2),maskXYZ(:,3),0);
TrainFC (s,Cluster) = nanmean(FC);

end
clear FC FCgm Cfil
for s= 1:numel(ValSub) 
Cfil = dir(fullfile(datadir{Seed},strcat(ValSub{s},ROIName{Seed})));
FCgm = spm_vol(fullfile(datadir{Seed},Cfil.name));
FC=spm_sample_vol(FCgm,maskXYZ(:,1),maskXYZ(:,2),maskXYZ(:,3),0);
ValFC (s,Cluster) = nanmean(FC);

end
SeedFCTrain=[SeedFCTrain,TrainFC];
SeedFCVal=[SeedFCVal,ValFC];
clear MASK MASKV ind maskXYZ FC FCgm Cfil TrainFC ValFC
end
AllSeedFCTrain=[AllSeedFCTrain,SeedFCTrain];
AllSeedFCVal=[AllSeedFCVal,SeedFCVal];
clear SeedFCTrain SeedFCVal
end
DesignMatrix1=x2fx(CovTr, 'linear',CovCatIdx);
[y_train, reg_y] = regress_confounds(yTrain, DesignMatrix1);
[x_train, reg_x] = regress_confounds(AllSeedFCTrain, DesignMatrix1);

[nn,d]=size(x_train);
[model.rvm, model.hyperparams, model.diagnostics] = SparseBayes('Gaussian', [x_train, ones(nn,1)], y_train);
 model.weights = zeros(d,1);

   model.b = 0;
   if model.rvm.Relevant(end)==(d+1) 
     model.b = model.rvm.Value(end);
     model.weights(model.rvm.Relevant(1:end-1)) = model.rvm.Value(1:end-1);
   else
     model.weights(model.rvm.Relevant) = model.rvm.Value;
   end 

NumCov=1:size(CovTr,2);
NumCov(:,CovSiteIdx)=[];
CovCate=intersect(NumCov,CovCatIdx);
 for Cate=1:length(CovCate)
   if CovCate(1,Cate)>CovSiteIdx
       CovCate(1,Cate)=CovCate(1,Cate)-1;
   end
 end
   
CovTr(:,CovSiteIdx)=[];
DesignMatrix2=x2fx(CovTr, 'linear',CovCate);
[~, reg_y1] = regress_confounds(yTrain, DesignMatrix2);
[~, reg_x1] = regress_confounds(AllSeedFCTrain, DesignMatrix2);

if unique(CovVal(:,CovSiteIdx))~=1
DesignMatrixVal=x2fx(CovVal(:,CovSiteIdx), 'linear',1);
[y_val, ~] = regress_confounds(yVal, DesignMatrixVal);
[x_val, ~] = regress_confounds(AllSeedFCVal, DesignMatrixVal);
else
    y_val=yVal;
    x_val=AllSeedFCVal;
end
CovVal(:,CovSiteIdx)=[];
[y_val,~] = regress_confounds(y_val, x2fx(CovVal, 'linear'), reg_y1);
x_val = regress_confounds(x_val,  x2fx(CovVal, 'linear'), reg_x1);

ypred_val = x_val*model.weights + model.b; 

results.YPred=ypred_val;
results.YTrue=y_val;
[results.R,results.P]=corr(ypred_val,y_val, 'type', 'Pearson', 'Rows', 'complete');
   
   







