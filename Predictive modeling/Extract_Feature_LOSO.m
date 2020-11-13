
function [AllFCROI, InddAll]=Extract_Feature_LOSO(Index,Allsub,Cov, CovName, CovSiteIdx, NumSites, ROIName, ROIMaskdir, SPMdir, datadir, writedir)

%   Last edited by Ji Chen on Oct-2020 

INDEX=find(Index==1); TrainSub=cell(0);
for i=1:numel(INDEX)
TrainSubX=Allsub{INDEX(i),1};
TrainSub=[TrainSub, TrainSubX];
end

TrainSub=TrainSub';
TrainCov=Cov(INDEX,1:end);

if min(TrainCov(:,CovSiteIdx))~=1
    TrainCov(:,CovSiteIdx)=TrainCov(:,CovSiteIdx)-1;
else
    U1=setdiff(1:NumSites,unique(TrainCov(:,CovSiteIdx)));
    if U1~=NumSites
    QW=find(TrainCov(:,CovSiteIdx)>U1);
    QW1=TrainCov(:,CovSiteIdx);
    QW1(QW,1)=QW1(QW,1)-1;
    TrainCov(:,CovSiteIdx)=QW1;
    end
        
end

TrainSiteIdx=TrainCov(:,CovSiteIdx);
NewCovIdx=setdiff(1:size(Cov,2),CovSiteIdx);
AllFCROI=[]; 
for ROI=1:length(ROIMaskdir)
    
for j=1:length(unique(TrainSiteIdx))
    a=find(TrainSiteIdx==j);
    a=a';
for i=a
fils1(i) = dir(fullfile(datadir{ROI},strcat(TrainSub{i},ROIName{ROI}))); 
end
fils1=fils1';
filsAll{j}=fils1(a,1);
clear fils1 a
end

%%Call the SPM functions to derive the statistical images for the seed regions
matlabbatch={};
matlabbatch{1}.spm.stats.factorial_design.dir = {SPMdir};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac.name = 'site';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac.dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac.variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac.gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac.ancova = 0;


for j=1:length(unique(TrainSiteIdx))
  file=filsAll{j};  
    aa=find(TrainSiteIdx==j);
    aa=aa';
for file1 = 1:size(file,1)
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(j).scans{file1,1} = fullfile(datadir{ROI},file(file1).name);

matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(j).conds = TrainCov(aa,CovSiteIdx);

end
clear file aa
end

matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.fmain.fnum = 1;
for i=1:length(NewCovIdx)
matlabbatch{1}.spm.stats.factorial_design.cov(i).c = TrainCov(:,NewCovIdx(i));

matlabbatch{1}.spm.stats.factorial_design.cov(i).cname = CovName{i};
matlabbatch{1}.spm.stats.factorial_design.cov(i).iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov(i).iCC = 1;
end
%%

matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 0;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {strcat(ROIMaskdir{ROI},',1')};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
spm_jobman('run',matlabbatch);

%% run estimate
matlabbatch ={};
matlabbatch{1}.spm.stats.fmri_est.spmmat = {strcat(SPMdir,'/SPM.mat')};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);

%% define contrasts 
matlabbatch = {};
matlabbatch{1}.spm.stats.con.spmmat = {strcat(SPMdir,'/SPM.mat')};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Neg';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [zeros(1,length(unique(TrainSiteIdx))),zeros(1,length(NewCovIdx)-1),-1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;

matlabbatch{1}.spm.stats.con.spmmat = {strcat(SPMdir,'/SPM.mat')};
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Pos';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [zeros(1,length(unique(TrainSiteIdx))),zeros(1,length(NewCovIdx)-1),1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);

load(strcat(SPMdir,'/SPM.mat'))
[K1,~] = CorrClusTh(SPM,0.001,0.05);
InddF=cell(0); AllFCA=[];
for Contrast=1:2
[data,~,head]=rest_readfile(strcat(SPMdir,strcat('/','spmT_000',num2str(Contrast), '.nii')));
A=[]; Index1=[]; b=[]; 
data1=data(:);
Ind1=find(data1~=0);
Ind2=1-(tcdf(data1(Ind1,1),SPM.nscan-8));
data1(Ind1,1)=Ind2;
data1(find(data1>=0.001),1)=0;
data1(find(data1<0.001&data1>0),1)=1;
data=reshape(data1,head.dim);
nozeropos=find(data~=0);
[i j k]=ind2sub(size(data),nozeropos);
cor=[i j k];
FC1=[];
L=cor';
dim = [max(L(1,:)) max(L(2,:)) max(L(3,:))];
vol = zeros(dim(1),dim(2),dim(3));
indx = sub2ind(dim,L(1,:)',L(2,:)',L(3,:)');
vol(indx) = 1;
[cci,~] = bwlabeln(vol,18);
A = cci(indx');

clusterID = unique(A);

for i=1:numel(clusterID)
    Z(1,i)= length(find(A(1,:)==i));
end

Index1=find(Z > K1); 
MASK=spm_vol(ROIMaskdir{ROI});
MASKV=spm_read_vols(MASK);
ind=find(MASKV>-1);
[maskXYZ(:,1), maskXYZ(:,2), maskXYZ(:,3)] = ind2sub(MASK.dim,ind);

if size(Index1,2)>1 % if more than one cluster is detected to be significant

    for i=1:numel(Index1)
        [~,b]=find(A==Index1(1,i));
        INDD=cor(b,:);
        Indd{i}=INDD;
    end

ForWrite=[1:902629]';
ForWrite(1:902629,1)=0;
Inddd=[];Indddd=[];
for nn=1:numel(Indd)

        XD=Indd{nn};
        for i=1:size(XD,1)
            [aa,~]=find((maskXYZ(:,1)==XD(i,1))&(maskXYZ(:,2)==XD(i,2))&(maskXYZ(:,3)==XD(i,3)));
            Indddd=[Indddd,aa];
            
            for s= 1:numel(TrainSub)
                Cfil = dir(fullfile(datadir{ROI},strcat(TrainSub{s},ROIName{ROI})));
                
                FCgm = spm_vol(fullfile(datadir{ROI},Cfil.name));
                EX=Indd{nn};
                FC=spm_sample_vol(FCgm,EX(:,1),EX(:,2),EX(:,3),0);
                FCEX  = nanmean(FC);
                FC1(s,nn)=FCEX;
            end
            
        end
    Inddd=[Inddd,Indddd];
    clear XD 
    Indddd=[];
end
Inddd=Inddd';
ForWrite(Inddd,1)=1;
    if Contrast==1
    dinfo = dir (strcat(writedir{ROI},'\Neg\'));
    dinfo([dinfo.isdir])=[];
    
    for i=1:numel(dinfo)
        WW=strsplit(dinfo(i).name,'.');
        WWW=WW{1};
        AllWW(i,1)=str2num(WWW);
    end
     MaxWW=max(AllWW);
     MASK.fname= strcat(writedir{ROI},'\Neg\',[num2str(MaxWW+1),'.nii']);
    else
    dinfo = dir (strcat(writedir{ROI},'\Pos\'));
    dinfo([dinfo.isdir])=[];
    
    for i=1:numel(dinfo)
        WW=strsplit(dinfo(i).name,'.');
        WWW=WW{1};
        AllWW(i,1)=str2num(WWW);
    end
     MaxWW=max(AllWW);
     MASK.fname= strcat(writedir{ROI},'\Pos\',[num2str(MaxWW+1),'.nii']);
    end
    ForWriteF=reshape(ForWrite,MASK.dim);
    spm_write_vol(MASK,ForWriteF);

elseif size(Index1,2)==1% if there is one cluster detected as significant
[~,b]=find(A==Index1);
Indd=cor(b,:);   

    Inddd=[];
for i = 1:size(Indd,1) 
    [aa,~]=find((maskXYZ(:,1)==Indd(i,1))&(maskXYZ(:,2)==Indd(i,2))&(maskXYZ(:,3)==Indd(i,3)));
    Inddd=[Inddd,aa];
end
Inddd=Inddd';

ForWrite=[1:902629]';
ForWrite(1:902629,1)=0;

for s= 1:numel(TrainSub)
    
    Cfil = dir(fullfile(datadir{ROI},strcat(TrainSub{s},ROIName{ROI})));
    FCgm = spm_vol(fullfile(datadir{ROI},Cfil.name));
    FC=spm_sample_vol(FCgm,Indd(:,1),Indd(:,2),Indd(:,3),0);
    FC1 (s,1) = nanmean(FC);
end
  
ForWrite(Inddd,1)=1;
    
     if Contrast==1
    dinfo = dir (strcat(writedir{ROI},'\Neg\'));
    dinfo([dinfo.isdir])=[];
    
    for i=1:numel(dinfo)
        WW=strsplit(dinfo(i).name,'.');
        WWW=WW{1};
        AllWW(i,1)=str2num(WWW);
    end
     MaxWW=max(AllWW);
     MASK.fname= strcat(writedir{ROI},'\Neg\',[num2str(MaxWW+1),'.nii']);
    else
    dinfo = dir (strcat(writedir{ROI},'\Pos\'));
    dinfo([dinfo.isdir])=[];
    
    for i=1:numel(dinfo)
        WW=strsplit(dinfo(i).name,'.');
        WWW=WW{1};
        AllWW(i,1)=str2num(WWW);
    end
     MaxWW=max(AllWW);
     MASK.fname= strcat(writedir{ROI},'\Pos\',[num2str(MaxWW+1),'.nii']);
    end
    ForWriteF=reshape(ForWrite,MASK.dim);
    spm_write_vol(MASK,ForWriteF);
         
elseif size(Index1,2)==0 % if there is no cluster detected as significant

    ForWrite=[1:902629]';
    ForWrite(1:902629,1)=0;
    
    if Contrast==1
    dinfo = dir (strcat(writedir{ROI},'\Neg\'));
    dinfo([dinfo.isdir])=[];
    
    for i=1:numel(dinfo)
        WW=strsplit(dinfo(i).name,'.');
        WWW=WW{1};
        AllWW(i,1)=str2num(WWW);
    end
     MaxWW=max(AllWW);
     MASK.fname= strcat(writedir{ROI},'\Neg\',[num2str(MaxWW+1),'.nii']);
    else
    dinfo = dir (strcat(writedir{ROI},'\Pos\'));
    dinfo([dinfo.isdir])=[];
    
    for i=1:numel(dinfo)
        WW=strsplit(dinfo(i).name,'.');
        WWW=WW{1};
        AllWW(i,1)=str2num(WWW);
    end
     MaxWW=max(AllWW);
     MASK.fname= strcat(writedir{ROI},'\Pos\',[num2str(MaxWW+1),'.nii']);
    end
    ForWriteF=reshape(ForWrite,MASK.dim);
    spm_write_vol(MASK,ForWriteF);
    
           FC1=[];  
           Indd=[];

end
if ~isempty(FC1)
    AllFCA=[AllFCA,FC1];
end
 
if ~isempty(Indd)&&iscell(Indd)==1
    InddF=[InddF,Indd];
elseif ~isempty(Indd)&&iscell(Indd)==0
    InddF1=mat2cell(Indd,size(Indd,1));
    InddF=[InddF,InddF1];
end
   
 clear AllWW MaxWW WW WWW maskXYZ MASK MASKV Indd Inddd InddF1 FC FCgm Cfil i j k cci A indx clusterID Z Index1 b dinfo ForWrite EX FCEX ind Index1
end
if ~isempty(AllFCA)
AllFCROI=[AllFCROI,AllFCA];
end

InddAll.(char(ROI+64))=InddF;

which_dir = SPMdir;
dinfo = dir (which_dir);
dinfo([dinfo.isdir])=[];
filenames = fullfile(which_dir,{dinfo.name});
delete(filenames{:})
clear matlabbatch AllFCA InddF
end
