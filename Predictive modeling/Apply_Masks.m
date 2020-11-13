  
function TFCAll=Apply_Masks(Index, Allsub, ROIName, datadir, InddAll)

%---- Last edited by Ji Chen on Nov-2020

%Extract the test subjects
INDEX1=find(Index==1); TestSub=cell(0);
for i=1:size(INDEX1)
    TestSubX=Allsub{INDEX1(i),1};
    TestSub=[TestSub, TestSubX];
end

TestSub=TestSub';
%Extract the rsFC values for the test set using the voxel position (mask)
%indices derived from the training set

    TFCAll=[];
    for ROI=1:length(ROIName)
        if isempty(InddAll.(char(ROI+64)))
            TFC1=[];
        else
            
            for s= 1:numel(TestSub)
                Tfil = dir(fullfile(datadir{ROI},strcat(TestSub{s},ROIName{ROI})));
                for ss=1:numel(InddAll.(char(ROI+64)))
                    TFCgm = spm_vol(fullfile(datadir{ROI},Tfil.name));
                    EX=InddAll.(char(ROI+64)){ss};
                    TFC=spm_sample_vol(TFCgm,EX(:,1),EX(:,2),EX(:,3),0);
                    TFCEX (1,ss) = nanmean(TFC);
                end
                TFC1(s,:)=TFCEX;
            end
        end
        if ~isempty(TFC1)
            TFCAll=[TFCAll,TFC1];
        end
        clear TFC TFCgm TFC1 TFCEX
    end
  
