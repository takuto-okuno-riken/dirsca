%%
% get ROI time-series (matrix) from NIfTI 4D volume.
% returns ROI time-series (X)
% input:
%  V            nifti 4D volume (X x Y x Z x frames)
%  atlasV       nifti 3D atlas (X x Y x Z)
%  operation    operation for each plane ('mode'(default),'max','min','mean','median','sum')

function X = getRoiTSFromNifti4D(V, atlasV, operation)
    if nargin < 3, operation = 'mode'; end
    
    roiIdx = unique(atlasV);
    roiIdx(roiIdx==0) = []; % remove 0

    X = single(zeros(length(roiIdx),size(V,4)));
    A = reshape(V,[],size(V,4));
    for i=1:length(roiIdx)
%    parfor i=1:roiMax
        j = roiIdx(i);
        B = A(atlasV(:)==j,:);
        switch(operation)
        case 'mode'
            m = mode(B,1);
        case 'min'
            m = nanmin(B,[],1);
        case 'max'
            m = nanmax(B,[],1);
        case 'mean'
            m = nanmean(B,1);
        case 'median'
            m = nanmedian(B,1);
        case 'sum'
            m = nansum(B,'all');
        end
        X(i,:) = m;
    end
end

