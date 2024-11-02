% get nuisance time-series of aCompCor (white matter & csf ROI) component
% time-series based on Y.Behzadi et al. (2007).
% returns csf <compNum> components, wm <compNum> components.
% input:
%  V            4D fMRI time-series
%  csfV         mask for CSF
%  wmV          mask for White matter
%  Sd           detrend signal (time-series x node) before PCA (optional)
%  maskTh       mask threshold (percentile) for ROI of aCompCor (default:99)
%  compNum      component number of aCompCor (default:6)

function Xn = getNuisanceaCompCor(V, csfV, wmV, Sd, maskTh, compNum)
    if nargin < 6, compNum = 6; end
    if nargin < 5, maskTh = 99; end
    if nargin < 4, Sd = []; end
    
    % get ROI index for aCompCor
    csfIdx = find(csfV>=prctile(csfV(:),maskTh)); % extract top 1% voxels
    wmIdx = find(wmV>=prctile(wmV(:),maskTh)); % extract top 1% voxels

    % detrend
    V = reshape(V, size(V,1)*size(V,2)*size(V,3),size(V,4));
    M = nanmean(V,2);
    V = V - M;
    csfX = V(csfIdx,:);
    wmX = V(wmIdx,:);
    if ~isempty(Sd)
        [~, ~, perm, RiQ, dR2i] = regressPrepare(Sd);
        % 1st step OLS regression
        for i=1:length(csfIdx)
            [~, csfX(i,:)] = regressLinear(csfX(i,:)', Sd, [], [], perm, RiQ, dR2i);
        end
        for i=1:length(wmIdx)
            [~, wmX(i,:)] = regressLinear(wmX(i,:)', Sd, [], [], perm, RiQ, dR2i);
        end
    end

    % get PCA time-series
    st = 1; ed = st+compNum-1; Xn = [];
    if ~isempty(csfX)
        [coeff,score,~,~,explained,mu] = pca(csfX'); % relation : X' == score * coeff.' + repmat(mu,size(score,1),1);
        Xn(:,st:ed) = score(:,1:compNum);
        st = 1+compNum; ed = st+compNum-1;
    end
    if ~isempty(wmX)
        [coeff,score,~,~,explained,mu] = pca(wmX'); % relation : X' == score * coeff.' + repmat(mu,size(score,1),1);
        Xn(:,st:ed) = score(:,1:compNum);
    end
end
