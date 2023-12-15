% get nuisance time-series of tCompCor (high tSTD) component
% time-series based on Y.Behzadi et al. (2007).
% returns high tSTD <compNum> components.
% input:
%  V            4D fMRI time-series
%  Sd           detrend signal (time-series x node) before PCA (optional)
%  vTh          threshold (voxel number) for high tSTD of tCompCor (default:1000)
%  compNum      component number of tCompCor (default:6)

function Xn = getNuisancetCompCor(V, Sd, vTh, compNum)
    if nargin < 6, compNum = 6; end
    if nargin < 3, vTh = 1000; end
    if nargin < 2, Sd = []; end
    
    % get ROI index for aCompCor
    V = reshape(V, size(V,1)*size(V,2)*size(V,3),size(V,4));
    M = nanmean(V,2);
    V = V - M;
    S = nanstd(V,1,2);
    [~,I] = sort(S,'descend');
    idx = I(1:vTh);

    % detrend
    hsX = V(idx,:);
    if ~isempty(Sd)
        [~, ~, perm, RiQ, dR2i] = regressPrepare(Sd);
        % 1st step OLS regression
        for i=1:length(idx)
            [~, hsX(i,:)] = regressLinear(hsX(i,:)', Sd, [], [], perm, RiQ, dR2i);
        end
    end

    % get PCA time-series
    Xn = nan(size(V,2),compNum);
    [coeff,score,~,~,explained,mu] = pca(hsX'); % relation : X' == score * coeff.' + repmat(mu,size(score,1),1);
    Xn(:,1:compNum) = score(:,1:compNum);
end
