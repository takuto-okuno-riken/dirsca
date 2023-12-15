%%
% Calculate Seed analysis of correlation (functional connectivity)
% This function is used for single session or multi-session (fixed-effects).
% returns correlation coefficient matrix (R), p-value matrix (RSS)
% input:
%  Y         ROI or voxel time series (time series x node)
%  S         ROI or voxel time series (time series x node)(optional)

function [R, P] = calcSeedCorr(Y, S)
    if nargin < 2, S = []; end
    disp('process seed analysis of correlation ...');
    if isempty(S)
        [R, P] = corr(Y);
    else
        [R, P] = corr(Y, S);
    end
end
