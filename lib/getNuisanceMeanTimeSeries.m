% get nuisance time-series of Global mean, Global signal, mean CSF, mean WM
% returns global mean, global signal, mean csf, mean white matter.
% input:
%  V            4D fMRI time-series
%  csfV         mask for CSF (optional)
%  wmV          mask for White matter (optional)
%  gsV          mask for Global signal (whole brain mask)(optional)

function Xn = getNuisanceMeanTimeSeries(V, csfV, wmV, gsV)
    Xn = nan(size(V,4),4);
    for i=1:size(V,4)
        Vi = V(:,:,:,i);
        Xn(i,1) = nanmean(Vi(:)); % global mean
        if ~isempty(gsV)
            V1 = Vi .* single(gsV); V1(V1<=0) = nan;
            Xn(i,2) = nanmean(V1(:)); % global signal
        end
        if ~isempty(csfV)
            V1 = Vi .* single(csfV); V1(V1<=0) = nan;
            Xn(i,3) = nanmean(V1(:)); % csf mean
        end
        if ~isempty(wmV)
            V1 = Vi .* single(wmV); V1(V1<=0) = nan;
            Xn(i,4) = nanmean(V1(:)); % wm mean
        end
    end
    Xn = Xn - mean(Xn,1);
    if ~isempty(gsV)
        Xn = Xn / (std(Xn(:,2),1)*4); % normalize around [-1 1] range
    end
    if isempty(wmV), Xn(:,4) = []; end
    if isempty(csfV), Xn(:,3) = []; end
    if isempty(gsV), Xn(:,2) = []; end
end
