% get nuisance time-series of ICA-AROMA component time-series based on R.H.R.Pruim et al. (2015).
% returns Xn time-series.
% input:
%  TR           TR of fMRI data
%  csfV         mask for CSF voxels
%  egV          mask for Edge voxels
%  outV         mask for Out voxels
%  M            6 head motion parameters
%  Md           first derivative of 6 head motion parameters
%  icaPath      path of single-subject MELODIC result
%  hfcTh        threshold for High-frequency criterion (default:0.4)
%  csfTh        threshold for CSF criterion (default:0.1)
%  msTh         threshold for Motion-specific criterion (default:[-0.3956, 1.1177])

function [Xn, ICt, idxs] = getNuisanceICAaroma(TR, csfV, egV, outV, M, Md, icaPath, hfcTh, csfTh, msTh)
    if nargin < 10, msTh = [-0.3956, 1.1177]; end
    if nargin < 9, csfTh = 0.1; end
    if nargin < 8, hfcTh = 0.4; end

    Xn = [];
    ICt = [];
    idxs = {};
    % read ICA component for ICA-AROMA
    compf = [icaPath '/melodic_IC.nii.gz'];
    if ~exist(compf,'file')
        disp(['file not found (skipped) : ' compf]);
        return;
    end
    disp(['loading : ' compf]);
    compinfo = niftiinfo(compf);
    IC = single(niftiread(compinfo));
    IC = adjustVolumeDir(IC, compinfo);

    % read ICA component time series for ICA-AROMA
    comptf = [icaPath '/melodic_mix'];
    if ~exist(comptf,'file')
        disp(['file not found (skipped) : ' comptf]);
        return;
    end
    disp(['loading : ' comptf]);
    ICt = readmatrix(comptf);

    compftf = [icaPath '/melodic_FTmix'];
    if ~exist(compftf,'file')
        disp(['file not found (skipped) : ' compftf]);
        return;
    end
    disp(['loading : ' compftf]);
    ICft = readmatrix(compftf);

    if size(ICt,2) ~= size(IC,4)
        disp(['IC component number does not match.' comptf]);
        return;
    end

    cnum = size(IC,4); % component number
    edgeFrac = [];
    csfFrac = [];

    % Edge fraction & CSF fraction.
    for c=1:cnum
        ICC = abs(IC(:,:,:,c));
        totsum = sum(ICC,'all');
        Vcsf = ICC .* csfV;
        Veg = ICC .* egV;
        Vout = ICC .* outV;
        csfsum = sum(Vcsf,'all');
        egsum = sum(Veg,'all');
        outsum = sum(Vout,'all');
        csfFrac(end+1) = csfsum / totsum;
        edgeFrac(end+1) = (outsum + egsum) / (totsum - csfsum);
    end

    % High-frequency content
    Fs = 1 / TR;
    Ny = Fs / 2;
    f = Ny * (0:size(ICft,1)-1) / size(ICft,1);
    ICft(f < 0.01,:) = [];
    ftcum = cumsum(ICft,1);
    ftcum = ftcum ./ repmat(ftcum(end,:),[size(ftcum,1), 1]);
    ftcum = ftcum - 0.5;
    [~, idx] = min(abs(ftcum),[],1); % find close to 0. (0.5 of ftcum)
    hfc = (idx-1)/(size(ICft,1)-1);

    % Maximum correlation with realignment parameters.
    hmp = [M, Md, M.^2, Md.^2];
    hmp = [hmp(2:end-1,:), hmp(1:end-2,:), hmp(3:end,:)]; % 72HMP
    R = corr(ICt(2:end-1,:), hmp);
    mrpCorr = max(abs(R),[],2).';

    % Classify the ICs as motion or non-motion
    idxs{1} = find(hfc > hfcTh);
    idxs{2} = find(csfFrac > csfTh);
    idxs{3} = find(edgeFrac > msTh(1)*mrpCorr+msTh(2));
    idx = unique([idxs{1}, idxs{2}, idxs{3}]); % 'or' conditions
    Xn = ICt(:,idx);
end
