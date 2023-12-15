%%
% Calculate 2nd-Level Mixed-Effects Seed correlation analysis
% if CS is not suppplied, CY correlation matrix is calculated.
% (one sample t-test/Wilcoxon version)
% input:
%  CY          cells of ROI or voxel time series (time series x node) or filename or matfile object
%  CS          seed cells of ROI or voxel time series (time series x node)(optional)
%  cachePath   path for cache files (optional)
%  cacheID     cache identifier (optional)

function [B2, RSS2, T2, df, R, Z2] = calcSeedCorrMixed(CY, CS, cachePath, cacheID)
    if nargin < 4, cacheID = []; end
    if nargin < 3, cachePath = []; end
    if nargin < 2, CS = []; end
    isOutR = (nargout >= 5);

    disp('process Mixed-Effects correlation ...');
    tc = tic;

    % to struct from memory or file name
    CY = toStructCells(CY);
    if ~isempty(CS), CS = toStructCells(CS); end

    idx = [];
    if isempty(CS)
        v = size(CY{1}.X,2);
        T = triu(ones(v,v,'single'));
        idx = uint32(find(T==1));
        clear T;
    end

    n = length(CY); % subject length
    df = n - 1;
    CR = cell(n,1);
%    for j=1:n
    parfor j=1:n
        tt = tic;
        if ~isempty(cachePath)
            gs = loadSeedCorrCache(cachePath,cacheID,j);
            if ~isempty(gs)
                CR{j} = gs;
                t = toc(tt);
                disp(['process session(' num2str(j) ') t=' num2str(t)]);
                continue;
            end
        end

        % 1st-level estimation
        if isempty(CS)
            R2 = calcSeedCorr(CY{j}.X);
        else
            R2 = calcSeedCorr(CY{j}.X,CS{j}.X);
        end
        R2(isnan(R2)) = 0; % node time-series might be zero

        % for 2nd-level Y vector
        if isempty(CS)
            R2 = R2(idx); % to shrink memory
        else
            R2 = R2(:);
        end

        if ~isempty(cachePath)
            CR{j} = saveSeedCorrCache(cachePath,cacheID,j,R2);
        else
            m = struct; m.R2 = R2; % on memory
            CR{j} = m;
        end
        t = toc(tt);
        disp(['process session(' num2str(j) ') t=' num2str(t)]);
    end

    % calc 2nd-level (group) estimation
    % Wilcoxon Rank Sum Test
    if nargout >= 6
        R3 = zeros(n,length(CR{1}.R2),'single');
        for j=1:n
            R3(j,:) = CR{j}.R2;
        end
        Z2 = calcSignrankZ(R3);
        if isempty(CS)
            Z3 = zeros(v,v,'single');
            Z3(idx) = Z2;
            Z2 = Z3' + triu(Z3,1);
        else
            Z2 = reshape(Z2,size(CY{1}.X,2),size(CS{1}.X,2));
        end
    end

    % calc 2nd-level estimation
    % need to care memory consumption
    B2 = zeros(size(CR{1}.R2,1),1,'single');
    for j=1:n
        B2 = B2 + CR{j}.R2;
    end
    B2 = B2 / n; % B2 = mean(R1,2);
    RSS2 = zeros(size(CR{1}.R2,1),1,'single');
    if isOutR, R = nan(size(CR{1}.R2,1),n,'single'); end
    for j=1:n
        r = CR{j}.R2 - B2;
        RSS2 = RSS2 + r .* r;
        if isOutR, R(:,j) = r; end
    end

    S2 = sqrt(RSS2 / df);
    SE2 = S2 / sqrt(n);
    T2 = B2 ./ SE2;
    % # / 0 -> inf, % 0 / 0 case should be 0 (not NaN)
    T2(B2==0&SE2==0) = 0;

    if isempty(CS)
        T3 = zeros(v,v,'single');
        T3(idx) = T2;
        T2 = T3' + triu(T3,1);
    else
        T2 = reshape(T2,size(CY{1}.X,2),size(CS{1}.X,2));
    end

    t = toc(tc);
    disp(['done t=' num2str(t) 'sec'])
end


function CY = toStructCells(CY)
    for j=1:length(CY)
        if ischar(CY{j})
            CY{j} = matfile(CY{j});
        elseif isnumeric(CY{j})
            st = struct;
            st.X = CY{j};
            CY{j} = st;
        else
            % if it is matfile or struct, ignore it.
        end
    end
end

function m = loadSeedCorrCache(cachePath,cacheID,j)
    cachef = [cachePath '/cache-corr-' cacheID '-' num2str(j) '.mat'];
    if exist(cachef,'file')
        m = matfile(cachef);
    else
        m = [];
    end
end

function m = saveSeedCorrCache(cachePath,cacheID,j,R2)
    cachef = [cachePath '/cache-corr-' cacheID '-' num2str(j) '.mat'];
    save(cachef,'R2','-v7.3');
    m = matfile(cachef);
end
