%%
% Calculate 2nd-Level Mixed-Effects Seed pairwise Granger causality
% analysis (one sample t-test/Wilcoxon version).
% Z2 is main output and T2 will be empty
% input:
%  CY          cells of ROI or voxel time series (time series x node) or filename or matfile object
%  CS          seed cells of ROI or voxel time series (time series x node)(optional)
%  lags        number of lags for autoregression (default:3)
%  isSpot      spot lag time flag (default:false)
%  cachePath   path for cache files (optional)
%  cacheID     cache identifier (optional)
%  method      Z score estimation, exact or approximate (default:'')
%  isNorm      using normalized pairwise Granger causality

function [G2, RSS2, T2, df, Z2] = calcSeedPGCMixed(CY, CS, lags, isSpot, cachePath, cacheID, method, isNorm)
    if nargin < 8, isNorm = true; end
    if nargin < 7, method = ''; end
    if nargin < 6, cacheID = []; end
    if nargin < 5, cachePath = []; end
    if nargin < 4, isSpot = false; end
    if nargin < 3, lags = 3; end
    if nargin < 2, CS = []; end
    if lags==1, isSpot = false; end

    disp('process Mixed-Effects Pairwise Granger causality estimation ...');
    tc = tic;

    % to struct from memory or file name
    CY = toStructCells(CY);
    isempCS = isempty(CS);
    if ~isempCS, CS = toStructCells(CS); else, CS = CY; end

    roinum = size(CY{1}.X,2);
    srcnum = size(CS{1}.X,2);  % seed voxels for Granger Causality
    n = length(CY);
    CG3 = cell(n,1);
%    for j=1:n
    parfor j=1:n % does not work with thread pool
        tt = tic;
        if ~isempty(cachePath)
            gs = loadSeedPGCCache(cachePath,cacheID,lags,isSpot,j);
            if ~isempty(gs)
                CG3{j} = gs;
                t = toc(tt);
                disp(['process session(' num2str(j) ') t=' num2str(t)]);
                continue;
            end
        end

        % 1st-level estimation
        Y = CY{j}.X;
        Z = CS{j}.X;
        RSS1 = nan(roinum,srcnum,'single');
        RSS2 = nan(roinum,srcnum,'single');
        if isSpot, ed=1; else, ed=lags; end
        for k=1:roinum
            y = Y(lags+1:end,k);
            X2 = [];
            for p=1:ed
                X2 = [X2, Y(p:end-lags+(p-1),k)];
            end
            if ~isNorm
                X2 = [X2, ones(size(Y,1)-lags,1)]; % might not be good to add bias
            end
            for i=1:srcnum
                if k==i && isempCS
                    continue;
                end
                X1 = [];
                for p=1:ed
                    X1 = [X1, Z(p:end-lags+(p-1),i)];
                end
                % F = ((RSS1 - RSS2) / (p2 - p1)) / (RSS2 / n - p2)
                if isNorm
                    [~, r1] = regressLinear(y,X1);
                    [~, r2] = regressLinear(y,[X1, X2]);
                else
                    [~, r1] = regressLinear(y,X2);
                    [~, r2] = regressLinear(y,[X2, X1]);
                end
                RSS1(k,i) = r1'*r1;  % p1 = p
                RSS2(k,i) = r2'*r2;  % p2 = 2*p
            end
        end
        Y = []; % clear memory

        % for 1st level (subject) estimation
%        D3 = RSS1 - RSS2;
        D3 = RSS1 ./ RSS2;
        RSS1 = []; RSS2 = []; % clear memory
%        D3(RSS1==0&RSS2==0) = 0; % # / 0 -> inf, % 0 / 0 case should be 0 (not NaN)
        if isNorm
            DN = repmat(nanmean(D3,2),[1 srcnum]);
            G3 = D3 ./ DN;
            D3 = []; DN = []; % clear memory
            G3 = -log(G3); % to center zero (and flip)
    %        G3 = 1 - G3; % to center zero (and flip)
    %        G3(D3==0&DN==0) = 0; % 0 / 0 case should be 0 (not NaN)
        else
            G3 = log(D3);
            D3 = []; % clear memory
        end

        if ~isempty(cachePath)
            CG3{j} = saveSeedPGCCache(cachePath,cacheID,lags,isSpot,j,G3);
        else
            m = struct; m.G3 = G3;
            CG3{j} = m;
        end
        t = toc(tt);
        disp(['process session(' num2str(j) ') t=' num2str(t)]);
    end

    % calc 2nd-level (group) estimation
    % Wilcoxon Rank Sum Test
    % need to care memory consumption
    if nargout >= 5
        sep = ceil(srcnum / 1750); %ceil(n / 12.5); % separate processing to care memory
        step = ceil(srcnum / sep);
        Z2 = zeros(roinum,srcnum,'single');
        for k=1:sep
            if k*step > srcnum
                st = ((k-1)*step+1); ed = srcnum; step = ed-st+1; % last step
            else
                st = ((k-1)*step+1); ed = k*step;
            end
            G3 = zeros(n,roinum,step,'single');
%            for j=1:n
            parfor j=1:n % does not work with thread pool
                tc = tic();
                G3(j,:,:) = CG3{j}.G3(:,st:ed); % this is slow
                t = toc(tc);
                disp([num2str(k) '-' num2str(j) ') load G3(:,' num2str(st) ':' num2str(ed) ') t=' num2str(t)]);
            end
            G3 = reshape(G3, [n roinum*step]);
            Z1 = calcSignrankZ(G3, method);
            Z1 = reshape(Z1, [roinum, step]);
            Z2(:,st:ed) = Z1;
        end
        G2 = []; T2 = []; RSS2 = []; df = n - 1;
    else
        % calc 2nd-level (group) estimation
        % need to care memory consumption
        G2 = zeros(roinum,roinum,'single');
        for j=1:n
            G2 = G2 + CG3{j}.G3; % this is slow
        end
        G2 = G2 / n;
        RSS2 = zeros(roinum,roinum,'single');
        for j=1:n
            R2 = CG3{j}.G3 - G2; % this is slow
            RSS2 = RSS2 + R2 .* R2;
        end
    
        df = n - 1;
        S2 = sqrt(RSS2 / df);
        SE2 = S2 / sqrt(n);
        T2 = G2 ./ SE2;
        % # / 0 -> inf, % 0 / 0 case should be 0 (not NaN)
        T2(G2==0&SE2==0) = 0;
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

function m = loadSeedPGCCache(cachePath,cacheID,lags,isSpot,j)
    if isSpot, sp='sp'; else, sp=''; end
    cachef = [cachePath '/cache-pgc' num2str(lags) sp '-' cacheID '-' num2str(j) '.mat'];
    if exist(cachef,'file')
        m = matfile(cachef);
    else
        m = [];
    end
end

function m = saveSeedPGCCache(cachePath,cacheID,lags,isSpot,j,G3)
    if isSpot, sp='sp'; else, sp=''; end
    cachef = [cachePath '/cache-pgc' num2str(lags) sp '-' cacheID '-' num2str(j) '.mat'];
    save(cachef,'G3','-v7.3');
    m = matfile(cachef);
end
