%SIGNRANK Wilcoxon signed rank test for zero median.
% originally, this code is caming from MATLAB signrank.m
% Because original code takes only vector data input, we expanded this to matrix input.
% returns Z value vector (1 x node)
% input:
%  X          matrix data for one sample Wilcoxon signed rank test (sample n x node)
%  method     Z score estimation, exact or approximate (default:'')

function Z = calcSignrankZ(X, method)
    if nargin < 2, method = ''; end
    n = size(X,1);
    [tie_rank, tieadj] = tiedrankMat(abs(X));
    TR = tie_rank;
    TR(X<=0) = 0;
    W = sum(TR, 1);
    clear TR; clear X;
    if strcmp(method, 'exact')
        P = statsrexactMat(tie_rank, W);
        P(P>1) = 1;
		P = 2 * P;
        Z = norminv(1-P/2);
        % Z == 0 is twice high. replace it by approximated value
        idx = find(Z==0);
        % get sign
        Zs = (W-n*(n+1)/4) ./ sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
        Zt = Zs;
        Zs(Zs>=0) = 1;
        Zs(Zs<0) = -1;
        Z = Z .* Zs;
        Z(idx) = Zt(idx);
    else
        % this is 'approximate' method of signrank.m
        Z = (W-n*(n+1)/4) ./ sqrt((n*(n+1)*(2*n+1) - tieadj)/24);
    end
end

% --------------------------------
%STATSREXACT Compute exact tail probability for signed rank statistic (matrix version).
function P = statsrexactMat(V, W)
    n = size(V,1);
    V = sort(V,1);

    maxw = n*(n+1)/2;
    folded = (W > maxw/2);
    W = W .* (1-folded) + (maxw-W) .* folded;

    % do not calculate full voxels, find unique pattern of inputs
    tc = tic;
    U = unique(W); U(isnan(U))=[];
    V2 = V - (1:size(V,1))';
    V2 = V2 .* (1:size(V,1))';
    W2 = sum(abs(V2),1);
    clear V2;

    CV = {};
    CW = {};
    CIdx = {};
    for i=1:length(U)
        idx = find(W==U(i)); % unique pattern of W
        W2i = W2(idx);
        U2 = unique(W2i); U2(isnan(U2))=[];
        for j=1:length(U2)
            idx2 = find(W2i==U2(j)); % unique pattern of V within W==U(i)

            % check first one
            CV{end+1} = V(:,idx(idx2(1)))'; % make sure it's a row;
            CW{end+1} = W(idx(idx2(1)));
            CIdx{end+1} = idx(idx2(:));
        end
    end
    clear V; clear W;
    CP = cell(1,length(CV));
    parfor i=1:length(CV)
        % check first one
        v = CV{i}; % make sure it's a row
        w = CW{i};

        % multiply by 2 to force everything to integer.
        if any(v~=floor(v))
            v = round(2*v);
            w = round(2*w);
        end
        C = zeros(w+1,1);
        C(1) = 1; top = 1;
        for vj=v(v<=w)
           newtop = min(top+vj,w+1);
           hi = min(vj,w+1)+1:newtop;
           lo = 1:length(hi);
        
           C(hi) = C(hi) + C(lo);
           top = newtop;
        end
        C = C / (2^n);
        CP{i} = sum(C); % Get tail probability
    end
    P = nan(1,length(W2));
    for i=1:length(CP)
        P(CIdx{i}) = CP{i};
    end
%{
    parfor i=1:length(P)
        v = V(:,i)'; % make sure it's a row
        w = W(i);
        if isnan(w), continue; end

        % multiply by 2 to force everything to integer.
        if any(v~=floor(v))
            v = round(2*v);
            w = round(2*w);
        end
        C = zeros(w+1,1);
        C(1) = 1; top = 1;
        for vj=v(v<=w)
           newtop = min(top+vj,w+1);
           hi = min(vj,w+1)+1:newtop;
           lo = 1:length(hi);
        
           C(hi) = C(hi) + C(lo);
           top = newtop;
        end
        C = C / (2^n);
        P(i) = sum(C); % Get tail probability
    end
%}
    t = toc(tc);
    disp(['statsrexactMat t=' num2str(t)]);
end

% --------------------------------
function [SX, tieadj] = tiedrankMat(X)
    %TR Local tiedrank function to compute results for one column
    % Sort, then leave the NaNs (which are sorted to the end) alone
    tc = tic;
    sz1 = size(X,1);
    [SX, rowidx] = sort(X, 1);
    if isa(X,'single') 
        SX = single(SX);
        rowidx = single(rowidx);
    end
    numNaNs = sum(isnan(X), 1);
    xLen = sz1 - numNaNs;

    CRnk = cell(1, size(X,2));
    tieadj = zeros(1, size(X,2), 'single');

    % Adjust for ties.  Avoid using diff(sx) here in case there are infs.
    ties = SX(1:sz1-1,:) >= SX(2:sz1,:);

%    for j = 1:size(X,2)
    parfor j = 1:size(X,2)
        idx = find(ties(:,j));
        if isempty(idx), continue; end

        tieloc = [idx; xLen(j)+2];
        maxTies = numel(tieloc);
    
        % Use ranks counting from low end
        ranks = [1:xLen(j) NaN(1,numNaNs(j))]';
        tiecount = 1;
        while (tiecount < maxTies)
            tiestart = tieloc(tiecount);
            ntied = 2;
            while(tieloc(tiecount+1) == tieloc(tiecount)+1)
                tiecount = tiecount+1;
                ntied = ntied+1;
            end
            tieadj(j) = tieadj(j) + ntied*(ntied-1)*(ntied+1)/2;
            
            % Compute mean of tied ranks
            ranks(tiestart:tiestart+ntied-1) = ...
                          sum(ranks(tiestart:tiestart+ntied-1)) / ntied;
            tiecount = tiecount + 1;
        end
        CRnk{j} = ranks;
    end
    for j = 1:size(X,2)
        % Broadcast the ranks back out, including NaN where required.
        if isempty(CRnk{j})
            SX(rowidx(:,j), j) = [1:xLen(j) NaN(1,numNaNs(j))];
        else
            SX(rowidx(:,j), j) = CRnk{j};
        end
    end
    t = toc(tc);
    disp(['tiedrankMat t=' num2str(t)]);
end
