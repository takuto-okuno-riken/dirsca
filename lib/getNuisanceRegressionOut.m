function Vout = getNuisanceRegressionOut(V, Xn, maskV)
    sz = size(V,1)*size(V,2)*size(V,3);
    V = reshape(V,sz,size(V,4));
    M = nanmean(V,2);
    V = V - M;
    [~, ~, perm, RiQ, dR2i] = regressPrepare(Xn);
    CR = cell(1, sz);
    idxs = find(maskV > 0);
    for i=1:length(idxs)
%    parfor i=1:length(idxs)
        idx = idxs(i);
        % 1st step OLS regression
        [~, CR{i}] = regressLinear(V(idx,:)', Xn, [], [], perm, RiQ, dR2i);
    end
    for i=1:length(idxs)
        idx = idxs(i);
        V(idx,:) = CR{i};
    end
    V = V + M;
    Vout = reshape(V,size(maskV,1),size(maskV,2),size(maskV,3),size(V,2));
end
