%%
% Plot multi seed analysis results
% returns
%   cells of thresholded T or Z value 3D voxels (seeds are source, plus side) (CVsp), cells of thresholded T or Z value 3D voxels (seeds are source, minus side) (CVsm)
%   cells of thresholded T or Z value 3D voxels (seeds are target, plus side) (CVtp), cells of thresholded T or Z value 3D voxels (seeds are target, minus side) (CVtm)
% input:
%  tgNames      cells of target region names
%  Tgs          cells of target region vectors (node x 1)
%  T2           T-value or Z-value matrix (node x node)
%  thParam      threshold params for T-value matrix {df, Pth, corrMeth}
%            df:       degree of freedom
%            Pth:      P-value threshold for T or Z value matrix
%            corrMeth: family-wise error rate (FWER) correction method for P-value ([default] 'none','bonf','sidak','holm-bonf','holm-sidak')
%  atlasV       3D atlas volume (or {0, 1} mask volume)
%  isMask       if atlasV is mask volume, then true. otherwise, atlasV has regional info, and also used as mask.
%  isRtoL       X axis is right to left (default: true)
%  backV        3D background volume for plotNifti3DAxes (default: [])
%  sessionName  session name for title (optional)
%  rangePlus    plus T-value range ([min max]) (optional)
%  rangeMinus   minus T-value range ([min max]) (optional)
%  cmap         color map for 3D plot (default: hot)
%  proj         projection type of T-values or Z-values of one voxel
%  zidx         show z-index slices (default: [])
%  flatXY       show functional flat map (default: [])

function [CVsp, CVsm, CVtp, CVtm, Tmaxs, Tcnts] = plotSeedCorrImage(tgNames, Tgs, T2, thParam, atlasV, isMask, isRtoL, backV, sessionName, rangePlus, rangeMinus, cmap, proj, zidx, flatXY)
    if nargin < 15, flatXY = []; end
    if nargin < 14, zidx = []; end
    if nargin < 13, proj = 'max'; end
    if nargin < 12, cmap = hot; end
    if nargin < 11, rangeMinus = []; end
    if nargin < 10, rangePlus = []; end
    if nargin < 9, sessionName = ''; end
    if nargin < 8, backV = true; end
    if nargin < 7, isRtoL = true; end

    % init threshold parameters
    df = thParam{1}; Pth = thParam{2}; corrMeth = 'none';
    if length(thParam)>=3, corrMeth = thParam{3}; end
    if isinf(df), Tstr='Z'; else Tstr='T'; end % T-value or Z-value
    if isempty(proj), proj='max'; end

    if isempty(rangePlus), rangePlus = [nan nan]; end
    if isempty(rangeMinus), rangeMinus = [nan nan]; end
    orgRp = rangePlus; orgRm = rangeMinus;

    Tmaxs = [];
    Tcnts = [];
    CVsp = {}; CVsm = {};
    CVtp = {}; CVtm = {};

    for j=1:length(tgNames)
        tg = Tgs{j};
        tg(tg==0) = nan;

        % 'target' is source contrast 
        T2c = T2 .* tg'; % size should match
        if strcmp(proj, 'median')
            T2p = median(T2c,2,'omitnan');
            T2m = median(-T2c,2,'omitnan');
        elseif strcmp(proj, 'min')
            T2p = min(T2c,[],2,'omitnan');
            T2m = min(-T2c,[],2,'omitnan');
        else
            T2p = max(T2c,[],2,'omitnan');
            T2m = max(-T2c,[],2,'omitnan');
        end

        % T-value threshold (Bonferroni correction / Šidák correction)
        if strcmp(corrMeth,'bonf') || strcmp(corrMeth,'sidak')
            % only masked ROI is used
            if size(T2,1) == size(T2,2) % square size
                m = size(T2,1) * sum(~isnan(tg));
            else
                m = size(T2,1); % basically, column ROIs are independent.
            end
            if strcmp(corrMeth,'bonf')
                BPth = Pth / m;
            else
                % Šidák correction - almost same as bonferroni
                BPth = 1 - power(1 - Pth, 1/m);
            end
            if isinf(df)
                Tth = abs(norminv(BPth));
            else
                Tth = abs(tinv(BPth,df));
            end
            disp(['P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), ' Tstr '-value=' num2str(Tth)])
        % Holm–Bonferroni method
        elseif strcmp(corrMeth,'holm-bonf') || strcmp(corrMeth,'holm-sidak')
            T2c2 = T2 .* tg'; % size should match
            T2c2(isnan(T2c2)) = []; % remove NaN
            T2s = sort(T2c2(:),'descend'); % full matrix T-test
            for k=1:length(T2s)
                if strcmp(corrMeth,'holm-bonf')
                    BPth = Pth / (length(T2s) + 1 - k);
                else
                    BPth = 1 - power(1 - Pth, 1/(length(T2s) + 1 - k));
                end
                if isinf(df)
                    Tth = abs(norminv(BPth));
                else
                    Tth = abs(tinv(BPth,df));
                end
                if T2s(k) <= Tth
                    break;
                end
            end
            disp(['P-value=' num2str(Pth) ' (' corrMeth ' corrected=' num2str(BPth) '), ' Tstr '-value=' num2str(Tth)])
        else
            if isinf(df)
                Tth = abs(norminv(Pth));
            else
                Tth = abs(tinv(Pth,df));
            end
            disp(['P-value=' num2str(Pth) ', ' Tstr '-value=' num2str(Tth)])
        end

        % set auto range
        if isnan(orgRp(1))
            rangePlus(1) = Tth;
        end
        if isnan(orgRp(2))
            m = max(T2p(~isinf(T2p)));
            if m>Tth, rangePlus(2) = m; else, rangePlus(2) = Tth + 1; end
        end
        if isnan(orgRm(1))
            rangeMinus(1) = Tth;
        end
        if isnan(orgRm(2))
            m = max(T2m(~isinf(T2m)));
            if m>Tth, rangeMinus(2) = m; else, rangeMinus(2) = Tth + 1; end
        end

        T2p(T2p<Tth) = nan;   % thresholded T-value
        T2m(T2m<Tth) = nan;   % thresholded T-value
        if isMask == 1
            aIdx = find(atlasV(:) > 0);
            Vp = single(atlasV);
            Vp(:) = nan;
            Vm = Vp;
            Vp(aIdx) = T2p;
            Vm(aIdx) = T2m;
        else
            Vp = getNifti4DFromRoiTS(T2p, atlasV);
            Vm = getNifti4DFromRoiTS(T2m, atlasV);
        end

        % show figure
        figure; plotNifti3DAxes(Vp,'max',rangePlus,cmap,backV,gray,isRtoL);
        sgtitle(['From Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');
        figure; plotNifti3DAxes(Vm,'max',rangeMinus,cmap,backV,gray,isRtoL);
        sgtitle(['From Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');

        if ~isempty(zidx)
            figure; plotNifti3DZSlice(Vp,zidx,rangePlus,cmap,backV,gray,isRtoL,0.15);
            sgtitle(['From Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');

            figure; plotNifti3DZSlice(Vm,zidx,rangeMinus,cmap,backV,gray,isRtoL,0.15);
            sgtitle(['From Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');
        end

        if ~isempty(flatXY)
            figure; plotNifti3Dflatmap(T2p, atlasV, isMask, flatXY, 10, rangePlus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['From Seeds (plus) of ' sessionName ' : ' tgNames{j}]);

            figure; plotNifti3Dflatmap(T2m, atlasV, isMask, flatXY, 10, rangeMinus, cmap, [0.1 0.1 0.1], [0 0 0]);
            title(['From Seeds (minus) of ' sessionName ' : ' tgNames{j}]);
        end
        
        % show Tmax
        tmax = max(Vp(:));
        tcnt = length(find(Vp>=Tth));
        disp([Tstr 'max (source) of ' sessionName ' : ' tgNames{j} ' tmax=' num2str(tmax) ', tcnt=' num2str(tcnt)])

        Vp(Vp > rangePlus(2)) = rangePlus(2);
        Vm(Vm > rangeMinus(2)) = rangeMinus(2);
        CVsp{j} = Vp;
        CVsm{j} = Vm;

        % 'target' is target voxels
        if size(T2,1) == size(tg,1) && ~isequal(T2,T2')
            T2c = T2 .* tg; % size should match
            if strcmp(proj, 'median')
                T2p = median(T2c,1,'omitnan');
                T2m = median(-T2c,1,'omitnan');
            elseif strcmp(proj, 'min')
                T2p = min(T2c,[],1,'omitnan');
                T2m = min(-T2c,[],1,'omitnan');
            else
                T2p = max(T2c,[],1,'omitnan');
                T2m = max(-T2c,[],1,'omitnan');
            end
            T2p = T2p(1:end-1)'; % remove intercept ??
            T2m = T2m(1:end-1)'; % remove intercept ??

            T2p(T2p<Tth) = nan;   % thresholded T-value
            T2m(T2m<Tth) = nan;   % thresholded T-value
            if isMask == 1
                aIdx = find(atlasV(:) > 0);
                Vp = single(atlasV);
                Vp(:) = nan;
                Vm = Vp;
                Vp(aIdx) = T2p;
                Vm(aIdx) = T2m;
            else
                Vp = getNifti4DFromRoiTS(T2p, atlasV);
                Vm = getNifti4DFromRoiTS(T2m, atlasV);
            end
    
            % show figure
            figure; plotNifti3DAxes(Vp,'max',rangePlus,cmap,backV,gray,isRtoL);
            sgtitle(['To Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');
            figure; plotNifti3DAxes(Vm,'max',rangeMinus,cmap,backV,gray,isRtoL);
            sgtitle(['To Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');
    
            if ~isempty(zidx)
                figure; plotNifti3DZSlice(Vp,zidx,rangePlus,cmap,backV,gray,isRtoL,0.15);
                sgtitle(['To Seeds (plus) of ' sessionName ' : ' tgNames{j}],'Color','white');
    
                figure; plotNifti3DZSlice(Vm,zidx,rangeMinus,cmap,backV,gray,isRtoL,0.15);
                sgtitle(['To Seeds (minus) of ' sessionName ' : ' tgNames{j}],'Color','white');
            end

            if ~isempty(flatXY)
                figure; plotNifti3Dflatmap(T2p, atlasV, isMask, flatXY, 10, rangePlus, cmap, [0.1 0.1 0.1], [0 0 0]);
                title(['To Seeds (plus) of ' sessionName ' : ' tgNames{j}]);
    
                figure; plotNifti3Dflatmap(T2m, atlasV, isMask, flatXY, 10, rangeMinus, cmap, [0.1 0.1 0.1], [0 0 0]);
                title(['To Seeds (minus) of ' sessionName ' : ' tgNames{j}]);
            end

            % show Tmax
            tmax = max(Vp(:));
            tcnt = length(find(Vp>=Tth));
            disp([Tstr 'max (target) of ' sessionName ' : ' tgNames{j} ' tmax=' num2str(tmax) ', tcnt=' num2str(tcnt)])

            Vp(Vp > rangePlus(2)) = rangePlus(2);
            Vm(Vm > rangeMinus(2)) = rangeMinus(2);
            CVtp{j} = Vp;
            CVtm{j} = Vm;
        end

        % output 
        Tmaxs(j) = tmax;
        Tcnts(j) = tcnt;
    end
end
