%%
% directional and non-direcional SCA command line tool

function dirsca(varargin)

    % set version number
    versionNumber = '0.1';

    % add script path
    if ~isdeployed % checking MATLAB mode or stand-alone mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind).file;
        [exedir,exename,ext] = fileparts(relpath);
        if exist([exedir '/util'],'dir')
            addpath([exedir '/util']);
            addpath([exedir '/lib']);
        end
    end

    % get exe file full path
    global exePath;
    global exeName;
    [exePath, exeName, ext] = exeFilename();

    % init command line input
    handles.commandError = 0;
    handles.full = 0;
    handles.dir = 0;
    handles.nondir = 0;
    handles.niiFiles = {};
    handles.atlasFile = 'data/human_2mm_cubeRoi2.nii.gz';
    handles.seedFile = [];
    handles.compSeed = 0;
    handles.maskFile = [];
    handles.listFile = [];
    handles.rmframe = 10;
    handles.nuisance = 'gmacomp';
    handles.compNum = 6;
    handles.nuibr = 'data/human_2mm_brain_mask.nii.gz';
    handles.nuiwm = 'data/human_2mm_wm.nii.gz';
    handles.nuicsf = 'data/human_2mm_csf.nii.gz';
    handles.highpass = 0;
    handles.smooth = [];
    handles.lags = 0;
    handles.rankmeth = 'exact';
    handles.poolnum = 0;
    handles.outpath = 'results';
    handles.cachepath = 'results/cache';

    handles.showSig = 0;
    handles.showRas = 0;
    handles.noCache = 0;

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-d'}
                handles.dir = 1;
            case {'-n'}
                handles.nondir = 1;
            case {'-f','--full'}
                handles.full = 1;
            case {'-s','--seed'}
                handles.seedFile = varargin{i+1};
                i = i + 1;
            case {'--compseed'}
                handles.compSeed = 1;
            case {'--atlas'}
                handles.atlasFile = varargin{i+1};
                i = i + 1;
            case {'--mask'}
                handles.maskFile = varargin{i+1};
                i = i + 1;
            case {'--list'}
                handles.listFile = varargin{i+1};
                i = i + 1;
            case {'--nui'}
                handles.nuisance = varargin{i+1};
                i = i + 1;
            case {'--compNum'}
                handles.compNum = str2num(varargin{i+1});
                i = i + 1;
            case {'--nuibr'}
                handles.nuibr = varargin{i+1};
                i = i + 1;
            case {'--nuiwm'}
                handles.nuiwm = varargin{i+1};
                i = i + 1;
            case {'--nuicsf'}
                handles.nuicsf = varargin{i+1};
                i = i + 1;
            case {'--highpass'}
                handles.highpass = str2num(varargin{i+1});
                i = i + 1;
            case {'--rmframe'}
                handles.rmframe = str2num(varargin{i+1});
                i = i + 1;
            case {'--smooth'}
                handles.smooth = [str2num(varargin{i+1}) str2num(varargin{i+2}) str2num(varargin{i+3})];
                i = i + 3;
            case {'--lags'}
                handles.lags = str2num(varargin{i+1});
                i = i + 1;
            case {'--rankmeth'}
                handles.rankmeth = varargin{i+1};
                i = i + 1;
            case {'--pool'}
                handles.poolnum = str2num(varargin{i+1});
                i = i + 1;
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--cachepath'}
                handles.cachepath = varargin{i+1};
                i = i + 1;
            case {'--showsig'}
                handles.showSig = 1;
            case {'--showras'}
                handles.showRas = 1;
            case {'--nocache'}
                handles.noCache = 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'-v','--version'}
                disp([exeName ' version : ' num2str(versionNumber)]);
                return;
            otherwise
                if strcmp(varargin{i}(1), '-')
                    disp(['bad option : ' varargin{i}]);
                    i = size(varargin, 2);
                    handles.commandError = 1;
                else
                    handles.niiFiles = [handles.niiFiles varargin{i}];
                end
        end
        i = i + 1;
    end
    
    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif handles.dir == 0 && handles.nondir == 0
        disp('Please specify directional and/or non-directional analysis')
        showUsage();
        return;
    elseif (isempty(handles.niiFiles) && isempty(handles.listFile)) || (handles.full==0 && isempty(handles.seedFile)) || (isempty(handles.atlasFile) && isempty(handles.maskFile))
        disp('Seed NIfTI and input NIfTI files are required. please specify nifti files.');
        showUsage();
        return;
    end

    % process input files
    processInputFiles(handles);
end

%%
% show usage function
function showUsage()
    global exePath;
    global exeName;
    
    disp(['usage: ' exeName ' [options][-d][-n][-f][-s seed.nii.gz] file1.nii.gz ...']);
    disp('  -d                  perform directional seed-based connectivity analysis');
    disp('  -n                  perform non-directional seed-based connectivity analysis');
    disp('  -f, --full          full voxel connectivity analysis');
    disp('  -s, --seed file     NIfTI <file> of seed ROI voxels');
    disp('  --compseed          seed is 4D component file');
    disp('  --atlas file        NIfTI <file> of target ROI voxels (default: data/human_2mm_cubeRoi2.nii.gz)');
    disp('  --mask file         NIfTI <file> of target mask (full voxel) for seed analysis (Preferred over atlas option)');
    disp('  --list file         text <file> of subject NIfTI file list');
    disp('  --rmframe num       remove first <num> frames (default: 10)');
    disp('  --smooth x y z      gaussian kernel smoothing FWHM=[x y z] voxels (default: off)');
    disp('  --nui algo          nuisance factor removal method (default: "gmacomp")');
    disp('  --compnum num       component number for aCompCor (default: 6)');
    disp('  --nuibr file        brain mask NIfTI <file> for global signal of nuisance factor removal (default: data/human_2mm_brain_mask.nii.gz)');
    disp('  --nuiwm file        white matter mask NIfTI <file> for wm-mean of nuisance factor removal (default: data/human_2mm_wm.nii.gz)');
    disp('  --nuicsf file       csf mask NIfTI <file> for csf-mean of nuisance factor removal (default: data/human_2mm_csf.nii.gz)');
    disp('  --highpass freq     high-pass NIfTI <freq> Hz (default: off)');
    disp('  --lags num          spot time lag for directional SCA (default: auto)');
    disp('  --rankmeth type     ranking method for directional SCA (default: "exact")');
    disp('  --outpath path      output files <path> (default:"results")');
    disp('  --cachepath path    cache path <path> (default:"results/cache")');
    disp('  --showsig           show processed time-series of input NIfTI file');
    disp('  --showras           show raster plot of processed time-series of input NIfTI file');
    disp('  --nocache           do not use cache file for conversion');
    disp('  --pool num          working pool number for parallel calculation');
    disp('  -v, --version       show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;

    % load list file
    if ~isempty(handles.listFile)
        fid = fopen(handles.listFile);
        tline = fgetl(fid);
        while ischar(tline)
            handles.niiFiles = [handles.niiFiles tline];
            tline = fgetl(fid);
        end
        fclose(fid);
    end
    N = length(handles.niiFiles);

    % load SCA target voxels (cube ROI or mask)
    if ~isempty(handles.maskFile) && ~exist(handles.maskFile,'file')
        disp(['mask file is not found : ' handles.maskFile]);
        return;
    end
    if ~isempty(handles.atlasFile) &&  ~exist(handles.atlasFile,'file')
        disp(['atlas file is not found : ' handles.atlasFile]);
        return;
    end
    ismask = ~isempty(handles.maskFile);
    if ismask
        atlasFile = handles.maskFile;
    else
        atlasFile = handles.atlasFile;
    end
    [~,name,~] = fileparts(atlasFile);
    atlasName = strrep(name,'.nii','');

    disp(['load SCA target voxels : ' atlasFile]);
    atlasInfo = niftiinfo(atlasFile);
    atlasV = niftiread(atlasInfo);
    atlasV = adjustVolumeDir(atlasV,atlasInfo);
    if ismask
        maskIdx = find(atlasV > 0);
        nodeNum = length(maskIdx);
    else
        R = unique(atlasV);
        R(R==0) = []; % remove 0
        nodeNum = length(R);
    end

    % load SCA seed voxels
    if handles.full == 0
        if ~isempty(handles.seedFile) && ~exist(handles.seedFile,'file')
            disp(['seed file is not found : ' handles.seedFile]);
            return;
        end
        disp(['load SCA seed voxels : ' handles.seedFile]);
        seedInfo = niftiinfo(handles.seedFile);
        seedV = niftiread(seedInfo);
        seedV = adjustVolumeDir(seedV,seedInfo);
        if size(seedV,4) > 1 % 4D type 
            seedNum = size(seedV,4);
        else
            R = unique(seedV);
            R(R==0) = []; % remove 0
            seedNum = length(R);
        end
        [~,name,~] = fileparts(handles.seedFile);
        seedName = strrep(name,'.nii','');
    else
        seedNum = nodeNum;
        seedName = 'full';
    end

    % load mask files for Nuisance Signal Regression
    csfV = []; wmV = []; brV = [];
    if ~isempty(handles.nuicsf)
        disp(['load CSF mask : ' handles.nuicsf]);
        csfinfo = niftiinfo(handles.nuicsf);
        csfV = single(niftiread(csfinfo));
        csfV = adjustVolumeDir(csfV,csfinfo);
        if max(csfV(:)) > 1
            csfV = single(csfV) / 255; % to [0 1] range
        end
    end
    if ~isempty(handles.nuiwm)
        disp(['load white matter mask : ' handles.nuiwm]);
        wminfo = niftiinfo(handles.nuiwm);
        wmV = single(niftiread(wminfo));
        wmV = adjustVolumeDir(wmV,wminfo);
        if max(wmV(:)) > 1
            wmV = single(wmV) / 255; % to [0 1] range
        end
    end
    if ~isempty(handles.nuibr)
        disp(['load whole brain mask : ' handles.nuibr]);
        brinfo = niftiinfo(handles.nuibr);
        brV = single(niftiread(brinfo));
        brV = adjustVolumeDir(brV,brinfo);
        brV(brV>=1) = 1;
        brV(brV<1) = 0;
    end

    % check and make cache directory
    if ~exist(handles.cachepath,'dir')
        mkdir(handles.cachepath);
    end

    % set smoothing sigma
    sigma = [];
    if ~isempty(handles.smooth)
        sigma = handles.smooth / sqrt(8*log(2));  % voxel size to sigma
        filterSize = 2*ceil(2*sigma)+1;
    end

    if ~isempty(handles.smooth), smstr=['-s' num2str(handles.smooth(1)) '-' num2str(handles.smooth(2)) '-' num2str(handles.smooth(3))]; else, smstr=''; end
    if ~isempty(handles.nuisance), nstr=['-' handles.nuisance]; else, nstr=''; end

    % process each file (extract ROI time-series)
    CX = {}; names = {}; TR = 0;
    if handles.full == 0, CS = {}; else, CS = []; end
    for i = 1:N
        % load NIfTI files
        argv = handles.niiFiles{i};
        flist = dir(argv);
        if isempty(flist)
            disp(['file is not found. ignoring : ' argv]);
            continue;
        end
        for k=1:length(flist)
            % init data
            X = [];

            fname = [flist(k).folder filesep flist(k).name];
            [path,name,ext] = fileparts(fname);
            name = strrep(name,'.nii','');
            pathes = split(flist(k).folder, filesep);
            cpath = '';
            if ~isempty(pathes)
                cpath = [pathes{end} '-'];
            end

            % read nifti file
            cachename = [handles.cachepath '/roisig-' atlasName '-r' num2str(handles.rmframe) smstr nstr '-' cpath name];
            if handles.noCache > 0 || ~exist([cachename '.mat'],'file') || (handles.full==0 && ~exist([cachename '-' seedName '.mat'],'file'))
                disp(['processing : ' name]);
                info = niftiinfo(fname);
                TR = info.PixelDimensions(4);
                V = single(niftiread(info));
                V = adjustVolumeDir(V,info);
                V(isnan(V)) = 0;

                if size(V,1) ~= size(atlasV,1) || size(V,2) ~= size(atlasV,2) || size(V,3) ~= size(atlasV,3)
                    disp(['atlas space and fMRI space does not match. ignoring : ' fname]);
                    continue;
                end

                % remove first ## frames
                if handles.rmframe > 0 
                    disp(['remove first ' num2str(handles.rmframe) ' frames.']);
                    V = V(:,:,:,handles.rmframe:end);
                end

                % apply gaussian kernel smoothing
                if ~isempty(sigma)
                    disp(['apply gaussian kernel smoothing [' num2str(handles.smooth(1)) ' ' num2str(handles.smooth(2)) ' ' num2str(handles.smooth(3)) '] voxels.']);
                    for n=1:size(V,4)
                        V(:,:,:,n) = imgaussfilt3(V(:,:,:,n), sigma, 'FilterSize', filterSize);
                    end
                end

                % nuisance factor removal
                if ~isempty(handles.nuisance)
                    disp(['apply nuisance factor removal (' handles.nuisance ')']);
                    Xn = [];
                    if contains(handles.nuisance, 'gm')
                        Sd = getNuisanceMeanTimeSeries(V, [], [], []);
                        Xn = [Xn, Sd(:,1)];
                    end
                    if contains(handles.nuisance, 'gs') || contains(handles.nuisance, 'csf') || contains(handles.nuisance, 'wm') 
                        % get Nuisance time-series (Global Mean, Global Signal, CSF, WM)
                        Sd = getNuisanceMeanTimeSeries(V, csfV, wmV, brV);
                        if contains(handles.nuisance, 'gs'), Xn = [Xn, Sd(:,2)]; end
                        if contains(handles.nuisance, 'csf'), Xn = [Xn, Sd(:,3)]; end
                        if contains(handles.nuisance, 'wm'), Xn = [Xn, Sd(:,4)]; end
                    end
                    if contains(handles.nuisance, 'acomp')
                        % get Nuisance time-series (aCompCor)
                        Sd = Xn;
                        aComp = getNuisanceaCompCor(V, csfV, wmV, Sd, 99, handles.compNum);
                        Xn = [Sd, aComp];
                    end
%                    figure; imagesc([Xn], [-4, 12]); colorbar;
                    V = getNuisanceRegressionOut(V, Xn, brV); % nuisance factor removal within brain mask
                end

                % ROI time-series from voxel data
                if ismask
                    Z = reshape(V,[],size(V,4));
                    X = Z(maskIdx,:);
                else
                    X = getRoiTSFromNifti4D(V, atlasV, 'mean');
                end
                X = X - nanmean(X,2);

                % Seed ROI time-series
                S = [];
                if handles.full > 0
                    % nothing to do
                elseif handles.compSeed > 0
                    S = getCompTSFromNifti4D(V, seedV);
                    S = S - nanmean(S,2);
                else
                    for n=1:size(seedV,4) % might be 4D ROI type
                        ss = getRoiTSFromNifti4D(V, seedV(:,:,:,n), 'mean');
                        ss = ss - nanmean(ss,2);
                        S = [S; ss];
                    end
                end
                X = X'; S = S'; % for corr function

                % high pass filter as preprocessing step (M.W.Woolrich, 2001)
                if handles.highpass > 0
                    disp(['apply highpass filter (' num2str(handles.highpass) ' Hz)']);
                    X = highpass(X,handles.highpass,1/TR);
                end

                % output X matrix
                if handles.noCache == 0
                    save([cachename '.mat'],'X','TR','-v7.3'); fx.X = X;
                    if ~isempty(S)
                        X=S; save([cachename '-' seedName '.mat'],'X','TR','-v7.3'); fs.X = X;
                    end
                end
            else
                fx = load([cachename '.mat']);
                if exist([cachename '-' seedName '.mat'],'file')
                    fs = load([cachename '-' seedName '.mat']);
                end
            end

            % show processed signals
            if handles.showSig > 0
                figure; plot(fx.X);
                title(['ROI Signals : ' strrep(name,'_','-')]);
                xlabel('Time Series');
                ylabel('Signal Value');
            end

            % show processed signals
            if handles.showRas > 0
                figure; imagesc(fx.X.');
                title(['Raster plot of ROI Signals : ' strrep(name,'_','-')]);
                xlabel('Time Series');
                ylabel('Node number');
                colorbar;
            end

            if handles.noCache == 0
                CX{end+1} = [cachename '.mat'];
                if handles.full == 0, CS{end+1} = [cachename '-' seedName '.mat']; end
            else
                CX{end+1} = fx.X; TR = fx.TR;
                if handles.full == 0, CS{end+1} = fs.X; end
            end
            names{end+1} = name;
        end
    end

    idStr = [atlasName '-' seedName '-r' num2str(handles.rmframe) smstr nstr '-' names{1}];

    % 1st and 2nd level non-directional SCA
    if handles.nondir > 0
        ndscaT2 = [handles.outpath '/nondir-' idStr '.mat'];
        if exist(ndscaT2,'file')
            % load T2 volumes
            disp(['loading ' ndscaT2]);
            load(ndscaT2);
        else
            delete(gcp('nocreate')); % shutdown pools
            if handles.poolnum > 0, parpool(handles.poolnum); end % set pool num
            [B2, RSS2, T2, df, R] = calcSeedCorrMixed(CX, CS, handles.cachepath, idStr);

            % currently recels estimation is supported only for full mask mode.
            recels = {}; FWHMs = {};
            if ismask 
                if isempty(CS) % full mode
                    T = triu(ones(nodeNum,nodeNum,'single'));
                    idx = uint32(find(T==1));
                    T(idx) = RSS2;
                    RSS2 = T' + triu(T,1);
                    R2 = zeros(nodeNum,nodeNum,length(CX),'single');
                    for i=1:length(CX)
                        R3 = zeros(nodeNum,nodeNum,'single');
                        R3(idx) = R(:,i);
                        R3 = R3' + triu(R3,1);
                        R2(:,:,i) = R3;
                    end
                    R=R2; clear R2; clear T; clear R3;
                else % seed mode
                    RSS2 = reshape(RSS2,nodeNum,seedNum);
                    R = reshape(R,nodeNum,seedNum,[]);
                end
                for i=1:length(CX)
                    [recels{i}, FWHMs{i}] = estimateSmoothFWHM(squeeze(R(:,i,:)), RSS2(:,i), df, atlasV);
                end
            end
        
            % output T2 matrix
            disp(['saving ' ndscaT2]);
            save(ndscaT2,'df','T2','recels','FWHMs','-v7.3');
        end
    
        T2(T2>1e+3) = Inf; % treat as Inf
        T2(T2<-1e+3) = -Inf; % treat as -Inf

        figure; imagesc(log10(abs(T2))); colorbar;
        title(['non-directional SCA result (log10(|T-value|) matrix) of ' names{1}]);
        xlabel('Seed ROIs'); ylabel('Target voxels'); colorbar;
    end

    % 1st and 2nd level directional SCA
    if handles.dir > 0
        dscaZ2 = [handles.outpath '/dir-' idStr '.mat'];
        if exist(dscaZ2,'file')
            % load Z2 volumes
            disp(['loading ' dscaZ2]);
            load(dscaZ2);
        else
            % find optimal spot time lag (3-4 sec for human)
            if handles.lags == 0
                if ~contains(atlasFile, 'human') || TR == 0
                    handles.lags = 1; % unknown species
                else 
                    handles.lags = floor(4 / TR); % human case
                end
            end

            % calc directional SCA
            delete(gcp('nocreate')); % shutdown pools
            if handles.poolnum > 0, parpool(handles.poolnum); end % set pool num
            [G2, RSS2, ~, df, Z2] = calcSeedPGCMixed(CX, CS, handles.lags, true, handles.cachepath, idStr, handles.rankmeth);

            % output Z2 matrix
            disp(['saving ' dscaZ2]);
            save(dscaZ2,'df','Z2','-v7.3');
        end
    
        figure; imagesc(Z2); colorbar;
        title(['directional SCA result (Z2 matrix) of ' names{1}]);
        xlabel('Source (Seed) voxels'); ylabel('Target voxels'); colorbar;
    end
end

