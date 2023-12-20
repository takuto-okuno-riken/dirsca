%%
% Plot functional flat mapping command line tool

function flatmap(varargin)

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
    handles.inFiles = {};
    handles.atlasFile = 'data/human_2mm_cubeRoi2.nii.gz';
    handles.maskFile = [];
    handles.flatmapFile = 'data/human_ffm_cubeRoi2.mat';
    handles.range = [];
    handles.cmap = '';
    handles.backdot = [];
    handles.backgr = [0 0 0];

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'--atlas'}
                handles.atlasFile = varargin{i+1};
                i = i + 1;
            case {'--mask'}
                handles.maskFile = varargin{i+1};
                i = i + 1;
            case {'--flatmap'}
                handles.flatmapFile = varargin{i+1};
                i = i + 1;
            case {'--range'}
                handles.range = [str2num(varargin{i+1}) str2num(varargin{i+2})];
                i = i + 2;
            case {'--cmap'}
                handles.cmap = varargin{i+1};
                i = i + 1;
            case {'--backdot'}
                handles.backdot = [str2num(varargin{i+1}) str2num(varargin{i+2}) str2num(varargin{i+3})];
                i = i + 3;
            case {'--backgr'}
                handles.backgr = [str2num(varargin{i+1}) str2num(varargin{i+2}) str2num(varargin{i+3})];
                i = i + 3;
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
                    handles.inFiles = [handles.inFiles varargin{i}];
                end
        end
        i = i + 1;
    end
    
    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif isempty(handles.inFiles)
        disp('input NIfTI files are required. please specify mat files.');
        showUsage();
        return;
    elseif isempty(handles.atlasFile) && isempty(handles.maskFile)
        disp('Atlas or mask NIfTI file is required. please specify nifti files.');
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
    
    disp(['usage: ' exeName ' [options] file1.nii.gz ...']);
    disp('  --atlas file        NIfTI <file> of target ROI voxels (default: data/human_2mm_cubeRoi2.nii.gz)');
    disp('  --mask file         NIfTI <file> of target mask (full voxel) (Preferred over atlas option)');
    disp('  --flatmap file      functional flat map definition <file>(default: data/human_ffm_cubeRoi2.mat)');
    disp('  --range min max     T-value/Z-score color range for image plot');
    disp('  --cmap type         color map type for image plot (default: "hot")');
    disp('  --cmap type         color map type for image plot (default: "hot")');
    disp('  --backdot r g b     background dot color (default: [])');
    disp('  --backgr r g b      background color (default: 0 0 0)');
    disp('  -v, --version       show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;

    % load target voxels (cube ROI or mask)
    if ~isempty(handles.flatmapFile) && ~exist(handles.flatmapFile,'file')
        disp(['functional flat map definition is not found : ' handles.maskFile]);
        return;
    end
    disp(['load functional flat map definition : ' handles.flatmapFile]);
    load(handles.flatmapFile);

    % load target voxels (cube ROI or mask)
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

    disp(['load voxel mask : ' atlasFile]);
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
    isRtoL = false; % volume is adjusted already

    % set colormap
    switch handles.cmap
    case 'hot', cmap = hot;
    case 'jet', cmap = jet;
    case 'turbo', cmap = turbo;
    otherwise, cmap = parula;
    end

    % read input files
    for i = 1:length(handles.inFiles)
        % load result mat files
        argv = handles.inFiles{i};
        flist = dir(argv);
        if isempty(flist)
            disp(['file is not found. ignoring : ' argv]);
            continue;
        end
        for k=1:length(flist)
            fname = [flist(k).folder filesep flist(k).name];
            [path,name,ext] = fileparts(fname);
            name = strrep(name,'.nii','');

            disp(['read NIfTI file : ' name]);
            info = niftiinfo(fname);
            V = single(niftiread(info));
            V = adjustVolumeDir(V,info);
            V(isnan(V)) = 0;

            figure; crange = plotNifti3Dflatmap(V, atlasV, ismask, reduction, 10, handles.range, cmap, handles.backdot, handles.backgr);
            title([name]); colorbar;
            disp(['color range=[' num2str(crange(1)) ' ' num2str(crange(2)) ']']);
        end
    end
end

