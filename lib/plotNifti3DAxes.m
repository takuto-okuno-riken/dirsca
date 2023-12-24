%%
% Plotting NIfTI volume
% input:
%  V            nifti 3D volume (X x Y x Z)
%  operation    operation for each plane ('max'(default), 'min','mode','mean','median')
%  range        plot color range (default: [min, max])
%  cmap         color map for V (default: parula)
%  backV        3D background volume (default: [])
%  backcmap     color map for backV (default: gray)
%  isRtoL       X axis is right to left (default: false)
%  bRate        background color rate (default: 0.3)
%  fRate        front color rate (default: 0.8)

function plotNifti3DAxes(V, operation, range, cmap, backV, backcmap, isRtoL, bRate, fRate)
    if nargin < 9, fRate = 0.8; end
    if nargin < 8, bRate = 0.3; end
    if nargin < 7, isRtoL = false; end
    if nargin < 6, backcmap = gray; end
    if nargin < 5, backV = []; end
    if nargin < 4, cmap = parula; end
    if nargin < 3, range = [min(V(:)),max(V(:))]; end
    if nargin < 2, operation = 'max'; end

    if isnan(range(1)) && isnan(range(2))
        range = [0 1];
    end

    % cmap size check
    if size(cmap,1) ~= 256
        % resampling
        step = (size(cmap,1)-1) / 256.0;
        idx = 1:step:size(cmap,1);
        idx = floor(idx);
        cmap = cmap(idx(1:256),:);
    end

    switch(operation)
    case 'mode'
        YZ = mode(V,1);
        XZ = mode(V,2);
        XY = mode(V,3);
    case 'min'
        YZ = nanmin(V,[],1);
        XZ = nanmin(V,[],2);
        XY = nanmin(V,[],3);
    case 'max'
        YZ = nanmax(V,[],1);
        XZ = nanmax(V,[],2);
        XY = nanmax(V,[],3);
    case 'mean'
        YZ = nanmean(V,1);
        XZ = nanmean(V,2);
        XY = nanmean(V,3);
    case 'median'
        YZ = nanmedian(V,1);
        XZ = nanmedian(V,2);
        XY = nanmedian(V,3);
    end
    % rotate
    YZ = rot90(squeeze(YZ)); % from back (right should be right)
    XZ = rot90(squeeze(XZ));
    if isRtoL, XZ = fliplr(XZ); end
    XY = rot90(squeeze(XY)); % from top (right should be bottom)
    if isRtoL, XY = fliplr(XY); end
    % fixed color range to [0 1]
    YZ = (YZ - range(1)) / (range(2) - range(1)); YZ(YZ<0) = 0; YZ(YZ>1) = 1; YZ(isnan(YZ)) = 0;
    XZ = (XZ - range(1)) / (range(2) - range(1)); XZ(XZ<0) = 0; XZ(XZ>1) = 1; XZ(isnan(XZ)) = 0;
    XY = (XY - range(1)) / (range(2) - range(1)); XY(XY<0) = 0; XY(XY>1) = 1; XY(isnan(XY)) = 0;
    % gray to index map
    [YZi,~] = gray2ind(YZ,256);
    [XZi,~] = gray2ind(XZ,256);
    [XYi,~] = gray2ind(XY,256);

    if ~isempty(backV)
        YZb = nanmax(backV,[],1);
        XZb = nanmax(backV,[],2);
        XYb = nanmax(backV,[],3);
        % rotate
        YZb = rot90(squeeze(YZb)); % from back (right should be right)
        XZb = rot90(squeeze(XZb));
        if isRtoL, XZb = fliplr(XZb); end
        XYb = rot90(squeeze(XYb)); % from top (right should be bottom)
        if isRtoL, XYb = fliplr(XYb); end
        % fixed color range to [0 0.33]
        mb = single(nanmax(backV,[],'all'));
        YZb = single(YZb) / mb; YZb(YZb<0.01) = 0; YZb(isnan(YZb)) = 0;
        XZb = single(XZb) / mb; XZb(XZb<0.01) = 0; XZb(isnan(XZb)) = 0;
        XYb = single(XYb) / mb; XYb(XYb<0.01) = 0; XYb(isnan(XYb)) = 0;
        % gray to index map
        [YZbi,~] = gray2ind(YZb,256);
        [XZbi,~] = gray2ind(XZb,256);
        [XYbi,~] = gray2ind(XYb,256);
    end

    % change window size
    pt2 = get(gcf,{'Position'});
    pt = pt2{1};
    rect1 = [ pt(1), pt(2) - (pt(3)-pt(4)), pt(3), pt(3)]; % same h and h
    set(gcf, 'Position', rect1);

    % show plot 
    rect = [size(XYi,2)/(size(YZi,2)+size(XYi,2)), ...
        size(XYi,1)/(size(XYi,1)+size(XZi,1)),...
        size(YZi,2)/(size(YZi,2)+size(XZi,2)),...
        size(YZi,1)/(size(YZi,1)+size(XYi,1))];
    axes('Position', rect);
    if ~isempty(backV)
        imb = ind2rgb(YZbi,backcmap);
        imf = ind2rgb(YZi,cmap);
        imshow(imb * bRate + imf * fRate);
    else
        imshow(ind2rgb(YZi,cmap));
    end

    rect = [0, ...
        size(XYi,1)/(size(XYi,1)+size(XZi,1)),...
        size(XZi,2)/(size(YZi,2)+size(XZi,2)),...
        size(XZi,1)/(size(XZi,1)+size(XYi,1))];
    axes('Position', rect);
    if ~isempty(backV)
        imb = ind2rgb(XZbi,backcmap);
        imf = ind2rgb(XZi,cmap);
        imshow(imb * bRate + imf * fRate);
    else
        imshow(ind2rgb(XZi,cmap));
    end

    rect = [0, 0,...
        size(XYi,2)/(size(XYi,2)+size(YZi,2)),...
        size(XYi,1)/(size(XYi,1)+size(XZi,1))];
    axes('Position', rect);
    if ~isempty(backV)
        imb = ind2rgb(XYbi,backcmap);
        imf = ind2rgb(XYi,cmap);
        imshow(imb * bRate + imf * fRate);
    else
        imshow(ind2rgb(XYi,cmap));
    end

    % color bar
    rect = [size(XYi,2)/(size(YZi,2)+size(XYi,2))+0.1, 0.05,...
        (size(YZi,2)/(size(YZi,2)+size(XZi,2))) / 2,...
        size(XYi,1)/(size(XYi,1)+size(XZi,1))-0.15];
    axes('Position', rect);
    colormap(cmap); caxis(range);
    set(gca,'visible','off');
    colorbar;
end
