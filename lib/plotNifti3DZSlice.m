%%
% Plotting NIfTI volume
% input:
%  V            nifti 3D volume (X x Y x Z)
%  zidx         plotting Z plane index
%  crange       plot color range (default: [min, max])
%  cmap         color map for V (default: parula)
%  backV        3D background volume 
%  backcmap     color map for backV (default: gray)
%  isRtoL       X axis is right to left (default: true)
%  bRate        background color rate (default: 0.3)
%  fRate        front color rate (default: 0.8)

function plotNifti3DZSlice(V, zidx, crange, cmap, backV, backcmap, isRtoL, bRate, fRate)
    if nargin < 9, fRate = 0.8; end
    if nargin < 8, bRate = 0.3; end
    if nargin < 7, isRtoL = true; end
    if nargin < 6, backcmap = gray; end
    if nargin < 5, backV = []; end
    if nargin < 4, cmap = parula; end
    if nargin < 3, crange = [min(V(:)),max(V(:))]; end
    if nargin < 2, zidx = [1:size(V,3)]; end

    if isnan(crange(1)) && isnan(crange(2))
        crange = [0 1];
    end
    zlen = length(zidx);
    step = ceil(sqrt(zlen));

    % change window size
    pt2 = get(gcf,{'Position'});
    pt = pt2{1};
    rect1 = [ pt(1), pt(2) - (pt(3)-pt(4)), floor(pt(3)*(size(V,1)/size(V,2))), pt(3)]; % same h and h
    set(gcf, 'Position', rect1);

    % show plot
    for i=1:zlen
        z = zidx(i);

        XY = V(:,:,z);
        % rotate
        XY = rot90(squeeze(XY)); % from top (right should be bottom)
        if isRtoL, XY = fliplr(XY); end
        % fixed color range to [0 1]
        XY = (XY - crange(1)) / (crange(2) - crange(1)); XY(XY<0) = 0; XY(XY>1) = 1; XY(isnan(XY)) = 0;
        % gray to index map
        [XYi,~] = gray2ind(XY,256);

        if ~isempty(backV)
            XYb = backV(:,:,z);
            % rotate
            XYb = rot90(squeeze(XYb));  % from top (right should be bottom)
            if isRtoL, XYb = fliplr(XYb); end
            % fixed color range to [0 0.33]
            mb = single(nanmax(backV,[],'all'));
            XYb = single(XYb) / mb; XYb(XYb<0.01) = 0; XYb(isnan(XYb)) = 0;
            % gray to index map
            [XYbi,~] = gray2ind(XYb,256);
        end

        % show plot 
        rect = [mod((i-1),step)/step,...
            1-ceil(i/step)/step,...
            1.02/step, 1.02/step];
        axes('Position', rect);
        if ~isempty(backV)
            imb = ind2rgb(XYbi,backcmap);
            imf = ind2rgb(XYi,cmap);
            imshow(imb * bRate + imf * fRate);
        else
            imshow(ind2rgb(XYi,cmap));
        end
    end
%{
    % color bar
    rect = [size(YZi,1)/(size(YZi,1)+size(XYi,1))+0.1, 0.05,...
        (size(YZi,2)/(size(YZi,2)+size(XZi,2))) / 2,...
        size(XYi,1)/(size(XYi,1)+size(XZi,1))-0.15];
    setPlot3DAxis(rect);
    colormap(cmap); caxis(crange);
    set(gca,'visible','off');
    colorbar;
%}
end
