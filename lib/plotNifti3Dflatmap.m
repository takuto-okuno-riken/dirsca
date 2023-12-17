%%
% Plotting NIfTI volume by flat mapping
% input:
%  V            nifti 3D volume (X x Y x Z) or index volume
%  maskV        3D mask volume (X x Y x Z)
%  isFullVoxel  full voxel atlas or not for 3D mask volume
%  flatXY       plotting flatmap plane index (index x 2)
%  gRate        integer grouping rate for coloring (default: 10)
%  crange       plot color range (default: [min, max])
%  cmap         color map for V (default: hot)
%  zeroColor    group 'zero' color (default: [])
%  bkColor      background color (default: [0 0 0])

function crange = plotNifti3Dflatmap(V, maskV, isFullVoxel, flatXY, gRate, crange, cmap, zeroColor, bkColor)
    if nargin < 9, bkColor = [0 0 0]; end
    if nargin < 8, zeroColor = []; end
    if nargin < 7, cmap = hot; end
    if nargin < 6, crange = []; end
    if nargin < 5, gRate = 10; end

    if isempty(crange)
        crange = [min(V(:)),max(V(:))];
    end
    zlen = length(flatXY);

    if size(V,2) == 1 && size(V,3) == 1
        % flat index input
        Vidx = V;
    else
        % 3D volume input
        if isFullVoxel == 1
            mIdx = find(maskV(:) > 0);
            Vidx = V(mIdx);
        else
            Vidx = getRoiTSFromNifti4D(V, maskV, 'mean');
        end
    end

    if zlen ~= length(Vidx)
        disp('3D volume (masked) size is not equal to flatXY size.')
        return;
    end
    Vidx(isnan(Vidx)) = 0;
    if crange(1) > 0
        Vidx(Vidx<crange(1)) = 0;
    else
        Vidx(Vidx<crange(1)) = crange(1);
    end
    Vidx(Vidx>crange(2)) = crange(2);
    Vidx = floor(gRate * Vidx);

    % change window size
    pt2 = get(gcf,{'Position'});
    pt = pt2{1};
    rect1 = [ pt(1), pt(2) - (pt(3)-pt(4)), pt(3), pt(3)]; % same h and h
    set(gcf, 'Position', rect1);

    % calc color map
    uni = unique(Vidx);
    clen = length(cmap);
    clr = zeros(length(uni),3,'single');
    for i=1:length(uni)
        col = floor((uni(i) - gRate * crange(1)) / (gRate * (crange(2) - crange(1))) * clen);
        if col < 1, col = 1; end
        if col > clen, col = clen; end
        clr(i,:) = cmap(col,:);
    end

    % show plot
    if ~isempty(zeroColor), clr(1,:) = zeroColor; % might be 0 group
    else, clr(1,:) = clr(2,:); end
    colormap(cmap);
    gscatter(flatXY(:,1),flatXY(:,2),Vidx,clr);
    set(gca,'Color',bkColor);
    legend('off'); daspect([1 1 1]);
    xlim([min(flatXY(:,1))*1.1 max(flatXY(:,1))*1.1]);
    ylim([min(flatXY(:,2))*1.1 max(flatXY(:,2))*1.1]);
end
