% Estimation of smoothness based on [residual] images
% refs.
% Matlab script of spm_est_smoothness.m (in SPM12)
%
% S.J. Kiebel, J.B. Poline, K.J. Friston, A.P. Holmes, and K.J. Worsley.
% Robust Smoothness Estimation in Statistical Parametric Maps Using
% Standardized Residuals from the General Linear Model. NeuroImage,
% 10:756-766, 1999.
% input:
%  R          all residual matrix (voxels x frames)
%  RSS        Residual Sum of Squares in each voxel (voxels x 1)
%  df         degree of freedom
%  maskV      3D mask volume

function [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, maskV)
    % init
    dim = 3; % 3 dimension
    frames = size(R,2);
    % frame sampling .. this may not be necesarry
    nRes = min(frames, 64); % resolution
    iRes = round(linspace(1,frames,nRes))'; % get indices of sampling frames

    L = zeros(size(R,1),dim,dim);
    ssq = zeros(size(R,1),1);
    V = single(maskV); V(:) = nan;
    idx = find(maskV > 0);
    for t=1:nRes % each resolution
        % get residual image
        Rt = R(:,iRes(t));
        sqRMS = sqrt(RSS/df); % sqrt(RMS) = sqrt(r'*r/edf)
        V(idx) = Rt ./ sqRMS; % set values to volume

        % calc graduate each direction (edge voxel will have nan.)
        [Vdx,Vdy,Vdz] = calcGrad__(V);
        dx = Vdx(idx); dy = Vdy(idx); dz = Vdz(idx);
        % sum of squares (to remove nan, later)
        ssq = ssq + dx.*dx + dy.*dy + dz.*dz;
        % covariance of finite differences
        L(:,1,1) = L(:,1,1) + dx.*dx;
        L(:,1,2) = L(:,1,2) + dx.*dy;
        L(:,2,2) = L(:,2,2) + dy.*dy;
        L(:,1,3) = L(:,1,3) + dx.*dz;
        L(:,2,3) = L(:,2,3) + dy.*dz;
        L(:,3,3) = L(:,3,3) + dz.*dz;
    end
    % In terms of standardized residuals e/sqrt(RMS) this is
    %  \hat\Lambda = (1/DF) * grad(e/sqrt(RMS))' * grad(e/sqrt(RMS))
    L = L / nRes; % Average
    L = L * (frames / df); % Scale
    reselXYZ = [L(:,1,1) L(:,2,2)  L(:,3,3)];
    reselImg = L(:,1,1).*L(:,2,2).*L(:,3,3) + ...
               L(:,1,2).*L(:,2,3).*L(:,1,3)*2 - ...
               L(:,1,1).*L(:,2,3).*L(:,2,3) - ...
               L(:,1,2).*L(:,1,2).*L(:,3,3) - ...
               L(:,1,3).*L(:,2,2).*L(:,1,3);
    reselImg(reselImg<0) = 0;
    % Convert det(Lambda) and diag(Lambda) to units of resels
    reselImg = sqrt(reselImg/(4*log(2))^dim);
    reselXYZ = sqrt(reselXYZ/(4*log(2)));
    % RESEL estimator and Global equivalent FWHM
    idx = isnan(ssq) | ssq < sqrt(eps);
    reselImg = mean(reselImg(~idx,:));
    reselXYZ = mean(reselXYZ(~idx,:));

    recel = reselImg^(1/dim)*(reselXYZ/prod(reselXYZ)^(1/dim));
    FWHM  = full(sparse(1,1:dim,1./recel,1,3));
    FWHM(isnan(FWHM)) = 0;
    FWHM(~FWHM) = 1;
end

function [Vdx,Vdy,Vdz] = calcGrad__(V)
    Vdx = V; Vdy = V; Vdz = V;
    Vdx(:) = nan; Vdy(:) = nan; Vdz(:) = nan;
    Vdx(2:end-1,:,:) = (V(3:end,:,:) - V(1:end-2,:,:)) / 2;
    Vdy(:,2:end-1,:) = (V(:,3:end,:) - V(:,1:end-2,:)) / 2;
    Vdz(:,:,2:end-1) = (V(:,:,3:end) - V(:,:,1:end-2)) / 2;
end
