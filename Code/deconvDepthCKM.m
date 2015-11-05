function [] = deconvDepthCKM(img, img180, map, recon)
%------------------------------------------------------------%
% This method of obtaining the depth map uses two convolved images.
% It deconvolves them with a certain PSF at a specific defocus
% level and then calculates the distance between the two points.
%------------------------------------------------------------%

NOISE = 0.0006;
W20 = 0:0.5:3; 
maxDefocus = size(W20, 2);
NoPts = 870;
XYrange = 0.05;
R = 0.02;
CKMfit = load('../CKM/CKMdata/CKMfit.mat'); CKMfit = CKMfit.CKMfit;
plotsX = []; plotsY = []; defocusVals = [];
%------------------------------------------------------------%
% Use dilation and thresholding to get the blobs containing the
% peaks. Use centroid of each blob to segment a region around
% the true peak.
%------------------------------------------------------------%
diskelem = strel('disk', 10);
closed = imdilate(recon, diskelem);

thresh = multithresh(closed) / 2;
BW = (imquantize(closed, thresh) - 1);

[L, num] = bwlabel(BW);
rp = regionprops(L, 'Centroid');
figure; imshow(BW, [])
figure; imshow(closed, [])
    

L = 15;
%------------------------------------------------------------%
% Locate the points and their corresponding depths.
% Deconvolve with PSF's corresponding to that defocus and
% measure the distance between the points.
%------------------------------------------------------------%
for i = 1 : num

    %------------------------------------------------------------%
    % Locate the points and extract regions. Then find their
    % position and depth (approximate).
    %------------------------------------------------------------%
    C = round(rp(i).Centroid);

    % Extract the region
    region = recon(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
    % Maxima in the region will be the position of the peak
    peakInt = max(region(:));
    [X, Y] = find(region == peakInt);

    sz = size(region, 1); CC = [X, Y];
    coords = round([C(2) - CC(2) + sz/2, C(1) - CC(1) + sz/2]);

    % Intensity
    I = recon(coords(1), coords(2));
    % Get the depth
    depth = map(coords(1), coords(2));
    
    %------------------------------------------------------------%
    % For each point use this depth to deconvolve the full images.
    %------------------------------------------------------------%
    [psf, psf180] = CPMpsf(XYrange, NoPts, R, depth);

    % Deconvolve
    decon = wienerCustom(img, psf, NOISE);
    decon180 = wienerCustom(img180, psf180, NOISE);
    
    % Extract the regions
    region = decon(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
    region180 = decon180(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
    figure; imshowpair(region, region180)
    
    %------------------------------------------------------------%
    % Fit a Gaussian using least squares.
    %------------------------------------------------------------%
    [n, m] = size(region); [X, Y] = meshgrid(1:n, 1:m);
    options = optimset('TolX', 1e-20); 
    peak = max(region(:)); [rc, cc] = find(peak == region);
    peak180 = max(region180(:)); [rc180, cc180] = find(peak180 == region180);

    guess = [peak, rc(1), cc(1), 1];
    LB = [0, 1, 1, -20]; UB = [peak, n, n, 20];
    guess180 = [peak180, rc180(1), cc180(1), 1];
    LB180 = [0, 1, 1, -20]; UB180 = [peak180, n, n, 20];

    % least square fit
    params = lsqnonlin(@(P) objfun2(P, X, Y, region), ...
                       guess, LB, UB, options);
    params180 = lsqnonlin(@(P) objfun2(P, X, Y, region180), ...
                          guess180, LB180, UB180, options);
    hold on; plot(params(3), params(2), '*')
    hold on; plot(params180(3), params180(2), '*')
    
    
    %------------------------------------------------------------%
    % Find the distance between the peaks.
    %------------------------------------------------------------%
    shift  = [params(1), params(2)];
    shift180 = [params180(1), params180(2)];
    dist = sqrt((shift(1) - shift180(1)) ^ 2 + ...
                (shift(2) - shift180(2)) ^ 2)    
    
    %------------------------------------------------------------%
    % Find the depth information, which acts as a correction for
    % the initial depth used for deconvolution. Also localize
    % in 3D.
    %------------------------------------------------------------%
    initialDepth = depth
    
    defocus = feval(CKMfit, dist + initialDepth)
    defocusVals = [defocusVals; defocus];
    
    
    plotsX = [plotsX, coords(2)];
    plotsY = [plotsY, coords(1)];
end

defocusVals

end