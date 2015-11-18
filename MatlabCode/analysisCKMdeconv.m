function [depth, coords] = analysisCKMdeconv(img, img180, psf, ...
                              psf180, camera, CKMfit, NOISE, sigmaRange)
%------------------------------------------------------------%
% This method of obtaining the depth map uses two convolved images.
% It deconvolves them with a certain PSF at a specific defocus
% level and then calculates the distance between the two points.
%------------------------------------------------------------%

%------------------------------------------------------------%
% Find the point using cross-correlation
%------------------------------------------------------------%
% Find the Fourier transform of both images
FTpsf = fftshift(fft2(psf));
FTimg = fftshift(fft2(img));
% Power spectrum
CC = double(FTimg .* conj(FTpsf));% .* H);
% Cross-Correlation
cc = abs(ifftshift(ifft2(CC)));
[ROW, COL] = find(max(cc(:)) == cc);
C = [COL, ROW];


[W20, maxDefocus, NoPts, XYrange, R, f] = camera{:};
defocusVals = [];
L = 10; samp = 1 ;

%------------------------------------------------------------%
% For each point use this depth to deconvolve the full images.
%------------------------------------------------------------%
% Deconvolve
decon = wienerCustom(img, psf, NOISE);
decon180 = wienerCustom(img180, psf180, NOISE);
% Extract the regions
region = decon(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
region180 = decon180(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
region = imresize(region, samp);
region180 = imresize(region180, samp);
% region(region <= thresh) = 0;
% region180(region180 <= thresh) = 0;
% figure; imshow(region, []);

%------------------------------------------------------------%
% Fit a Gaussian using least squares.
%------------------------------------------------------------%
[n, m] = size(region); [X, Y] = meshgrid(1:n, 1:m);
options = optimset('TolX', 1e-20, 'Display', 'off'); 
peak = max(region(:)); [rc, cc] = find(peak == region);
peak180 = max(region180(:)); [rc180, cc180] = find(peak180 == region180);

guess = [peak, rc(1), cc(1), 1 * samp, 1 * samp];
LB = [0, 1, 1, -20, -20]; UB = [peak, n, n, 5 * samp, 5 * samp];
guess180 = [peak180, rc180(1), cc180(1), 1 * samp, 1 * samp];
LB180 = [0, 1, 1, -20, -20]; UB180 = [peak180, n, n, 5* samp,5* samp];

% least square fit
params = lsqnonlin(@(P) objfun2(P, X, Y, region), ...
                   guess, LB, UB, options);
params180 = lsqnonlin(@(P) objfun2(P, X, Y, region180), ...
                      guess180, LB180, UB180, options);
% 
% figure; imshowpair(region, region180)
% hold on; plot(params(2), params(3), '*')
% hold on; plot(params180(2), params180(3), '*')
%------------------------------------------------------------%
% Find the distance between the peaks.
%------------------------------------------------------------%
shift  = double([params(2), params(3)]) ./ samp;
shift180 = double([params180(2), params180(3)]) ./ samp;
dist = sqrt((shift(1) - shift180(1)) ^ 2 + ...
            (shift(2) - shift180(2)) ^ 2);

%------------------------------------------------------------%
% Find the depth information using the calibration map.
%------------------------------------------------------------%
defocus = feval(CKMfit, dist);

%------------------------------------------------------------%
% Use that defocus in the CKM method.
%------------------------------------------------------------%
camera{1} = defocus; % Defocus value for CKM
[imgCKM, mapCKM] = CKM(img, img180, camera, sigmaRange, NOISE);

%------------------------------------------------------------%
% Localize in 2D using CKM reconstructed images
%------------------------------------------------------------%
[ROW, COL] = find(max(imgCKM(:)) == imgCKM); L = 15;
region = imgCKM(ROW-L:ROW+L, COL-L:COL+L);
% figure; imshow(region, [])

[n, m] = size(region); [X, Y] = meshgrid(1:n, 1:m);
options = optimset('TolX', 1e-20, 'Display', 'off'); 
peak = max(region(:)); [rc, cc] = find(peak == region);
guess = [peak, rc(1), cc(1), 1, 1];
LB = [0, 1, 1, 0, 0]; UB = [peak, n, n, 5, 5];

% least square fit
params = lsqnonlin(@(P) objfun2(P, X, Y, region), ...
                   guess, LB, UB, options);
               
coords = [COL + (n/2 - params(3)) ./ 1, ROW + (m/2 - params(2)) ./ 1];
depth = mapCKM(round(coords(2)), round(coords(1)));
defocusRange = camera{1};
depth = defocusRange(depth);
                
end