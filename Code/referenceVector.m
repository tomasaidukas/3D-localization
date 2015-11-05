function [pairVect, p1, p2] = referenceVector(NoPts, XYrange, R, W20range)
%------------------------------------------------------------%
% Perform the methods as in the DHPSFmap to obtain the angle
% for a vector at 0 defocus. Use this as a reference.
%------------------------------------------------------------%

n = 1; %NB: effective 'n' of mode is: (n - abs(m))/2
m = 1;
Z = 0;
lambda = 550*10^(-9); %wavelength
Wo = 4.18*10^(-3); %beam waist

[x,y] = meshgrid(linspace(-XYrange,XYrange,NoPts),linspace(-XYrange,XYrange,NoPts));
[PHI,RHO]=cart2pol(x,y);

% The experimental SLM
U = load('U.mat');
U = U.U;

% Get the PSF from the measured SLM stored in U
[tet,p]=cart2pol(x,y);
p = p./R;
pup = U .* (p<=1);
psf = fftshift(fft2(pup)).*conj(fftshift(fft2(pup)));
psf = psf ./ sum(psf(:));

% Segment the region around the center that contains the peaks
c = size(psf) / 2;  L = 30;
box = psf(round(c(2))-L:round(c(2))+L,...
          round(c(1))-L:round(c(1))+L);


% Least square fitting routine for a double Gaussian
peaks = zeros(size(psf, 1), size(psf, 2));
peakCoords = []; boxpeaks = zeros(L*2 + 1, L*2 + 1); 
boxC = size(box, 1);

% Initial guess for the centre locations using Hough circle
% detector
[centers, radii] = imfindcircles(box, [1 10],'Sensitivity', 1);
temp1 = centers(1, :); temp2 = centers(2, :);
xc1 = temp1(1); yc1 = temp1(2); xc2 = temp2(1); yc2 = temp2(2);


[n, m] = size(box); [X, Y] = meshgrid(1:n, 1:m);
options = optimset('TolX', 1e-20);

% guess [normalization, xc, yc, sigma,
%        normalization, xc, yc, sigma]
guess = [max(box(:)), xc1, yc1, 1, ...
         max(box(:)), xc2, yc2, 1];
LB = [max(box(:))/4, 1, 1, -5, ...
      max(box(:))/4, 1, 1, -5];
UB = [max(box(:)), n, n, 5, ...
      max(box(:)), n, n, 5];

% least square fit
params = lsqnonlin(@(P) objfun(P, X, Y, box), guess, LB, UB, options);

% Shift the absolute co-ordinates
coords1 = [c(1) - boxC / 2 + params(2), ...
          c(2) - boxC / 2 + params(3)];
coords2 = [c(1) - boxC / 2 + params(6), ...
           c(2) - boxC / 2 + params(7)];

X = coords1(1) - coords2(1); Y = coords1(2) - coords2(2);

% Normalized vector between the lobes
pairVect = [X, Y];
pairLen = sqrt( X^2 + Y^2 );
pairVect = pairVect ./ pairLen
p1 = [coords1(1), coords2(1)];
p2 = [coords1(2), coords2(2)];

end





