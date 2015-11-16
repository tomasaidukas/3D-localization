clear all
close all

% Defocus <1
% No noise: M1 - fails
%           M2 - very good (some numerical precision errors)
%           DHPSF flawless
% High SNR: M1 - failure
%           M2 - fails close to W20 = 0
%           DHPSF very good

% Defocus >1 <2
% No noise: M1 - pretty good
%           M2 - very good (some numerical precision errors)
%           DHPSF slightly worse than M2 due to poor map
% High SNR: M1 - failure
%           M2 - good, better than DHPSF
%           DHPSF - fails due to poor map

% Need to fix DHPSF map. Can't distinguish which peak is which. After a
% certain W20 value the map gets reversed and does not follow the expected
% trend.

%------------------------------------------------------------%
% Load simulation data
%------------------------------------------------------------%
img = load('rand.mat');
img = img.finalimg;
img180 = load('randC.mat');
img180 = img180.finalimg180;
original = load('randO.mat');
original = original.original;
imgDH = load('randDH.mat');
imgDH = imgDH.finalimgDH;
NOISE = 0.0006;
sigmaRange = 10:15;

figure; imshow(img, [])
% figure; imshow(imgDH, [])

%------------------------------------------------------------%
% Parameters of the image/camera and initial defocus values 
% to be used for map generation.
%------------------------------------------------------------%
W20 = 0:0.01:3; 
maxDefocus = size(W20, 2);
NoPts = 870;
XYrange = 0.05; % from -5 mm to 5 mm
R = 0.03; % 2 mm           
f = 0.1; % 100 mm
% NOISE = 0.00001 % 0.00001 perfect for the no noise image
NOISE = 0.0005;
camera = {W20; maxDefocus; NoPts; XYrange; R; f};

% [ref, ref180] = CPMpsf(XYrange, NoPts, R, 1);
% [psf, psf180] = CPMpsf(XYrange, NoPts, R, 1);    
% % Deconvolve
% decon = wienerCustom(psf, ref, NOISE);
% figure; imshow(decon, [])
% asdasdsa


%------------------------------------------------------------%
% Generates depth maps for CKM and DH-PSF.
%------------------------------------------------------------%
% [PSFs, angles, W20fit] = DHPSFmap(camera);
angles = load('./DHPSFdata/angles.mat'); angles = angles.angles;
PSFs = load('./DHPSFdata/templates.mat'); PSFs = PSFs.PSFs;
W20fit = load('./DHPSFdata/ang2defocus.mat'); angle2defocus = W20fit.W20fit;

% [dist, defoc] = CKMmap(camera, NOISE);      % Problem with this method is the available precision. At very small defocus the distance is not measurable even by interpolating the image to make it 20 times larger.

% Crop size does not matter - interpolation does, hence select a lot of
% points.
[dist, defoc] = CKMmap2(camera, NOISE);    % This method is sensitive to noise and the SNR parameter used for deconvolution. This leads to a distorted deconvolved points with some systematic errors and the centroid(found by Gaussian fitting) vary slightly - especially for small defocus values. So, if incorrect method is used to deconvolve the points, the obtained results are not very precise.
asdasd
%------------------------------------------------------------%
% DH-PSF localisation.
%------------------------------------------------------------%
% Load up the reference vector according to which the angles were
% measured. This vector corresponds to the angle at 0 defocus and is a
% unit vetor.
[DHPSF2D, depthDHPSF, Xdhpsf Ydhpsf] = ...
                    DHPSF(imgDH, angles, PSFs, angle2defocus);

%------------------------------------------------------------%
% CKM localisation.
%------------------------------------------------------------%
camera = {0:0.6:2; maxDefocus; NoPts; XYrange; R; f};
% Get the approximate point locations
[imgCKM1, mapCKM1] = CKM(img, img180, camera, sigmaRange, NOISE);
                       
% depths = correlationDepthCKM(img, img180, imgCKM1, camera)
% camera = {depths; maxDefocus; NoPts; XYrange; R; f};
% [imgCKM, mapCKM] = CKM(img, img180, camera, sigmaRange, NOISE);
% [CKM3D, CKM2D, depthCKM, Xckm, Yckm] = localCKM(mapCKM, imgCKM, depths);


% optimize if it works well (pull out the psf calc before the loop)
depths = deconvDepthCKM(img, img180, mapCKM1, imgCKM1, camera, NOISE)
camera = {depths; maxDefocus; NoPts; XYrange; R; f};
[imgCKM, mapCKM] = CKM(img, img180, camera, sigmaRange, NOISE);
[CKM3D, CKM2D, depthCKM2, Xckm, Yckm] = localCKM(mapCKM, imgCKM, depths);

%------------------------------------------------------------%
% Visualize
%------------------------------------------------------------%
originalBW = imregionalmax(original); [row, col] = find(originalBW);
figure; imshow(original, [])
        hold on; plot(col, row, '.r')
        hold on; plot(Xckm, Yckm, '.b')
        hold on; plot(Xdhpsf, Ydhpsf, '.g')

%------------------------------------------------------------%
% Display depth values to check 3D localisation
%------------------------------------------------------------%
% depthCKM
depthCKM2
depthDHPSF