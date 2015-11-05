function [] = CKMshift(img, img180, W20, NoPts, XYrange,...
                                        R, sigmaRange)
%------------------------------------------------------------%
% Plots the shift induced by different defocus values.
%------------------------------------------------------------%
NOISE = 0.0006;
sigmaRange = 10;

% figure; imshow(img, [])
% figure; imshow(imgDH, [])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defocus
W20 = 0:0.5:2;
maxDefocus = size(W20, 2);
% Camera parameters
NoPts = 870;
XYrange = 0.05;
R = 0.02;
% Initialize arrays
maxDefocus = size(W20, 2);
maxSigma = size(sigmaRange, 2);
distance = []; defocus = []; NOISE = 0.00000001;

% Create an image with a random convolved point inside
[refPsf, refPsf180] = CPMpsf(XYrange, NoPts, R, 0);
canvas = zeros(NoPts, NoPts);
randx = randi([100, NoPts - 100], 1, 1);
randy = randi([100, NoPts - 100], 1, 1);
canvas(randx, randy) = 1;
img = abs(fftshift(ifft2(fftshift(fft2(refPsf)) .* fftshift(fft2(canvas)))));
figure; imshow(img, [])

figure; 
for i = 1:maxDefocus
    [psf, psf180] = CPMpsf(XYrange, NoPts, R, W20(i));    
    decon = wienerCustom(img, psf, NOISE);
%         figure; imshow(decon, [])

    % Extract the peak profile
    [row, col] = find(max(decon(:)) == decon);
    region = decon(row-10:row+10, col-10:col+10);
    sampling = 10;
    [n, m] = size(region);
    [XC, YC] = meshgrid(1:m);
    [XI, YI] = meshgrid(1:(1/sampling):m);
    region = interp2(XC, YC, region, XI, YI);
%         figure; imshow(region, [])

    c = round(size(region, 2)/2); profile = region(:, c);
    cc = hsv(12); plot(profile, 'color', cc(i,:))
    hold on
end

end
