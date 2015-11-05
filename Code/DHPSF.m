function [img2D, DHPSFmap, plotX, plotY] = DHPSF(img, angles, PSFs, angle2defocus)
%------------------------------------------------------------%
% DHPSF method as developed by Standford university people.
% Locate peaks using correlation -> segment them out ->
% -> fit a 2D double Gaussian -> measure angle and distance
% between the peaks -> use a map to get defocus and point position.
%------------------------------------------------------------%

%------------------------------------------------------------%
% Filter noise. These are usually just corrupted pixels 1 pixel wide,
% while the PSF peaks of interest are much larger than that. Median
% filter removes the single pixels.
%------------------------------------------------------------%
% imgF = double(medfilt2(img, [3,3]));
imgF = img;
img2D = zeros(size(img)); depths = zeros(size(img)); midPts = [];
DHPSFmap = []; plotX = []; plotY = [];

%------------------------------------------------------------%
% Perform Phase correlation to locate the rough peak positions
%------------------------------------------------------------%
corrstack = 0;
for i = 1 : size(PSFs, 3)

    psf = PSFs(:, :, i);
    % Low pass filter in the frequency domain
    h = fspecial('gaussian', [6, 6], 1.5);
    H = fftshift(fft2(h, size(img, 1), size(img, 2)));

    % Find the Fourier transform of both images
    FTpsf = fftshift(fft2(psf));
    FTimg = fftshift(fft2(imgF));
    
    % Power spectrum
    CC = double(FTimg .* conj(FTpsf) .* H);
    AC = double(conj(FTimg) .* (FTimg) .* H);

    % Cross-Correlation
    cc = abs(ifftshift(ifft2(CC)));

    corrstack = corrstack + cc;
end

% Threshold the image
thresh = multithresh(corrstack);
binary = (imquantize(corrstack, thresh) - 1);

% Use morphological operations to obtain connected regions, which
% contain the peaks.
diskelem = strel('disk', 6);
closed = imdilate(binary, diskelem);
figure; imshow(closed)

% Segment these peaks 
[L, num] = bwlabel(closed);
rp = regionprops(L, 'Centroid');
K = 25;

%------------------------------------------------------------%
% For each peak pair find their midpoints and angles to indicate the
% depth by fitting a 2D double Gaussian.
%------------------------------------------------------------%
for i = 1 : num
    prop = rp(i);
    C = prop.Centroid;

    % Extract the region
    box = imgF(C(2)-K:C(2)+K, C(1)-K:C(1)+K);
    
    % Least square fitting routine for a double Gaussian
    % Arrays used for result storage
    peakCoords = []; c = size(psf) / 2;
    boxC = size(box, 1);

    %------------------------------------------------------------%
    % Initial guess for the centre locations using Hough circle
    % detector
    %------------------------------------------------------------%
    [centers, radii] = imfindcircles(box, [1 10],'Sensitivity', 1);
    % Take non-overlaping circles
    temp1 = centers(1, :); temp2 = centers(2, :);
    xc1 = temp1(1); yc1 = temp1(2); xc2 = temp2(1); yc2 = temp2(2);

    for k = 2 : size(radii, 1)
        temp1 = centers(1, :); temp2 = centers(k, :);
        xc1 = temp1(1); yc1 = temp1(2); xc2 = temp2(1); yc2 = temp2(2);
        vect1 = xc1 - xc2; vect2 = yc1 - yc2;

        if sqrt(vect1^2 + vect2^2) >= (radii(1) + radii(k))
            xc1 = temp1(1); yc1 = temp1(2);
            xc2 = temp2(1); yc2 = temp2(2);
            break
        end
    end
    
    
    %------------------------------------------------------------%
    % Initial guess for the peak heights using the peak values at the
    % centroid positions. Fit the Gaussian.
    %------------------------------------------------------------%

    [n, m] = size(box); [X, Y] = meshgrid(1:n);
    options = optimset('TolX', 1e-20);
    [n, m] = size(box); 
    
    % guess [normalization, xc, yc, sigma,
    %        normalization, xc, yc, sigma]
    guess = [max(box(:)), xc1, yc1, 1, ...
             max(box(:)), xc2, yc2, 1];
    LB = [max(box(:))/4, 1, 1, -5, ...
          max(box(:))/4, 1, 1, -5];
    UB = [max(box(:)), n, n, 5, ...
          max(box(:)), n, n, 5];

    % least square fit
    params = lsqnonlin(@(P) objfun(P, X, Y, box), guess, LB, UB, options)

    % Shift the absolute co-ordinates
    coords1 = [C(1) - boxC / 2 + params(2), ...
               C(2) - boxC / 2 + params(3)];
    coords2 = [C(1) - boxC / 2 + params(6), ...
               C(2) - boxC / 2 + params(7)];        

    % Plot peaks superimposed inside the box for testing
    midpt = [params(2) + params(6), params(3) + params(7)] ./ 2;
    
    figure('units', 'normalized', 'position', [0 0 1 1]); 
    imshow(box, [])
    hold on; plot(params(2), params(3), '*')
%             hold on; plot(params(6), params(7), '*')
    viscircles([centers(1,:); centers(k,:)], ...
               [radii(1); radii(k)], 'EdgeColor','b');    
    hold on; plot(midpt(1), midpt(2), '*')
    hold on; plot([midpt(1), params(2)], [midpt(2), params(3)])
    hold on; plot([midpt(1), params(2)], [midpt(2), midpt(2)])

    %------------------------------------------------------------%
    % Take each peak pair and find the distance between them as well as
    % the angle they make with the horizontal axis to determine the
    % depth. The midpoint will be the localized co-ordinate in 2D.
    %------------------------------------------------------------%
    ang = atand((params(2) - params(6)) / (params(3) - params(7)))
    midpt = [coords1(1) + coords2(1), coords1(2) + coords2(2)] ./ 2;

    %------------------------------------------------------------%
    % Store the angle between the pair corresponding to each defocus
    % value. Also store the midpoint between them, which corresponds to
    % the localised 2D position of the particle.
    %------------------------------------------------------------%
    angles = [angles; ang];
    midPts = [midPts; midpt];
    
    %------------------------------------------------------------%
    % Use the created map to map angle -> defocus.
    %------------------------------------------------------------%
    defocus = feval(angle2defocus, ang);
    DHPSFmap = [DHPSFmap ; defocus];

    img2D(round(midpt(2)), round(midpt(1))) = 1;
    plotX = [plotX, midpt(1)]; plotY = [plotY, midpt(2)];
end
DHPSFmap
end

