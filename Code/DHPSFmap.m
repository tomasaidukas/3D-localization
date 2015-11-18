function [PSFs, angles, W20fit] = DHPSFmap(camera)
%------------------------------------------------------------%
% Create a map that maps angle between the lobes of the DH-PSF
% to defocus (depth). Same procedure as described in DHPSF is used.
%------------------------------------------------------------%

[W20range, maxDefocus, NoPts, XYrange, R, f] = camera{:};
canvas = zeros(NoPts, NoPts);
canvas(NoPts / 2, NoPts / 2) = 1; samp = 1;
store1 = []; store2 = [];

[x,y] = meshgrid(linspace(-XYrange,XYrange,NoPts),linspace(-XYrange,XYrange,NoPts));
[tet,p]=cart2pol(double(x),double(y));

%------------------------------------------------------------%
% Load the experimental SLM representing the phase mask.
%------------------------------------------------------------%
U = load('U.mat');
U = U.U;

%------------------------------------------------------------%
% Loop through a range of defocus values and simulate the DH-PSF in
% order to apply different defocus aberration and record its peak
% locations using a 2D Gaussian fit using least squares. This fits a
% double Gaussian for a given region containing two peaks of the
% DH-PSF.
%------------------------------------------------------------%
indices = 1:size(W20range, 2);  angles = []; midPts = []; oldang = 0;
midpoints = [];
for index = indices

    W20 = W20range(index);
    % Aperture containing ones and zeros
    apt = double(p <= R);
    % pupil must be normalized from 1 to 0
    p = apt .* p;
    p = p ./ max(p(:));

    % Defocus
    deltaW20 = W20 .* p .^2;
    defocus = exp(1i .* 2 .* pi .* deltaW20);

    % Pupil
    pup = U .* apt .* defocus;

    psf = fftshift(fft2(pup)).*conj(fftshift(fft2(pup)));
    psf = psf ./ sum(psf(:));
    
    %------------------------------------------------------------%
    % Store the angle between the pair corresponding to each defocus
    % value. Also store the midpoint between them, which corresponds to
    % the localised 2D position of the particle.
    %------------------------------------------------------------%
    PSFs(:, :, index) = psf;
    
    psf = abs(fftshift(ifft2(fftshift(fft2(psf)) .* fftshift(fft2(canvas)))));
    
    % Segment the region around the center that contains the peaks
    c = size(psf) ./ 2;  L = 25;
    box = psf(round(c(2))-L:round(c(2))+L,...
              round(c(1))-L:round(c(1))+L);
    box = imresize(box, samp);


    %------------------------------------------------------------%
    % Least square fitting routine for a double Gaussian
    %------------------------------------------------------------%
    peaks = zeros(size(psf, 1), size(psf, 2));
    peakCoords = [];
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
%                 figure; imshow(box, [])
%                         viscircles([centers(1,:); centers(k,:)], ...
%                                    [radii(1); radii(k)], 'EdgeColor','b');
            break
        end
    end
    
    %------------------------------------------------------------%
    % Fit the Gaussian.
    %------------------------------------------------------------%
    [n, m] = size(box); [X, Y] = meshgrid(1:n, 1:m);
    options = optimset('TolX', 1e-20);

    % guess [normalization, xc, yc, sigma,
    %        normalization, xc, yc, sigma]
    guess = [max(box(:)), xc1, yc1, 1*samp, ...
             max(box(:)), xc2, yc2, 1*samp];
    LB = [max(box(:))/4, 1, 1, 0, ...
          max(box(:))/4, 1, 1, 0];        
    UB = [max(box(:)), n, n, 5*samp, ...
          max(box(:)), n, n, 5*samp];

    % Least square fit.
    params = lsqnonlin(@(P) objfun(P, X, Y, box), guess, LB, UB, options);

    % Shift the absolute co-ordinates.
    coords1 = [c(1) + (boxC / 2 - params(2)) ./ samp, ...
               c(2) + (boxC / 2 - params(3)) ./ samp];
    coords2 = [c(1) + (boxC / 2 - params(6)) ./ samp, ...
               c(2) + (boxC / 2 - params(7)) ./ samp];
%     midpt = [params(2) + params(6), params(3) + params(7)] ./ 2;
% 
%     figure('units', 'normalized', 'position', [0 0 1 1]); 
%     imshow(box, [])
%     hold on; plot(params(2), params(3), '*')
% %             hold on; plot(params(6), params(7), '*')
%     viscircles([centers(1,:); centers(k,:)], ...
%                [radii(1); radii(k)], 'EdgeColor','b');    
%     hold on; plot(midpt(1), midpt(2), '*')
%     hold on; plot([midpt(1), params(2)], [midpt(2), params(3)])
%     hold on; plot([midpt(1), params(2)], [midpt(2), midpt(2)])
    %------------------------------------------------------------%
    % [X, Y] is a vector pointing from one peak to another.
    %------------------------------------------------------------%    
    par1 = (params(2) - params(6)) ./ samp, par2 = (params(3) - params(7)) ./ samp
    store1 = [store1; par1]; store2 = [store2; par2];
    
    midpt = [coords1(1) + coords2(1), coords1(2) + coords2(2)] ./ 2;
    midpoints = [midpoints; midpt];


    
end

%------------------------------------------------------------%
% This gives back a map: angles -> defocus
% polyfit(x, y, power)
%------------------------------------------------------------%

angles = atand(store2 ./ store1)

W20fit = fit(angles, W20range', 'linearinterp');

% The psf is displaced slightly.
correction = midpoints - NoPts / 2;

% Given an angle map it maps them to defocus
% polyval(yFit, x)
angle2defocus = feval(W20fit, angles);

save('./DHPSFdata/ang2defocus2.mat', 'W20fit')
save('./DHPSFdata/angles.mat', 'angles');
save('./DHPSFdata/templates.mat', 'PSFs');



figure; plot(angle2defocus, angles, 'r')
        hold on
        plot(W20range, angles, '*')

end






