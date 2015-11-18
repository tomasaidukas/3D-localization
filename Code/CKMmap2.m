function [distance, defocus] = CKMmap2(camera, NOISE)
%------------------------------------------------------------%
% Generates a range of PSF's and then deconvolves each PSF
% with a PSF at 0 defocus. This will produce a distance between
% the restored points depending on which PSF was used for
% convolution.
%------------------------------------------------------------%
[W20, maxDefocus, NoPts, XYrange, R, f] = camera{:};
samp = 1; c = NoPts / 2; L = 10; thresh = 0;
distance = []; defocus = [];


[ref, ref180] = CPMpsf(XYrange, NoPts, R, 0);
%------------------------------------------------------------%
% Get the point position using perfect reconstruction and
% calibrate with respect to it.
%------------------------------------------------------------%
for i = 1:maxDefocus
    
    %------------------------------------------------------------%
    % A method which fits a 2D Gaussian to the PSF's and measures
    % the displacement between deconvolved points.
    % ------------------------------------------------------------%
    % Get the PSF's and close up to the center peak
    [psf, psf180] = CPMpsf(XYrange, NoPts, R, W20(i));        

    % Deconvolve
    decon = wienerCustom(psf, ref, NOISE);
    decon180 = wienerCustom(psf180, ref180, NOISE);   
    % Extract the regions
    region = decon(c-L:c+L, c-L:c+L);
    region180 = decon180(c-L:c+L, c-L:c+L);
    region = imresize(region, samp);
    region180 = imresize(region180, samp);
%     region(region <= thresh) = 0;
%     region180(region180 <= thresh) = 0;
%     figure; imshow(region, [])
    
    % Fit a Gaussian fn.
    [n, m] = size(region); [X, Y] = meshgrid(1:n, 1:m);
    options = optimset('TolX', 1e-20, 'Display', 'off'); 
    peak = max(region(:)); [rc, cc] = find(peak == region);
    peak180 = max(region180(:)); [rc180, cc180] = find(peak180 == region180);

    guess = [peak, rc(1), cc(1), 1 * samp, 1 * samp];
    LB = [0, 1, 1, 0, 0]; UB = [peak, n, n, 5 * samp, 5 * samp];
    guess180 = [peak180, rc180(1), cc180(1), 1 * samp, 1 * samp];
    LB180 = [0, 1, 1, 0, 0]; UB180 = [peak180, n, n, 5 * samp, 5 * samp];

    % least square fit
    params = lsqnonlin(@(P) objfun2(P, X, Y, region), guess, LB, UB, options);
    params2 = lsqnonlin(@(P) objfun2(P, X, Y, region180), guess180, LB180, UB180, options)

    % Shift away from the centre
    shift = double([params(2), params(3)]) ./ samp;
    shift2 = double([params2(2), params2(3)]) ./ samp;
    % Distance      
    dist = sqrt((shift(1) - shift2(1)) ^ 2 + ...
                (shift(2) - shift2(2)) ^ 2);
    
%     figure; imshow(region, [])
%     figure; imshowpair(region, region180)
%             hold on
%             plot(shift(1), shift(2), '*')
%             hold on 
%             plot(shift2(1), shift2(2), '*')

    defocus = [defocus; W20(i)];
    distance = [distance; dist];
end

%------------------------------------------------------------%
% Plot the results and use linear interpolation to fit 
% the data. This function is then used for the CKM technique.
%------------------------------------------------------------%

% polynomial fit
figure; plot(defocus, distance)    
ylabel('distance'); xlabel('W20')

% CKMfit = fit(defocus, distance, 'poly2')
% fitted = feval(CKMfit, defocus);
% figure; plot(defocus, distance, '*')
%         hold on
%         plot(defocus, fitted)
%         ylabel('distance'); xlabel('W20')
% save('../CKM/CKMdata/CKMfitPoly.mat', 'CKMfit')

CKMfit = fit(distance, defocus, 'linearinterp');
fitted = feval(CKMfit, distance);

figure; plot(distance, defocus, '*')
        hold on
        plot(distance, fitted)
        xlabel('distance'); ylabel('W20')


save('../CKM/CKMdata/CKMfitLinNOSQRT2.mat', 'CKMfit')

end
