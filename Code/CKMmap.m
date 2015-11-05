function [distance, defocus] = CKMmap(img, img180, original, W20, NoPts, XYrange,...
                     R, sigmaRange, NOISE)
%------------------------------------------------------------%
% Use phase correlation to
% measure the displacement between the PSF's in the normal
% and conjugate images. (convovled with a conjugate PSF)
%------------------------------------------------------------%
maxDefocus = size(W20, 2);
maxSigma = size(sigmaRange, 2);
distance = []; defocus = [];

for i = 1:maxDefocus

    %------------------------------------------------------------%
    % A method which fits a 2D Gaussian to the PSF's and measures
    % the displacement.
    % ------------------------------------------------------------%
%     % Get the PSF's and close up to the center peak
%     [psf, psf180] = CPMpsf(XYrange, NoPts, R, W20(i));    
%     c = NoPts / 2;
%     psf = psf(c-20:c+20, c-20:c+20);
%     psf180 = psf180(c-20:c+20, c-20:c+20);
% 
%     % Fit a Gaussian fn.
%     [n, m] = size(psf); [X, Y] = meshgrid(1:n, 1:m);
%     options = optimset('TolX', 1e-20); 
%     peak = max(psf(:)); [xc, yc] = find(peak == psf);
%     peak180 = max(psf180(:)); [xc180, yc180] = find(peak180 == psf180);
% 
%     guess = [peak, xc(1), yc(1), 1];
%     LB = [0, 1, 1, -20]; UB = [peak, n, n, 20];
%     guess180 = [peak180, xc180(1), yc180(1), 1];
%     LB180 = [0, 1, 1, -20]; UB180 = [peak180, n, n, 20];
% 
%     % least square fit
%     params = lsqnonlin(@(P) objfun2(P, X, Y, psf), guess, LB, UB, options);
%     params180 = lsqnonlin(@(P) objfun2(P, X, Y, psf180), guess180, LB180, UB180, options);
%     % Shift away from the centre
%     shift = [params(2), params(3)];
%     shift180 = [params180(2), params180(3)];
%     % Distance      
%     dist = sqrt((shift(1) - shift180(1)) ^ 2 + ...
%                 (shift(2) - shift180(2)) ^ 2)
% 
%     figure; imshow(psf, [])
%             hold on
%             plot(shift(1), shift(2), '*')
%             hold on
%             plot(shift180(1), shift180(2), '*')
    
    %------------------------------------------------------------%
    % Get the PSF's for a range of defocus values and segment them out
    %------------------------------------------------------------%
    % Get the PSF's
    [psf, psf180] = CPMpsf(XYrange, NoPts, R, W20(i));    
    % Center location and box crop parameters
    loc = [NoPts / 2, NoPts / 2]; K = 50; L = 50;
    % Extract the regions
    region = psf(loc(2)-L:loc(2)+K, ...
                 loc(1)-L:loc(1)+K);
    region180 = psf180(loc(2)-K:loc(2)+L, ...
                       loc(1)-K:loc(1)+L);
    
    %------------------------------------------------------------%
    % Interpolate regions for higher precision. 'sampling' defines
    % by how much the image size should be increased
    %------------------------------------------------------------%
    sampling = 20;
    region = imresize(region, sampling);
    region180 = imresize(region180, sampling);

    %------------------------------------------------------------%
    % Perform correlation and auto-correlation to obtain the
    % shift between the PSF's
    %------------------------------------------------------------%
    FT1 = fftshift(fft2(region)); FT2 = fftshift(fft2(region180));
    % Power spectrum
    CC = double(FT1 .* conj(FT2));
    AC = double(FT1 .* conj(FT1));
    % Cross-Correlation
    cc = abs(ifftshift(ifft2(CC)));
    ac = abs(ifftshift(ifft2(AC)));
    % Peak locations
    [row1, col1] = find(max(cc(:)) == cc);
    [row2, col2] = find(max(ac(:)) == ac);
    % K and L should be the same. They were used for testing.
    shift = [row1 - row2 , col1 - col2] ./ sampling + (K - L);
    % Find the Euclidean distance.
    if shift(1) > 0
        dist = -sqrt(shift(1)^2 + shift(2)^2)
    else
        dist = sqrt(shift(1)^2 + shift(2)^2)
    end

    defocus = [defocus; W20(i)];
    distance = [distance; dist];
end

%------------------------------------------------------------%
% Plot the results and use linear interpolation to fit 
% the data. This function is then used for the CKM technique.
%------------------------------------------------------------%
figure; plot(defocus, distance)    
ylabel('distance'); xlabel('W20')

CKMfit = fit(distance, defocus, 'linearinterp');
fitted = feval(CKMfit, distance);

figure; plot(distance, defocus, '*')
        hold on
        plot(distance, fitted)
        xlabel('distance'); ylabel('W20')

save('../CKM/CKMdata/CKMfit.mat', 'CKMfit')
save('../CKM/CKMdata/distance.mat', 'distance');

end
