function [distance, defocus] = CKMmap(camera, NOISE)
%------------------------------------------------------------%
% Use phase correlation to
% measure the displacement between the PSF's in the normal
% and conjugate images. (convovled with a conjugate PSF)
%------------------------------------------------------------%
[W20, maxDefocus, NoPts, XYrange, R, f] = camera{:};

distance = []; defocus = [];

for i = 1:maxDefocus    
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
    shift = [row1 - row2 , col1 - col2] ./ sampling;
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
% %------------------------------------------------------------%
% figure; plot(defocus, distance)    
% ylabel('distance'); xlabel('W20')
% 
% CKMfit = fit(distance, defocus, 'linearinterp');
% fitted = feval(CKMfit, distance);
% 
% CKMinvfit = fit(defocus, distance, 'linearinterp');
% ayylmao = feval(CKMinvfit, defocus);
% 
% figure; plot(distance, defocus, '*')
%         hold on
%         plot(distance, fitted)
%         xlabel('distance'); ylabel('W20')
% 
% save('../CKM/CKMdata/CKMfit.mat', 'CKMfit')
% save('../CKM/CKMdata/distance.mat', 'distance');
% save('../CKM/CKMdata/CKMinvfit.mat', 'CKMinvfit')

CKMfit = fit(defocus, distance, 'poly2')

figure; plot(defocus, distance)    
ylabel('distance'); xlabel('W20')

fitted = feval(CKMfit, defocus);
figure; plot(defocus, distance, '*')
        hold on
        plot(defocus, fitted)
        ylabel('distance'); xlabel('W20')
        
save('../CKM/CKMdata/CKMfit1.mat', 'CKMfit')

end
