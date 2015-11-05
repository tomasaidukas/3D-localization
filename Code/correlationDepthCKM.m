function [defocusVals] = correlationDepthCKM(img, img180, W20, NoPts, XYrange, R, S)
%------------------------------------------------------------%
% This method uses phase correlation to find the locations
% of the points to be reconstructed and then segments regions
% around them. These regions are phase correlated to obtain the 
% displacement between the PSF and conjugate PSF. This displacement
% is then used to find the defocus by using previously calculated
% depth map passed to this function.
%------------------------------------------------------------%

    depthMap = zeros(size(img));
    CKMfit = load('../CKM/CKMdata/CKMfit.mat'); 
    CKMfit = CKMfit.CKMfit;
    distances = []; defocusVals = [];
    maxDefocus = size(W20, 2);

    %------------------------------------------------------------%
    % Get the sum of correlations with different defocus values
    % and use them for obtaining the approximate location of the
    % points.
    %------------------------------------------------------------%
    sumCorr = 0; sumCorr180 = 0;
    for i = 1:maxDefocus

        [psf, psf180] = CPMpsf(XYrange, NoPts, R, W20(i));
        
        [Corr, ac] = PhaseCorrelation(img, psf);
        [Corr180, ac180] = PhaseCorrelation(img180, psf180);
        sumCorr = sumCorr + Corr;
        sumCorr180 = sumCorr180 + Corr180;
    end
    
    %------------------------------------------------------------%
    % Threshold the image and use binary closing to obtain blobs
    % containing points of interest.
    %------------------------------------------------------------%
    thresh = multithresh(sumCorr)*2; 
    binary = (imquantize(sumCorr, thresh) - 1);
    diskelem = strel('disk', 3);
    binary = imclose(binary, diskelem);
    figure; imshow(sumCorr, [])
    figure; imshow(binary, [])

    
    %------------------------------------------------------------%
    % Label each blob, find the centroid and segment a region
    % around the centroid to find the regions containing PSF's
    %------------------------------------------------------------%
    [lbl, num] = bwlabel(binary);
    props = regionprops(lbl, 'Centroid');
    
    %------------------------------------------------------------%
    % Measure distance between PSF's.
    %------------------------------------------------------------%
    for i = 1 : num
        % Location of the peak.
        loc = props(i).Centroid; loc = round(fliplr(loc)); K = 70; L = 70;        
        region = img(loc(1)-L:loc(1)+K, ...
                     loc(2)-L:loc(2)+K);
        region180 = img180(loc(1)-K:loc(1)+L, ...
                           loc(2)-K:loc(2)+L);
        
        %------------------------------------------------------------%
        % Interpolate
        %------------------------------------------------------------%
        sampling = 20;
        region = imresize(region, sampling);
        region180 = imresize(region180, sampling);
    
        %------------------------------------------------------------%
        % Correlate
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
        shift = [row1 - row2 , col1 - col2] ./ sampling + (K - L);
        % Find the Euclidean distance
        if shift(1) > 0 | shift(2) > 0
            dist = -sqrt(shift(1)^2 + shift(2)^2);
        else
            dist = sqrt(shift(1)^2 + shift(2)^2);
        end
        
        
        %------------------------------------------------------------%
        % Obtain the defocus for each point
        %------------------------------------------------------------%
        defocus = feval(CKMfit, dist);
        if (defocus < max(W20)) & (defocus > -max(W20))
            defocusVals = [defocusVals; defocus];
        end
    end
end
