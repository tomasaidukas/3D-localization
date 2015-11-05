function [imgavg, depthavg] = CKM(img, img180, original, W20, NoPts, XYrange,...
                                  R, sigmaRange, NOISE)
%------------------------------------------------------------%
% Performs CKM reconstruction technique on two images of the same scene,
% but different defocus values. It needs to be supplied with the image
% parameters/dimensions and a depth map W20 obtained by previous method of
% phase correlation and distance between the points
%------------------------------------------------------------%

% Median filter the images
% img = double(medfilt2(img, [3,3]));
% img180 = double(medfilt2(img180, [3,3]));


%------------------------------------------------------------%
% Generate a range of PSF with different defocus values
%------------------------------------------------------------%
[psf, psf180] = CPMpsf(XYrange, NoPts, R, 0);
maxDefocus = size(W20, 1);
maxSigma = size(sigmaRange, 2);

defocusmap = zeros(NoPts, NoPts, maxSigma);
stackofPSF = zeros(NoPts, NoPts, maxDefocus);
stackofPSF180 = zeros(NoPts, NoPts, maxDefocus);

for i = 1:maxDefocus
    
    % NOTE THAT i=1 is the defocus strength at 0!
    [psf, psf180] = CPMpsf(XYrange, NoPts, R, W20(i));

    stackofPSF(:, :, i) = psf;
    stackofPSF180(:, :, i) = psf180;
end

%------------------------------------------------------------%
% Use Wiener filter to deconvolve the images with PSF's
% from the stack.
%------------------------------------------------------------%
restImg = zeros(NoPts, NoPts, maxDefocus);
restImg180 = zeros(NoPts, NoPts, maxDefocus);

for i = 1:maxDefocus
    psf = stackofPSF(:, :, i);
    psf180 = stackofPSF180(:, :, i);

    decon = wienerCustom(img, psf, NOISE);
    decon180 = wienerCustom(img180, psf180, NOISE);

    restImg(:, :, i) = decon;
    restImg180(:, :, i) = decon180;
end


%------------------------------------------------------------%
% Filter each difference with a different standard 
% deviation of the Gaussian filter
%------------------------------------------------------------%
for sigma = 1 : maxSigma

    S = sigmaRange(sigma);
    gaussfilt = fspecial('gaussian',[2*ceil(2*S)+1 2*ceil(2*S)+1], S);
    filtdif = zeros(NoPts, NoPts, maxDefocus);

    for i = 1:maxDefocus
        subtr = restImg(:, :, i) - restImg180(:, :, i);
        squared = abs(subtr);
        differ = double(imfilter(squared, gaussfilt, 'conv'));
        filtdif(:, :, i) = differ;
    end



    %------------------------------------------------------------%
    % From every filtered difference image select a neighborhood
    % and find which one has the minimum value
    % Loop through all reconstructed image differences
    % and look which reconstruction produced the best
    % result for each defocusing value
    %------------------------------------------------------------%
    for i = 1 : NoPts
        for j = 1 : NoPts
            temp = 999;
            for w = 1:maxDefocus
                % Compare the pixels/neighbourhoods
                % for each reconstruction difference and find which one
                % gives the smallest pixel value. This corresponds to the best
                % defocus value.
                pix = filtdif(i, j, w);
                
                if pix < temp
                    temp = pix;
                    bestW = w;
                end
            end

            defocusmap(i, j, sigma) = bestW;
        end
    end

    %------------------------------------------------------------%
    % Perform a spatially varying deconvolution.
    %------------------------------------------------------------%
    for i = 1 : NoPts
        for j = 1 : NoPts
            str = defocusmap(i, j, sigma);
            reconIMG(i, j, sigma) = restImg(i, j, str);
            reconIMG180(i, j, sigma) = restImg180(i, j, str);
        end
    end
end


%------------------------------------------------------------%
% Take the average between different images for different Gaussian
% kernel sizes. Also compute the average depth map used for
% reconstruction.
%------------------------------------------------------------%
imgavg = 0;
for i = 1:maxSigma
    imgavg = imgavg + reconIMG(:,:,i) + reconIMG180(:,:,i);
end
imgavg = double(imgavg ./ (2*maxSigma));


depthavg = 0;
for i = 1:maxSigma
    depthavg = depthavg + defocusmap(:,:,i);
end
depthavg = round(depthavg ./ (maxSigma));

end
