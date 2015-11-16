clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up the images
img = double(imread('test.jpg'));
img180 = double(imread('test180.jpg'));

% The original images are now perfectly aligned
% due to PSF shifting
figure; imshowpair(img180, img)
% figure; imshow((img180 - img), [])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a range of PSF with different defocus values
[psf, psf180] = psfdefocus(0);
[xdim, ydim] = size(psf);
stackofPSF = zeros(xdim, ydim, 5);
stackofPSF180 = zeros(xdim, ydim, 5);

% Use the PSF at zero defocus as a reference for
% cross-correlation dimension shifting, since
% this was used to create the images. (Non conjugate)
shiftRef = psf;

for i = 1:5
    
    % NOTE THAT i=1 is the defocus strength at 0!!!!!!!!!!!!!!!!!!!!!!!!
    [psf, psf180] = psfdefocus(i-1);
    
    stackofPSF(:, :, i) = cortrans(shiftRef, psf);
    stackofPSF180(:, :, i) = cortrans(shiftRef, psf180);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Wiener filter to deconvolve the images with PSF's
% from the stack.
restImg = zeros(size(img, 1), size(img, 2), 5);
restImg180 = zeros(size(img180, 1), size(img180, 2), 5);

for i = 1:5
    psf = stackofPSF(:, :, i);
    psf180 = stackofPSF180(:, :, i);
    
    estimated_snr = 0.1;
%     decon = double(deconvwnr(img, psf, estimated_snr));
%     decon180 = double(deconvwnr(img180, psf180, estimated_snr));
    decon = double(deconvlucy(img, psf));
    decon180 = double(deconvlucy(img180, psf180));
    restImg(:, :, i) = decon;
    restImg180(:, :, i) = decon180;
    
%     figure; imshow(restImg(:, :, i), [])
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take neighbourhoods from each image and check which patch
% has the minimized value for each PSF
% 1. Subtract both imags
% 2. Apply a Gaussian low pass filter
% 3. Take a small neighbourhood from the difference
% 4. Create a map of the errors for each pixel
% 5. Repeat for other PSF
% 6. Once each patch has the best defocus value estimate
%    use Wiener filter to reconstruct each patch

gaussfilt = fspecial('gaussian', [5 5], 1);
filtdif = zeros(size(img, 1), size(img, 2), 5);

for i = 1:5
    subtr = restImg(:, :, i) - restImg180(:, :, i);
    squared = subtr .^ 2;
    
%     figure; imshow(squared, [])
    filtdif(:, :, i) = double(imfilter(squared, gaussfilt, 'conv'));
%     filtdif(:, :, i) = double(squared);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now compare the neighbourhoods of the images
img = filtdif(:, :, 1);
indices = size(img, 1) * size(img, 2);
defocusmap = double(zeros(size(img, 1), size(img, 2)));
minmap = double(zeros(size(img, 1), size(img, 2)));

minmap(:) = 999999999999;
defocusmap(:) = 999;
% Loop through all reconstructed image differences
% and look which reconstruction produced the best
% result for each defocusing value

% For each image pair...
for i = 1:4
    img1 = filtdif(:, :, i);
    for j = (i+1):5
        % Compare the pixels/neighbourhoods
        % for each reconstruction and find which is
        % the best for each value
        img2 = filtdif(:, :, j);

        % Subtract the two images
        % -ve values indicate that img1 has smaller values
        % +ve values indicate that img2 has smaller values
        subtr = img1 - img2;

        % All values bigger than 0 (+ve) indicate that img2
        % has smaller values. Set the -ve values to 0.
        % Use find() to get the non-zero elements
        copySubtr = subtr;
        copySubtr(copySubtr > 0) = 0;
        [x2 y2, num2] = find(copySubtr);
        [index2] = sub2ind(size(copySubtr), x2, y2);
        
        copySubtr = subtr;
        copySubtr(copySubtr < 0) = 0;
        [x1 y1, num1] = find(copySubtr);
        [index1] = sub2ind(size(copySubtr), x1, y1);
        
        for idx = index1
            if img1(idx) < minmap(idx)
                minmap(idx) = img1(idx);
                defocusmap(idx) = i;
            end
        end

        for idx = index2
            if img2(idx) < minmap(idx)
                minmap(idx) = img2(idx);
                defocusmap(idx) = j;
            end
        end
    end
end

defocusmap(defocusmap == 999) = 0;
figure; imshow(defocusmap, []), colorbar

for i = 1:5
    temp = filtdif(:, :, i);
    sum(temp(:))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform a spatially varying deconvolution.
% Use different PSF's to deconvolve the full images
% Then combine these images with a certain weighting factor
img = double(imread('test.jpg'));
img180 = double(imread('test180.jpg'));

unique(defocusmap)

% xdim = size(img, 1);
% ydim = size(img, 2);
% strs = unique(defocusmap);
% 
% foravg = zeros(xdim, ydim, size(strs, 1));
% foravg180 = zeros(xdim, ydim, size(strs, 1));
% 
% The zeroth element does not count
% length = size(strs, 1) - 1;
% 
% for i = 1:length
%     
%     if strs(i) ~= 0
%         % Loop through the best PSF's
%         psf = stackofPSF(:, :, strs(i));
%         psf180 = stackofPSF180(:, :, strs(i));
% 
%         % Use these PSF's for deconvolution
%         foravg(:, :, i) = deconvlucy(img, psf);
%         foravg180(:, :, i) = deconvlucy(img180, psf180);
% %         foravg(:, :, i) = double(deconvwnr(img, psf, estimated_snr));
% %         foravg180(:, :, i) = double(deconvwnr(img180, psf180, estimated_snr));
%     end
% end
% 
% imagesum = 0;
% for i = 1:length
%     imagesum = imagesum + foravg(:, :, i) + foravg180(:, :, i);
% end
% 
% avgimg = imagesum ./ (2*length);
% figure; imshow(avgimg, [])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Ideal low pass filter
% M = size(avgimg, 1); N = size(avgimg, 2);
% [X, Y] = meshgrid(linspace(-4, 4, N), linspace(-4, 4, M));
% [theta, p] = cart2pol(X, Y);
% lowpass = double(p <= 1);
% 
% % Apply the filter
% avgimgFFT = fftshift(fft2(avgimg));
% filtered = real(ifft2(fftshift(double(lowpass .* avgimgFFT))));
% figure; imshow(double(filtered), [])