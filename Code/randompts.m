function [img, imgConj, canvas, imgDH] = randompts(XYrange, NoPts, R, sz, str, numpts)
%------------------------------------------------------------%
% Generate points with different defocus value
% Their positions are random and must be sparse
%
% dimx, dimy specify the number of pixels along each axis
% R is the radius of the aperture in pixels
% sz is the size of the points
% str is the defocus strength
% numpts is how many random points will be generated
%------------------------------------------------------------%

[blurpsf, blurpsfConj] = CPMpsf(XYrange, NoPts, R, str);    

% ------------------------------------------------------------ %
% Create the image with random points.
%------------------------------------------------------------%
canvas = zeros(NoPts, NoPts);
plotX = []; plotY = [];

for i = 1 : numpts
    randx = randi([50, NoPts - 50], 1, 1);
    randy = randi([50, NoPts - 50], 1, 1);
    pixel = rand;
    canvas(randx:sz+randx, randy:sz+randy) = pixel;
end

%------------------------------------------------------------%
% Convolve the image.
%------------------------------------------------------------%
img = abs(fftshift(ifft2(fftshift(fft2(blurpsf)) .* fftshift(fft2(canvas)))));
imgConj = abs(fftshift(ifft2(fftshift(fft2(blurpsfConj)) .* fftshift(fft2(canvas)))));

blurpsf = simulateDHPSF(NoPts, XYrange, R, str);
imgDH = abs(fftshift(ifft2(fftshift(fft2(blurpsf)) .* fftshift(fft2(canvas)))));

end
