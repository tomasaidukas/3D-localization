function [cc, ac] = PhaseCorrelation(img, psf)
%------------------------------------------------------------%
% Performs phase correlation and auto correlation. Returns both
% correlations.
%------------------------------------------------------------%

% Find the Fourier transform of both images
FTpsf = fftshift(fft2(psf));
FTimg = fftshift(fft2(img));

% Power spectrum
CC = double(FTimg .* conj(FTpsf));
AC = double(FTimg .* conj(FTimg));

% Cross-Correlation
cc = abs(ifftshift(ifft2(CC)));
ac = abs(ifftshift(ifft2(AC)));
      
  end
