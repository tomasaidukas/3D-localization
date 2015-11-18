function psf = simulateDHPSF(NoPts, XYrange, R, W20)
%------------------------------------------------------------%
% Create a DHPSF psf at a certain defocus value.
%------------------------------------------------------------%

[x,y] = meshgrid(linspace(-XYrange, XYrange, NoPts),...
                 linspace(-XYrange, XYrange, NoPts));
[tet,p]=cart2pol(double(x),double(y));

% DH phase mask, normalized
U = load('U.mat');
U = U.U;

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

% figure; imshow(psf, [])
% figure; imshow(U .* (p<=R), [])
end

