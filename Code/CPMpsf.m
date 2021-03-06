function [psf, psf180] = CPMpsf(XYrange, NoPts, R, str)
%------------------------------------------------------------%
% Create PSF's using a cubic phase mask.
%------------------------------------------------------------%

[X,Y] = meshgrid(linspace(-XYrange, XYrange, NoPts),...CPM = double(exp(1i .* 2 .* pi .* alpha .* (X.^3 + Y.^3)));
                 linspace(-XYrange, XYrange, NoPts));
[tet, p] = cart2pol(X, Y);
apt = double(p <= R);
p = apt .* p;
p = p ./ max(p(:));

defocus = double(exp(1i .* 2 .* pi .* str .* p.^2));

alpha = 4.0; % Peak aberration
CPM = double(exp(1i .* 2 .* pi .* alpha .* (X.^3 + Y.^3)));

% figure; imshow(p, [])
% figure; imshow(CPM, [])
% figure; imshow(CPM .* apt, [])
% figure; imshow(apt, [])

pupilFunc = apt .* CPM .* defocus;
pupilFunc180 = conj(pupilFunc);

psf = fftshift(fft2(pupilFunc)) .* conj(fftshift(fft2(pupilFunc)));
psf180 = fftshift(fft2(pupilFunc180)) .* conj(fftshift(fft2(pupilFunc180)));

psf = psf ./ sum(psf(:));
psf180 = psf180 ./ sum(psf180(:));

end
