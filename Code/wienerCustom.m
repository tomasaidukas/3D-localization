function restored = wienerCustom(img, psf, K)
%------------------------------------------------------------%
% Wiener filter.
%------------------------------------------------------------%
    
FT = fftshift(fft2(img));
OTF = fftshift(fft2(psf));

OTF(OTF == 0) = 0.0000000001;
modOTF = OTF .* conj(OTF);
restored = FT .* conj(OTF) ./ (modOTF + K);
restored = abs(ifftshift(ifft2(restored)));

end