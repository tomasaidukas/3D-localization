function psf = simulateDHPSF(NoPts, XYrange, R, W20)
%------------------------------------------------------------%
% Create a DHPSF psf at a certain defocus value.
%------------------------------------------------------------%

n = 1; %NB: effective 'n' of mode is: (n - abs(m))/2
m = 1;
Z = 0;
lambda = 550*10^(-9); %wavelength
Wo = 4.18*10^(-3); %beam waist

[x,y] = meshgrid(linspace(-XYrange,XYrange,NoPts),linspace(-XYrange,XYrange,NoPts));
[PHI,RHO]=cart2pol(x,y);
%     
%     Zo = pi*(Wo^2)/lambda; %Ryleigh distance
%     Zh = Z/Zo;
%     Wzh = Wo.*sqrt(1.+Zh.^2);
%     PSIzh = atan(Zh);
%     RHOh = RHO./Wzh;
%     G_rhoh_zh = (Wo./Wzh).*exp(-(RHOh.^2)).*exp(1i.*(RHOh.^2).*Zh).*exp(-1i.*PSIzh);
% 
%     U = zeros(NoPts,NoPts);
% 
%     for u = 1:4
%         Rnm_rhoh = ((sqrt(2).*RHOh).^abs(m)).*Ln_k((n-abs(m))/2,abs(m),2.*RHOh.^2);%Ln_m((n-abs(m))/2,abs(m),2.*RHOh.^2);
%         PHIm_phi = exp(1i.*m.*PHI);
%         Zn_zh = exp(-1i.*n.*PSIzh);
%         Unm = ((G_rhoh_zh.*Rnm_rhoh).*PHIm_phi).*Zn_zh; %GL mode n,m
%         Unm = Unm ./ max(abs(Unm(:)));
%         U = U + Unm;
%         n = n+4;
%         m = m+2;
%     end

U = load('U.mat');
U = U.U;

[tet,p]=cart2pol(x,y);
p = p./R;
pup = U .* (p<=1) .* exp(1i .* 2 .* pi .* W20 .* p.^2);
%pup = exp(1i.*angle(U)).*(p<=1).*(exp(1i.*2.*pi.*W20.*p.^2)); %%PHASE ONLY
psf = fftshift(fft2(pup)).*conj(fftshift(fft2(pup)));
psf = psf ./ sum(psf(:));


end

