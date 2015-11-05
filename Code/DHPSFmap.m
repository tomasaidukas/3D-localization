function [PSFs, angles, W20fit] = DHPSFmap(NoPts, XYrange, R, W20range)
%------------------------------------------------------------%
% Create a map that maps angle between the lobes of the DH-PSF
% to defocus (depth). Same procedure as described in DHPSF is used.
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
% 
%------------------------------------------------------------%
% Load the experimental SLM representing the phase mask.
%------------------------------------------------------------%
U = load('U.mat');
U = U.U;
[refVect, p1, p2] = referenceVector(NoPts, XYrange, R, W20range);
referenceMidPt = [];

%------------------------------------------------------------%
% Loop through a range of defocus values and simulate the DH-PSF in
% order to apply different defocus aberration and record its peak
% locations using a 2D Gaussian fit using least squares. This fits a
% double Gaussian for a given region containing two peaks of the
% DH-PSF.
%------------------------------------------------------------%
indices = 1:size(W20range, 2);  angles = []; midPts = []; oldang = 0;
for index = indices

    W20 = W20range(index);

    % Get the PSF from the measured SLM stored in U
    [tet,p]=cart2pol(x,y);
    p = p./R;
    pup = U .* (p<=1) .* exp(1i .* 2 .* pi .* W20 .* p.^2);
    % pup = exp(1i.*angle(U)).*(p<=1).*(exp(1i.*2.*pi.*W20.*p.^2));
    psf = fftshift(fft2(pup)).*conj(fftshift(fft2(pup)));
    psf = psf ./ sum(psf(:));

    % Segment the region around the center that contains the peaks
    c = size(psf) / 2;  L = 25;
    box = psf(round(c(2))-L:round(c(2))+L,...
              round(c(1))-L:round(c(1))+L);


    %------------------------------------------------------------%
    % Least square fitting routine for a double Gaussian
    %------------------------------------------------------------%
    peaks = zeros(size(psf, 1), size(psf, 2));
    peakCoords = []; boxpeaks = zeros(L*2 + 1, L*2 + 1); 
    boxC = size(box, 1);

    %------------------------------------------------------------%
    % Initial guess for the centre locations using Hough circle
    % detector
    %------------------------------------------------------------%
    [centers, radii] = imfindcircles(box, [1 10],'Sensitivity', 1);
    
    % Take non-overlaping circles
    temp1 = centers(1, :); temp2 = centers(2, :);
    xc1 = temp1(1); yc1 = temp1(2); xc2 = temp2(1); yc2 = temp2(2);

    for k = 2 : size(radii, 1)
        temp1 = centers(1, :); temp2 = centers(k, :);
        xc1 = temp1(1); yc1 = temp1(2); xc2 = temp2(1); yc2 = temp2(2);
        vect1 = xc1 - xc2; vect2 = yc1 - yc2;

        if sqrt(vect1^2 + vect2^2) >= (radii(1) + radii(k))
            xc1 = temp1(1); yc1 = temp1(2);
            xc2 = temp2(1); yc2 = temp2(2);
%                 figure; imshow(box, [])
%                         viscircles([centers(1,:); centers(k,:)], ...
%                                    [radii(1); radii(k)], 'EdgeColor','b');
            break
        end
    end
    
    %------------------------------------------------------------%
    % Fit the Gaussian.
    %------------------------------------------------------------%
    [n, m] = size(box); [X, Y] = meshgrid(1:(1/5):n);

    [XC, YC] = meshgrid(1:m);
    box = interp2(XC, YC, box, X, Y);
    options = optimset('TolX', 1e-20);

    % guess [normalization, xc, yc, sigma,
    %        normalization, xc, yc, sigma]
    guess = [max(box(:)), xc1, yc1, 1, ...
             max(box(:)), xc2, yc2, 1];
    LB = [max(box(:))/4, 1, 1, -5, ...
          max(box(:))/4, 1, 1, -5];        
    UB = [max(box(:)), n, n, 5, ...
          max(box(:)), n, n, 5];

    % Least square fit.
    params = lsqnonlin(@(P) objfun(P, X, Y, box), guess, LB, UB, options);

    % Shift the absolute co-ordinates.
    coords1 = [c(1) - boxC / 2 + params(2), ...
              c(2) - boxC / 2 + params(3)];
    coords2 = [c(1) - boxC / 2 + params(6), ...
               c(2) - boxC / 2 + params(7)];
    midpt = [params(2) + params(6), params(3) + params(7)] ./ 2;

%     figure('units', 'normalized', 'position', [0 0 1 1]); 
%     imshow(box, [])
%     hold on; plot(params(2), params(3), '*')
% %             hold on; plot(params(6), params(7), '*')
%     viscircles([centers(1,:); centers(k,:)], ...
%                [radii(1); radii(k)], 'EdgeColor','b');    
%     hold on; plot(midpt(1), midpt(2), '*')
%     hold on; plot([midpt(1), params(2)], [midpt(2), params(3)])
%     hold on; plot([midpt(1), params(2)], [midpt(2), midpt(2)])
    %------------------------------------------------------------%
    % [X, Y] is a vector pointing from one peak to another.
    %------------------------------------------------------------%    
    ang = atand((params(2) - params(6))/ (params(3) - params(7)))

    %------------------------------------------------------------%
    % Store the angle between the pair corresponding to each defocus
    % value. Also store the midpoint between them, which corresponds to
    % the localised 2D position of the particle.
    %------------------------------------------------------------%
    angles = [angles; ang];

    PSFs(:, :, index) = psf;

    oldang = ang;
end

%------------------------------------------------------------%
% This gives back a map: angles -> defocus
% polyfit(x, y, power)
%------------------------------------------------------------%
W20fit = fit(angles, W20range', 'linearinterp');

% Given an angle map it maps them to defocus
% polyval(yFit, x)
angle2defocus = feval(W20fit, angles);

save('./DHPSFdata/ang2defocus.mat', 'W20fit')
save('./DHPSFdata/angles.mat', 'angles');
save('./DHPSFdata/templates.mat', 'PSFs');


figure; plot(angle2defocus, angles, 'r')
        hold on
        plot(W20range, angles, '*')
        
end






