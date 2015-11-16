function [] = makeplots(pupil, imgtit, imgname)
% Takes in the pupil function and find the MTF, PSF and MTF in 2D.
% Also takes in the title and filename for the images to be saved.
%
% The point spread function from Fourier optics theory is equal to the
% Fourier transform of the Pupil function, which in this case is
% represented by apt. It has 1's within the radius 1 and 0 outside.
% Moreover, in real life we obsorve not the amplitude, but the intensity of
% waves, hence the modulus squared (=fn * conj(fn)) is done.
% The PSF is also normalized.
%
% The Fourier transform of the PSF gives the optical transfer function
% (OTF). The modulus of it is the MTF and the MTF goes from the centre
% towards the edge along the Y values

% Find the PSF
psf = fftshift(fft2(pupil)) .* conj(fftshift(fft2(pupil)));
% Normalize
psf = psf ./ sum(psf(:));
% Crop for display
centre = size(psf, 1) / 2 + 0.5
display = psf(centre-centre / 2 : centre + centre / 2 ,...
          centre-centre / 2 : centre + centre / 2 );

% Plot the PSF
figure; imshow(display, []); 
        title(strcat(imgtit, ' PSF'));
%         img = imadjust(display, [min(display(:)), max(display(:))]);
%         imwrite(img, strcat(imgname,'PSF.jpg'), 'jpg')

% Find the MTF
MTF = abs(fftshift(fft2(psf)));
centre = size(MTF, 1) / 2 + 0.5

% Find the OTF
PSF = fftshift(fft2(pupil)) .* conj(fftshift(fft2(pupil)));
OTF = fftshift(fft2(PSF));
PTF = angle(OTF);

figure; imshow(PTF, []);
figure; plot(unwrap(PTF(size(PTF,1)/2+1,:)))
        set(gca,'position',[0 0 1 1],'units','normalized');
        img = getframe(gcf);
        imwrite(img.cdata, strcat(imgname,'PTF.jpg'), 'jpg');
        
% % Get PTF along the horizontal direction
% x = linspace(0, size(PTF, 1), 1024);
% y = round(size(PTF, 1) / 2);
% interpPTF = unwrap(interp2(PTF, x, y, 'spline'));
% figure; plot(interpPTF)
        

% For a small aperture, the PSF is big and MTF is small,
% therefore, the MTF is croped a bit around the edges for
% better visibility
% MTF = MTF(centre-centre / 2 : centre + centre / 2 ,...
%           centre-centre / 2 : centre + centre / 2 );
%       
% centre and length values for the MTF matrix
centre = size(MTF, 1) / 2 + 0.5
length = size(MTF, 1)


% The index arrays used to access horizontal and vertical
% values for the MTF matrix

radius = length - centre;

% Interpolate values along the lines at different angles
[x1, y1] = getlines(centre, radius, 0);
interpMTF1 = interp2(MTF, x1, y1);

[x2, y2] = getlines(centre, radius, pi/4);
interpMTF2 = interp2(MTF, x2, y2);

[x3, y3] = getlines(centre, radius, pi/2);
interpMTF3 = interp2(MTF, x3, y3);

[x4, y4] = getlines(centre, radius, pi/3);
interpMTF4 = interp2(MTF, x4, y4);

% Plot MTF along horizontal, vertical and 45 degree directions
figure; plot(interpMTF1, 'r'); 
        hold on; 
        plot(interpMTF2, 'g');
        hold on;
        plot(interpMTF3, 'b');
        hold on;
        plot(interpMTF4, 'y');
        
        %title(strcat(imgtit, ' MTF'))
%           ylim([0, 0.3]);
%           xlim([0, 255]);
         saveas(gcf, strcat(imgname, 'MTF.jpg'));        
        
        
% Plot 2D MTF
figure; imshow(MTF, []);
        hold on
        plot(x1, y1, 'r');
        hold on
        plot(x2, y2, 'g');
        hold on
        plot(x3, y3, 'b');
        hold on
        plot(x4, y4, 'y');
        title(strcat(imgtit, ' MTF in 2D Space'));
        
        set(gca,'position',[0 0 1 1],'units','normalized')
        img = getframe(gcf);
        imwrite(img.cdata, [strcat(imgname, 'MTFspace'), '.jpg']);
        %saveas(gca, strcat(imgname, 'MTFspace.jpg'));

end