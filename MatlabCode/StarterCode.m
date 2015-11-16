clear

% Create a mesh, convert to polar. The created mesh represents an array of
% length 4 split into 512 pieces. Converting to polar allows to express it
% in angles and radius. Since we want a circular aperture of radius 1 only
% values <= are taken from the array of radii values.

[X, Y] = meshgrid(linspace(-2, 2, 2*512), linspace(-2, 2, 2*512));
[theta, p] = cart2pol(X, Y);
apt = double(p <= 1);

figure; imshow(apt, []);
        title('Pupil Function')
        imwrite(apt, 'pupil.jpg', 'jpg')
        
% Zernike polynomials representing aberrations (in the exponential form)
astigstr = sqrt(6);      % As seen from an online reference
defocusstr = sqrt(3);
comastr = sqrt(8);

defocus = exp(1i .* 2 .* pi .* defocusstr * (2 .* (p.^2) - 1));
astigmatism = exp(1i .* 2 .* pi .* astigstr * (p.^2) .* sin(2 .* theta));
coma = exp(1i .* 2 .* pi .* comastr * (3 * p.^3 - 2 * p) .* sin(theta));


% PUPIL FUNCTIONS WITH VARIOUS MASKS AND ABERRATIONS
pupil = apt;
pupilDefocus = apt .* defocus;
pupilAstig = apt .* astigmatism;
pupilComa = apt .* coma;

% CUBIC PHASE MASK
alpha = 5;      % Strength
CPM = exp(1i .* 2 .* pi .* alpha .* (X.^3 + Y.^3));

pupilCPM = apt .* CPM;
pupilDefocusCPM = apt .* CPM .* defocus;
pupilAstigCPM = apt .* CPM .* astigmatism;
pupilComaCPM = apt .* CPM .* coma;

% Find the PSF and MTF for all pupil functions
% makeplots(pupil, 'Diffraction Limited', 'normal')
% makeplots(pupilDefocus, 'Defocused', 'defocus')
% makeplots(pupilAstig, 'Astigmatic', 'astig')
% makeplots(pupilComa, 'Coma', 'coma')
% makeplots(pupilCPM, 'CPM', 'cpm')
% makeplots(pupilDefocusCPM, 'Defocused CPM', 'cpmdefocus')
% makeplots(pupilAstigCPM, 'Astigmatic CPM', 'cpmastig')
% makeplots(pupilComaCPM, 'Coma CPM', 'cpmcoma')

