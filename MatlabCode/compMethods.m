clear all
close all

%------------------------------------------------------------%
% Generate a single image at the centre.
%------------------------------------------------------------% 

for W20 = 3:0.2:4
    maxDefocus = size(W20, 2);
    NoPts = 870;
    XYrange = 0.05;
    R = 0.03;
    W20range = 1;
    Location = NoPts / 2;
    f = 0.1;
    NOISE = 0.0005;
    camera = {W20; maxDefocus; NoPts; XYrange; R; f};
    thresh = 0;
    sigmaRange = 5:10;


    [blurpsf, blurpsfConj] = CPMpsf(XYrange, NoPts, R, W20);   
    blurpsfDH = simulateDHPSF(NoPts, XYrange, R, W20);
    [deconpsf, deconpsf180] = CPMpsf(XYrange, NoPts, R, 0);
    deconpsfDH = simulateDHPSF(NoPts, XYrange, R, 0);

    CKMfit = load('../CKM/CKMdata/CKMfitLinNOSQRT.mat'); CKMfitCKM = CKMfit.CKMfit;
    shiftsX = load('./DHPSFdata/shiftX.mat'); shiftsX = shiftsX.shiftX;
    shiftsY = load('./DHPSFdata/shiftY.mat'); shiftsY = shiftsY.shiftY;
    W20fit = load('./DHPSFdata/ang2defocus.mat'); angle2defocus = W20fit.W20fit;

    %------------------------------------------------------------%
    % Create the convolved images.
    %------------------------------------------------------------%
    canvas = zeros(NoPts, NoPts);
    % canvas(Location, Location) = 1;
%     randx = randi([50, NoPts - 50], 1, 1);
%     randy = randi([50, NoPts - 50], 1, 1);
    randx = 400; randy = 400;
    canvas(randx, randy) = 1;

    img = abs(fftshift(ifft2(fftshift(fft2(blurpsf)) .* fftshift(fft2(canvas)))));
    img180 = abs(fftshift(ifft2(fftshift(fft2(blurpsfConj)) .* fftshift(fft2(canvas)))));
    imgDH = abs(fftshift(ifft2(fftshift(fft2(blurpsfDH)) .* fftshift(fft2(canvas)))));

    Imax = max(img(:));
    DHmax = max(imgDH(:));


    % DHPSF corrections
    [X, Y, Z] = analysisDHPSF(imgDH, deconpsfDH, ...
                              angle2defocus, shiftsX, shiftsY);
    correct = [X - randx, Y - randy];


    % Store the final results
    ckmVSdepth = []; ckmVSX = []; ckmVSY = [];
    dhpsfVSdepth = []; dhpsfVSX = []; dhpsfVSY = [];

    errorX = []; errorY = []; errorZ = [];
    errorXDH = []; errorYDH = []; errorZDH = [];

    SNR = [];

    for j = 1:30

    % Store the intermediate results
    ckmX = []; ckmY = []; ckmZ = [];
    XDH = []; YDH = []; ZDH = [];

    var = 0.0000000001 * 1.5^(j);
    % PSNR = 10*log(Imax^2 / MSE)
    % MSE = variance for unbiased noise (mean = 0)
    imgPSNR = 10 * log10(Imax^2 / var);

        for i = 1:1

            %------------------------------------------------------------%
            % Add noise
            %------------------------------------------------------------%
            imgN = imnoise(img, 'gaussian', 0, var);
            img180N = imnoise(img180, 'gaussian', 0, var);
            imgDHN = imnoise(imgDH, 'gaussian', 0, var);
            NOISE = 0.0005;

            %------------------------------------------------------------%
            % CKM method.
            %------------------------------------------------------------%
            try
                [defocus, breaker] = analysisCKMdeconv(imgN, img180N, deconpsf, deconpsf180, ...
                                     camera, CKMfitCKM, NOISE, randx, randy);
                camera{1} = defocus; % Defocus value for CKM
            
                [imgCKM, mapCKM] = CKM(imgN, img180N, camera, sigmaRange, NOISE);

                [ROW, COL] = find(max(imgCKM(:)) == imgCKM); L = 15;
                region = imgCKM(ROW-L:ROW+L, COL-L:COL+L);

                [n, m] = size(region); [X, Y] = meshgrid(1:n, 1:m);
                options = optimset('TolX', 1e-20, 'Display', 'off'); 
                peak = max(region(:)); [rc, cc] = find(peak == region);
                guess = [peak, rc(1), cc(1), 5, 5];
                LB = [0, 1, 1, 0, 0]; UB = [peak, n, n, 25, 25];

                % least square fit
                params = lsqnonlin(@(P) objfun2(P, X, Y, region), ...
                                   guess, LB, UB, options);
                coords = [COL + (n/2 - params(3)) ./ 1, ROW + (m/2 - params(2)) ./ 1];
                depth = mapCKM(round(coords(2)), round(coords(1)));
                defocusRange = camera{1};
                depth = defocusRange(depth);

                %------------------------------------------------------------%
                % Errors in x, y, z
                %------------------------------------------------------------%
                ckmX = [ckmX, coords(2) + 0.5];
                ckmY = [ckmY, coords(1) + 0.5];
                ckmZ = [ckmZ, depth];
            catch
                ckmX = [ckmX, 0];
                ckmY = [ckmY, 0];
                ckmZ = [ckmZ, 0];
            end
            %------------------------------------------------------------%
            % Do the same for DHPSF
            %------------------------------------------------------------%
            try
                [X, Y, Z] = analysisDHPSF(imgDHN, deconpsfDH, ...
                                          angle2defocus, shiftsX, shiftsY);
                XDH = [XDH, X - correct(1)];
                YDH = [YDH, Y - correct(2)];
                ZDH = [ZDH, Z];
            catch
                XDH = [XDH, 0];
                YDH = [YDH, 0];
                ZDH = [ZDH, 0];
            end


        end

        % Means
        temp1 = mean(ckmX); temp2 = mean(ckmY); temp3 = mean(ckmZ);
        temp4 = mean(XDH); temp5 = mean(YDH); temp6 = mean(ZDH);

        % Errors as standard deviations
        error1 = std(ckmX); error2 = std(ckmY); error3 = std(ckmZ);
        error4 = std(XDH); error5 = std(YDH); error6 = std(ZDH);

        % Means for each data point
        ckmVSX = [ckmVSX, temp1];
        dhpsfVSX = [dhpsfVSX, temp4];

        ckmVSY = [ckmVSY, temp2];
        dhpsfVSY = [dhpsfVSY, temp5];

        ckmVSdepth = [ckmVSdepth, temp3];
        dhpsfVSdepth = [dhpsfVSdepth, temp6];

        % Errors for each data point
        errorX = [errorX, error1];
        errorXDH = [errorXDH, error4];

        errorY = [errorY, error2];
        errorYDH = [errorYDH, error5];

        errorZ = [errorZ, error3];
        errorZDH = [errorZDH, error6];

        SNR = [SNR, imgPSNR];
    end


    save(strcat(strcat('analysis/ckmVSX' , num2str(W20)), '.mat'), 'ckmVSX')
    save(strcat(strcat('analysis/ckmVSY' , num2str(W20)), '.mat'), 'ckmVSY')
    save(strcat(strcat('analysis/ckmVSdepth' , num2str(W20)), '.mat'), 'ckmVSdepth')

    save(strcat(strcat('analysis/dhpsfVSX' , num2str(W20)), '.mat'), 'dhpsfVSX')
    save(strcat(strcat('analysis/dhpsfVSY' , num2str(W20)), '.mat'), 'dhpsfVSY')
    save(strcat(strcat('analysis/dhpsfVSdepth' , num2str(W20)), '.mat'), 'dhpsfVSdepth')

    save(strcat(strcat('analysis/errorX', num2str(W20)), '.mat'), 'errorX')
    save(strcat(strcat('analysis/errorY', num2str(W20)), '.mat'), 'errorY')
    save(strcat(strcat('analysis/errorZ', num2str(W20)), '.mat'), 'errorZ')

    save(strcat(strcat('analysis/errorXDH', num2str(W20)), '.mat'), 'errorXDH')
    save(strcat(strcat('analysis/errorYDH', num2str(W20)), '.mat'), 'errorYDH')
    save(strcat(strcat('analysis/errorZDH', num2str(W20)), '.mat'), 'errorZDH')

end



% clear all
% 
% ckmVSX = open('analysis/ckmVSX3.mat'); ckmVSX = ckmVSX.ckmVSX;
% ckmVSY = open('analysis/ckmVSY3.mat'); ckmVSY = ckmVSY.ckmVSY;
% ckmVSdepth = open('analysis/ckmVSdepth3.mat'); ckmVSdepth = ckmVSdepth.ckmVSdepth;
% 
% dhpsfVSX = open('analysis/dhpsfVSX3.mat'); dhpsfVSX = dhpsfVSX.dhpsfVSX;
% dhpsfVSY = open('analysis/dhpsfVSY3.mat'); dhpsfVSY = dhpsfVSY.dhpsfVSY;
% dhpsfVSdepth = open('analysis/dhpsfVSdepth3.mat'); dhpsfVSdepth = dhpsfVSdepth.dhpsfVSdepth;
% 
% errorX = open('analysis/errorX3.mat'); errorX = errorX.errorX;
% errorY = open('analysis/errorY3.mat'); errorY = errorY.errorY;
% errorZ = open('analysis/errorZ3.mat'); errorZ = errorZ.errorZ;
% 
% errorXDH = open('analysis/errorXDH3.mat'); errorXDH = errorXDH.errorXDH;
% errorYDH = open('analysis/errorYDH3.mat'); errorYDH = errorYDH.errorYDH;
% errorZDH = open('analysis/errorZDH3.mat'); errorZDH = errorZDH.errorZDH;
% 
% W20 = 3;
% maxDefocus = size(W20, 2);
% NoPts = 870;
% XYrange = 0.05;
% R = 0.03;
% W20range = 1;
% Location = NoPts / 2;
% f = 0.1;
% NOISE = 0.0005;
% camera = {W20; maxDefocus; NoPts; XYrange; R; f};
% thresh = 0;
% sigmaRange = 5:10;%:15;
% canvas = zeros(NoPts, NoPts);
% % randx = 811, randy = 291; % W20 = 3
% canvas(randx, randy) = 1;
% 
% [blurpsf, blurpsfConj] = CPMpsf(XYrange, NoPts, R, W20);
% 
% SNR = [];
% img = abs(fftshift(ifft2(fftshift(fft2(blurpsf)) .* fftshift(fft2(canvas)))));
% Imax = max(img(:));
% 
% for j = 1:30
% 
%     % Store the intermediate results
%     ckmX = []; ckmY = []; ckmZ = [];
%     XDH = []; YDH = []; ZDH = [];
%     var = 0.0000000001 * 1.5^(j);
%     imgPSNR = 10 * log10(Imax^2 / var);
%     SNR = [SNR, imgPSNR];
% end
% 
% figure; errorbar(fliplr(SNR), W20 - ckmVSdepth, errorZ, 'b*'), hold on
%         errorbar(fliplr(SNR), W20 - dhpsfVSdepth, errorZDH, 'r*')
%         title('Depth error vs PSNR'); xlabel('PSNR'); ylabel('Depth error')
%         legend('CKM','DH-PSF')
% 
% figure; errorbar(fliplr(SNR), randx - ckmVSX, errorX, 'b*'), hold on
%         errorbar(fliplr(SNR), randx - dhpsfVSX, errorXDH, 'r*')
%         title('X error vs PSNR'); xlabel('PSNR'); ylabel('X error')
%         legend('CKM','DH-PSF')
%         
% figure; errorbar(fliplr(SNR), randy - ckmVSY, errorY, 'b*')
%         hold on, errorbar(fliplr(SNR), randy - dhpsfVSY, errorYDH, 'r*')
%         title('Y error vs PSNR'); xlabel('PSNR'); ylabel('Y error')
%         legend('CKM','DH-PSF')
%         
%        