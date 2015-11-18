clear all

%------------------------------------------------------------%
% Generate a single image at the centre.
%------------------------------------------------------------% 
SNRrange = 10:30;
SampleRange = 1:50;
DefocusRange = 0:2:2;

for W20 = DefocusRange
    % Camera and mesh grid parameters
    maxDefocus = size(W20, 2);
    % Data used to simulate the PSF mesh
    NoPts = 870;    XYrange = 2;    R = 1;    Location = NoPts / 2;
    f = 0.1;    NOISE = 0.0001;     sigmaRange = 10:15;

    camera = {W20; maxDefocus; NoPts; XYrange; R; f};
    
    %------------------------------------------------------------%
    % DH psf's
    %------------------------------------------------------------%
    deconpsfDH = simulateDHPSF(NoPts, XYrange, R, 0);
    blurpsfDH = simulateDHPSF(NoPts, XYrange, R, W20);
    
    %------------------------------------------------------------%
    % CKM psf's
    %------------------------------------------------------------%
    [deconpsf, deconpsf180] = CPMpsf(XYrange, NoPts, R, 0);
    [blurpsf, blurpsfConj] = CPMpsf(XYrange, NoPts, R, W20); 
    
    %------------------------------------------------------------%
    % Calibration maps
    %------------------------------------------------------------%
    CKMfit = load('../CKM/CKMdata/CKMfitLinNOSQRT2.mat'); CKMfitCKM = CKMfit.CKMfit;
    W20fit = load('./DHPSFdata/ang2defocus2.mat'); angle2defocus = W20fit.W20fit;

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

    %------------------------------------------------------------%
    % DHPSF corrections
    %------------------------------------------------------------%
    [X, Y, Z] = analysisDHPSF(imgDH, deconpsfDH, angle2defocus);
    correct = [X - randx, Y - randy];
%     correct = [0, 0];

    % Store the final results to these arrays
    ckmVSdepth = []; ckmVSX = []; ckmVSY = [];
    dhpsfVSdepth = []; dhpsfVSX = []; dhpsfVSY = [];

    errorX = []; errorY = []; errorZ = [];
    errorXDH = []; errorYDH = []; errorZDH = [];

    SNR = [];

    for PSNR = SNRrange

        % Store the intermediate results
        ckmX = []; ckmY = []; ckmZ = [];
        XDH = []; YDH = []; ZDH = [];

        %------------------------------------------------------------%
        % Noise and SNR
        %------------------------------------------------------------%
        % PSNR = 10*log(Imax^2 / MSE)
        % MSE = variance for unbiased noise (mean = 0)
    %     imgPSNR = 10 * log10(Imax^2 / var);
    %     SNR = 20*log10(mean(img(:)) / sqrt(var))
        var = Imax^2 / 10^(PSNR / 10);

        
        for i = SampleRange

            %------------------------------------------------------------%
            % Add noise
            %------------------------------------------------------------%
            imgN = imnoise(img, 'gaussian', 0, var);
            img180N = imnoise(img180, 'gaussian', 0, var);
            imgDHN = imnoise(imgDH, 'gaussian', 0, var);
            
%             figure; imshow(imgN, [])
%             figure; imshow(imgDHN, [])
            %------------------------------------------------------------%
            % CKM method.
            %------------------------------------------------------------%
            try
                [depth, coords] = analysisCKMdeconv(imgN, img180N, ...
                                     deconpsf, deconpsf180, camera, ...
                                     CKMfitCKM, NOISE, sigmaRange);
                %------------------------------------------------------------%
                % Localized co-ordinates for CKM (add 0.5 to remove the
                % constant error due to some poor calibration)
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
                                          angle2defocus);
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

        SNR = [SNR, PSNR];
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




% % clear all
% 
% ckmVSX = open('analysis/ckmVSX2.mat'); ckmVSX = ckmVSX.ckmVSX;
% ckmVSY = open('analysis/ckmVSY2.mat'); ckmVSY = ckmVSY.ckmVSY;
% ckmVSdepth = open('analysis/ckmVSdepth2.mat'); ckmVSdepth = ckmVSdepth.ckmVSdepth;
% 
% dhpsfVSX = open('analysis/dhpsfVSX2.mat'); dhpsfVSX = dhpsfVSX.dhpsfVSX;
% dhpsfVSY = open('analysis/dhpsfVSY2.mat'); dhpsfVSY = dhpsfVSY.dhpsfVSY;
% dhpsfVSdepth = open('analysis/dhpsfVSdepth2.mat'); dhpsfVSdepth = dhpsfVSdepth.dhpsfVSdepth;
% 
% errorX = open('analysis/errorX2.mat'); errorX = errorX.errorX;
% errorY = open('analysis/errorY2.mat'); errorY = errorY.errorY;
% errorZ = open('analysis/errorZ2.mat'); errorZ = errorZ.errorZ;
% 
% errorXDH = open('analysis/errorXDH2.mat'); errorXDH = errorXDH.errorXDH;
% errorYDH = open('analysis/errorYDH2.mat'); errorYDH = errorYDH.errorYDH;
% errorZDH = open('analysis/errorZDH2.mat'); errorZDH = errorZDH.errorZDH;
% 
% 
% 
% SNR = 10:30;    W20 = 2;   maxDefocus = size(W20, 2);
% NoPts = 870;    XYrange = 2;      R = 1; % 3 mm  
% Location = NoPts / 2;   f = 0.1;    NOISE = 0.0006;
% 
% 
% camera = {W20; maxDefocus; NoPts; XYrange; R; f};
% sigmaRange = 5:10;  randx = 400; randy = 400;


%     figure; errorbar( (SNR), W20 - ckmVSdepth, errorZ, 'b*'), hold on
%             errorbar( (SNR), W20 - dhpsfVSdepth, errorZDH, 'r*'), hold on
figure;     plot( (SNR), W20 - ckmVSdepth, 'b'), hold on
            plot( (SNR), W20 - dhpsfVSdepth, 'r')
            set(gca, 'xdir','reverse')
            title('Depth error vs PSNR'); xlabel('PSNR'); ylabel('Depth error')
            legend('CKM','DH-PSF')
            saveas(gcf, strcat(strcat('analysis/X' , num2str(W20)), '.png'))
            savefig(strcat(strcat('analysis/X' , num2str(W20))))

%     figure; errorbar( (SNR), randx - ckmVSX, errorX, 'b*'), hold on
%             errorbar( (SNR), randx - dhpsfVSX, errorXDH, 'r*'), hold on
figure;     plot( (SNR), randx - ckmVSX, 'b'), hold on
            plot( (SNR), randx - dhpsfVSX, 'r')
            set(gca, 'xdir','reverse')
            title('X error vs PSNR'); xlabel('PSNR'); ylabel('X error')
            legend('CKM','DH-PSF')
            saveas(gcf, strcat(strcat('analysis/Y' , num2str(W20)), '.png'))
            savefig(strcat(strcat('analysis/Y' , num2str(W20))))

%     figure; errorbar( (SNR), randy - ckmVSY, errorY, 'b*'), hold on
%             errorbar( (SNR), randy - dhpsfVSY, errorYDH, 'r*'), hold on
figure;     plot( (SNR), randy - ckmVSY, 'b'), hold on
            plot( (SNR), randy - dhpsfVSY, 'r')
            set(gca, 'xdir','reverse')
            title('Y error vs PSNR'); xlabel('PSNR'); ylabel('Y error')
            legend('CKM','DH-PSF')
            saveas(gcf, strcat(strcat('analysis/Z' , num2str(W20)), '.png'))
            savefig(strcat(strcat('analysis/Z' , num2str(W20))))

    
    
figure;     plot( (SNR), errorZ, 'b'), hold on
            plot( (SNR), errorZDH, 'r')
            set(gca, 'xdir','reverse')
            title('Depth stdev vs PSNR'); xlabel('PSNR'); ylabel('Depth stdev')
            legend('CKM','DH-PSF')
            saveas(gcf, strcat(strcat('analysis/Z' , num2str(W20)), '.png'))
            savefig(strcat(strcat('analysis/Z' , num2str(W20))))


figure;     plot( (SNR), errorX, 'b'), hold on
            plot( (SNR), errorXDH, 'r')
            set(gca, 'xdir','reverse')
            title('X stedv vs PSNR'); xlabel('PSNR'); ylabel('X stdev')
            legend('CKM','DH-PSF')
            saveas(gcf, strcat(strcat('analysis/X' , num2str(W20)), '.png'))
            savefig(strcat(strcat('analysis/X' , num2str(W20))))


figure;     plot( (SNR), errorY, 'b'), hold on
            plot( (SNR), errorYDH, 'r')
            set(gca, 'xdir','reverse')
            title('Y stdev vs PSNR'); xlabel('PSNR'); ylabel('Y stdev')
            legend('CKM','DH-PSF')
            saveas(gcf, strcat(strcat('analysis/Y' , num2str(W20)), '.png'))
            savefig(strcat(strcat('analysis/Y' , num2str(W20))))

        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% 
% %------------------------------------------------------------%
% % Generate a single image at the centre.
% %------------------------------------------------------------% 
% SNRrange = 10:30;
% SampleRange = 1:50;
% DefocusRange = 0:2;
% 
% for W20 = DefocusRange
%     % Camera and mesh grid parameters
%     maxDefocus = size(W20, 2);
%     % Data used to simulate the PSF mesh
%     NoPts = 870;    XYrange = 2;    R = 1;    Location = NoPts / 2;
%     f = 0.1;    NOISE = 0.0006;     sigmaRange = 10:15;
% 
%     camera = {W20; maxDefocus; NoPts; XYrange; R; f};
%     
%     
%     %------------------------------------------------------------%
%     % DH psf's
%     %------------------------------------------------------------%
%     deconpsfDH = simulateDHPSF(NoPts, XYrange, R, 0);
%     blurpsfDH = simulateDHPSF(NoPts, XYrange, R, W20);
%     
%     %------------------------------------------------------------%
%     % CKM psf's
%     %------------------------------------------------------------%
%     [deconpsf, deconpsf180] = CPMpsf(XYrange, NoPts, R, 0);
%     [blurpsf, blurpsfConj] = CPMpsf(XYrange, NoPts, R, W20); 
%     
%     %------------------------------------------------------------%
%     % Calibration maps
%     %------------------------------------------------------------%
%     CKMfit = load('../CKM/CKMdata/CKMfitLinNOSQRT2.mat'); CKMfitCKM = CKMfit.CKMfit;
%     W20fit = load('./DHPSFdata/ang2defocus2.mat'); angle2defocus = W20fit.W20fit;
% 
%     %------------------------------------------------------------%
%     % Create the convolved images.
%     %------------------------------------------------------------%
%     canvas = zeros(NoPts, NoPts);
%     % canvas(Location, Location) = 1;
% %     randx = randi([50, NoPts - 50], 1, 1);
% %     randy = randi([50, NoPts - 50], 1, 1);
%     randx = 400; randy = 400;
%     canvas(randx, randy) = 1;
% 
%     img = abs(fftshift(ifft2(fftshift(fft2(blurpsf)) .* fftshift(fft2(canvas)))));
%     img180 = abs(fftshift(ifft2(fftshift(fft2(blurpsfConj)) .* fftshift(fft2(canvas)))));
%     imgDH = abs(fftshift(ifft2(fftshift(fft2(blurpsfDH)) .* fftshift(fft2(canvas)))));
% 
%     Imax = max(img(:));
%     DHmax = max(imgDH(:));
% 
%     %------------------------------------------------------------%
%     % DHPSF corrections
%     %------------------------------------------------------------%
%     [X, Y, Z] = analysisDHPSF(imgDH, deconpsfDH, angle2defocus);
%     correct = [X - randx, Y - randy];
% %     correct = [0, 0];
% 
%     % Store the final results to these arrays
%     ckmVSdepth = []; ckmVSX = []; ckmVSY = [];
%     dhpsfVSdepth = []; dhpsfVSX = []; dhpsfVSY = [];
% 
%     errorX = []; errorY = []; errorZ = [];
%     errorXDH = []; errorYDH = []; errorZDH = [];
% 
%     SNR = [];
% 
%     for PSNR = SNRrange
% 
%         % Store the intermediate results
%         ckmX = []; ckmY = []; ckmZ = [];
%         XDH = []; YDH = []; ZDH = [];
% 
%         %------------------------------------------------------------%
%         % Noise and SNR
%         %------------------------------------------------------------%
%         % PSNR = 10*log(Imax^2 / MSE)
%         % MSE = variance for unbiased noise (mean = 0)
%     %     imgPSNR = 10 * log10(Imax^2 / var);
%         var = Imax^2 / 10^(PSNR / 10);
% 
%         for i = SampleRange
% 
%             %------------------------------------------------------------%
%             % Add noise
%             %------------------------------------------------------------%
%             imgN = imnoise(img, 'gaussian', 0, var);
%             img180N = imnoise(img180, 'gaussian', 0, var);
%             imgDHN = imnoise(imgDH, 'gaussian', 0, var);
%             
% %             figure; imshow(imgN, [])
% %             figure; imshow(imgDHN, [])
%             %------------------------------------------------------------%
%             % CKM method.
%             %------------------------------------------------------------%
%             try
%                 [depth, coords] = analysisCKMdeconv(imgN, img180N, ...
%                                      deconpsf, deconpsf180, camera, ...
%                                      CKMfitCKM, NOISE, sigmaRange);
%                 %------------------------------------------------------------%
%                 % Localized co-ordinates for CKM (add 0.5 to remove the
%                 % constant error due to some poor calibration)
%                 %------------------------------------------------------------%
%                 ckmX = [ckmX, coords(2) + 0.5];
%                 ckmY = [ckmY, coords(1) + 0.5];
%                 ckmZ = [ckmZ, depth];
%             catch
%                 ckmX = [ckmX, 0];
%                 ckmY = [ckmY, 0];
%                 ckmZ = [ckmZ, 0];
%             end
%             %------------------------------------------------------------%
%             % Do the same for DHPSF
%             %------------------------------------------------------------%
%             try
%                 [X, Y, Z] = analysisDHPSF(imgDHN, deconpsfDH, ...
%                                           angle2defocus);
%                 XDH = [XDH, X - correct(1)];
%                 YDH = [YDH, Y - correct(2)];
%                 ZDH = [ZDH, Z];
%             catch
%                 XDH = [XDH, 0];
%                 YDH = [YDH, 0];
%                 ZDH = [ZDH, 0];
%             end
% 
% 
%         end
% 
%         % Means
%         temp1 = mean(ckmX); temp2 = mean(ckmY); temp3 = mean(ckmZ);
%         temp4 = mean(XDH); temp5 = mean(YDH); temp6 = mean(ZDH);
% 
%         % Errors as standard deviations
%         error1 = std(ckmX); error2 = std(ckmY); error3 = std(ckmZ);
%         error4 = std(XDH); error5 = std(YDH); error6 = std(ZDH);
% 
%         % Means for each data point
%         ckmVSX = [ckmVSX, temp1];
%         dhpsfVSX = [dhpsfVSX, temp4];
% 
%         ckmVSY = [ckmVSY, temp2];
%         dhpsfVSY = [dhpsfVSY, temp5];
% 
%         ckmVSdepth = [ckmVSdepth, temp3];
%         dhpsfVSdepth = [dhpsfVSdepth, temp6];
% 
%         % Errors for each data point
%         errorX = [errorX, error1];
%         errorXDH = [errorXDH, error4];
% 
%         errorY = [errorY, error2];
%         errorYDH = [errorYDH, error5];
% 
%         errorZ = [errorZ, error3];
%         errorZDH = [errorZDH, error6];
% 
%         SNR = [SNR, PSNR];
%     end
% 
% 
%     save(strcat(strcat('analysis2/ckmVSX' , num2str(W20)), '.mat'), 'ckmVSX')
%     save(strcat(strcat('analysis2/ckmVSY' , num2str(W20)), '.mat'), 'ckmVSY')
%     save(strcat(strcat('analysis2/ckmVSdepth' , num2str(W20)), '.mat'), 'ckmVSdepth')
% 
%     save(strcat(strcat('analysis2/dhpsfVSX' , num2str(W20)), '.mat'), 'dhpsfVSX')
%     save(strcat(strcat('analysis2/dhpsfVSY' , num2str(W20)), '.mat'), 'dhpsfVSY')
%     save(strcat(strcat('analysis2/dhpsfVSdepth' , num2str(W20)), '.mat'), 'dhpsfVSdepth')
% 
%     save(strcat(strcat('analysis2/errorX', num2str(W20)), '.mat'), 'errorX')
%     save(strcat(strcat('analysis2/errorY', num2str(W20)), '.mat'), 'errorY')
%     save(strcat(strcat('analysis2/errorZ', num2str(W20)), '.mat'), 'errorZ')
% 
%     save(strcat(strcat('analysis2/errorXDH', num2str(W20)), '.mat'), 'errorXDH')
%     save(strcat(strcat('analysis2/errorYDH', num2str(W20)), '.mat'), 'errorYDH')
%     save(strcat(strcat('analysis2/errorZDH', num2str(W20)), '.mat'), 'errorZDH')
% 
% 
% 
% 
% 
% % clear all
% % 
% % ckmVSX = open('analysis/ckmVSX3.mat'); ckmVSX = ckmVSX.ckmVSX;
% % ckmVSY = open('analysis/ckmVSY3.mat'); ckmVSY = ckmVSY.ckmVSY;
% % ckmVSdepth = open('analysis/ckmVSdepth3.mat'); ckmVSdepth = ckmVSdepth.ckmVSdepth;
% % 
% % dhpsfVSX = open('analysis/dhpsfVSX3.mat'); dhpsfVSX = dhpsfVSX.dhpsfVSX;
% % dhpsfVSY = open('analysis/dhpsfVSY3.mat'); dhpsfVSY = dhpsfVSY.dhpsfVSY;
% % dhpsfVSdepth = open('analysis/dhpsfVSdepth3.mat'); dhpsfVSdepth = dhpsfVSdepth.dhpsfVSdepth;
% % 
% % errorX = open('analysis/errorX3.mat'); errorX = errorX.errorX;
% % errorY = open('analysis/errorY3.mat'); errorY = errorY.errorY;
% % errorZ = open('analysis/errorZ3.mat'); errorZ = errorZ.errorZ;
% % 
% % errorXDH = open('analysis/errorXDH3.mat'); errorXDH = errorXDH.errorXDH;
% % errorYDH = open('analysis/errorYDH3.mat'); errorYDH = errorYDH.errorYDH;
% % errorZDH = open('analysis/errorZDH3.mat'); errorZDH = errorZDH.errorZDH;
% % 
% % SNR = 10:30;    W20 = 2;   maxDefocus = size(W20, 2);
% % NoPts = 870;    XYrange = 2;      R = 1; % 3 mm  
% % Location = NoPts / 2;   f = 0.1;    NOISE = 0.0006;
% % 
% % camera = {W20; maxDefocus; NoPts; XYrange; R; f};
% % sigmaRange = 5:10;  randx = 400; randy = 400;      
% 
% 
% 
%     figure; errorbar( (SNR), W20 - ckmVSdepth, errorZ, 'b*'), hold on
%             errorbar( (SNR), W20 - dhpsfVSdepth, errorZDH, 'r*'), hold on
%             plot( (SNR), W20 - ckmVSdepth, 'b'), hold on
%             plot( (SNR), W20 - dhpsfVSdepth, 'r')
%             set(gca, 'xdir','reverse')
%             title('Depth error vs PSNR'); xlabel('PSNR'); ylabel('Depth error')
%             legend('CKM','DH-PSF')
%     saveas(gcf, strcat(strcat('analysis2/X' , num2str(W20)), '.png'))
%     savefig(strcat(strcat('analysis2/X' , num2str(W20))))
% 
% %     figure; errorbar( (SNR), randx - ckmVSX, errorX, 'b*'), hold on
% %             errorbar( (SNR), randx - dhpsfVSX, errorXDH, 'r*'), hold on
% figure;     plot( (SNR), randx - ckmVSX, 'b'), hold on
%             plot( (SNR), randx - dhpsfVSX, 'r')
%             set(gca, 'xdir','reverse')
%             title('X error vs PSNR'); xlabel('PSNR'); ylabel('X error')
%             legend('CKM','DH-PSF')
%     saveas(gcf, strcat(strcat('analysis2/Y' , num2str(W20)), '.png'))
%     savefig(strcat(strcat('analysis2/Y' , num2str(W20))))
% 
%     figure; errorbar( (SNR), randy - ckmVSY, errorY, 'b*'), hold on
%             errorbar( (SNR), randy - dhpsfVSY, errorYDH, 'r*'), hold on
%             plot( (SNR), randy - ckmVSY, 'b'), hold on
%             plot( (SNR), randy - dhpsfVSY, 'r')
%             set(gca, 'xdir','reverse')
%             title('Y error vs PSNR'); xlabel('PSNR'); ylabel('Y error')
%             legend('CKM','DH-PSF')
%     saveas(gcf, strcat(strcat('analysis2/Z' , num2str(W20)), '.png'))
%     savefig(strcat(strcat('analysis2/Z' , num2str(W20))))
% 
% end