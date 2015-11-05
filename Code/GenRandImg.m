% Create random image
clear
close all

NoPts = 870;
sz = 0;
XYrange = 0.05;
R = 0.02;
W20range = 0:1:2;

% array to store all the images
finalimg = 0;
finalimg180 = 0;
original = 0;
finalimgDH = 0;

for W20 = W20range
    [img, imgC, imgO, imgDH] = randompts(XYrange, NoPts, R, sz, W20, 1);
    finalimg = finalimg + img;
    finalimg180 = finalimg180 + imgC;
    finalimgDH = finalimgDH + imgDH;
    original = original + imgO;
end

var = 0.0000001;

finalimg = imnoise(finalimg, 'gaussian', 0, var);
finalimg180 = imnoise(finalimg180, 'gaussian', 0, var);
finalimgDH = imnoise(finalimgDH, 'gaussian', 0, var);
finalimg = imnoise(finalimg, 'poisson');
finalimg180 = imnoise(finalimg180, 'poisson');
finalimgDH = imnoise(finalimgDH, 'poisson');

SNR(finalimg, var), SNR(finalimgDH, var), SNR(finalimg180, var)

figure; imshow(finalimg, [])
figure; imshow(finalimg180, [])
figure; imshow(finalimgDH, [])
figure; imshow(original, [])

save('randDH.mat', 'finalimgDH');
save('rand.mat', 'finalimg');
save('randC.mat', 'finalimg180');
save('randO.mat', 'original');
% save('pointsTHETA.mat', 'angles');