function [img3D, img2D, depths, plotsX, plotsY] = localCKM(map, recon, W20)
%------------------------------------------------------------%
% Takes in a defocus map and a reconstructed image.
% Based on the defocus map, the pixels in the reconstructed image
% are associated a depth value.
%------------------------------------------------------------%

NoPts = size(recon, 1);
img3D = zeros(NoPts, NoPts, size(W20, 1));
img2D = zeros(NoPts, NoPts);    
depths = []; plotsX = []; plotsY = [];

%------------------------------------------------------------%
% Use dilation and thresholding to get the blobs containing the
% peaks. Use centroid of each blob to segment a region around
% the true peak.
%------------------------------------------------------------%
diskelem = strel('disk', 10);
closed = imdilate(recon, diskelem);

thresh = multithresh(closed) / 2;
BW = (imquantize(closed, thresh) - 1);

[L, num] = bwlabel(BW);
rp = regionprops(L, 'Centroid');
figure; imshowpair(BW, closed)

L = 15;
%------------------------------------------------------------%
% Find the exact location of the maxima for each reconstructed
% peak. This maxima is the 2D localisation co-ordinate.
%------------------------------------------------------------%
for i = 1 : num

    C = round(rp(i).Centroid);

    % Extract the region
    region = recon(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
    % Maxima in the region will be the position of the peak
    peakInt = max(region(:));
    [X, Y] = find(region == peakInt);

    sz = size(region, 1); CC = [X, Y];
    coords = round([C(2) - CC(2) + sz/2, C(1) - CC(1) + sz/2]);

    % Intensity
    I = recon(coords(1), coords(2));
    img2D(coords(1), coords(2)) = I;
    plotsX = [plotsX, coords(2)];
    plotsY = [plotsY, coords(1)];

    % Get the depth
    depth = map(coords(1), coords(2));
    depths = [depths; W20(depth)];
end

depths = sort(depths)
figure; imshow(map, [])
        hold on
        plot(plotsX, plotsY, '*')
end
