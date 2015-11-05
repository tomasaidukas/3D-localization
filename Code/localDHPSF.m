function [img2D, depths] = localDHPSF(img, original, anglemap, depthmap)
%------------------------------------------------------------%
% Unused function.
%------------------------------------------------------------%
    % Filter noise. These are usually just corrupted pixels 1 pixel wide,
    % while the PSF peaks of interest are much larger than that. Median
    % filter removes the single pixels.
    [Gmag,Gdir] = imgradient(img);
    figure; imshow(Gmag.^2, [])
    
    imgF = medfilt2(img, [3,3]);
    imgF = medfilt2(imgF, [5,5]);

    % Threshold
    thresh = multithresh(imgF);
    binary = (imquantize(imgF, thresh) - 1);        
    figure; imshow(binary, [])
    
    % Merge noisy pixels and open them to be left only with the peaks
    diskelem = strel('disk', 2);
    closed = imclose(binary, diskelem);
    diskelem = strel('disk', 9);
    opened = imopen(closed, diskelem);  
    figure; imshow(opened, [])
    
    [L, num] = bwlabel(opened);
    rp = regionprops(L, 'Centroid');
    
    % Extract the DH-PSF's from the filtered image
    % and locate the peaks using circular hough transform
    L = 30;
    peaks = zeros(size(img, 1), size(img, 2));
    peakCoords = [];
    for i = 1 : num
        
        % Absolute co-ordinate
        C = round(rp(i).Centroid);
        
        % Extract the region and filter it with a gaussian
        box = imgF(C(2)-L:C(2)+L, C(1)-L:C(1)+L);
        sz = size(box, 1);

%         figure; imshowpair(boxpeaks, box)
%         % Threshold and find peaks
%         thresh = multithresh(region);
%         BWregion = (imquantize(region, thresh) - 1);
%         [centers, radii] = imfindcircles(region, [1 10]);
%         
%         figure; imshow(region, [])
%         
%         % If two peaks were found
%         if numel(radii) == 2
%             % Peak co-ordinates
%             peak1 = centers(1, :)
%             peak2 = centers(2, :)
%             
%             % Shift from the centroid for each peak
%             shift1 = round([C(2) + peak1(2) - sz/2, C(1) + peak1(1) - sz/2]);
%             shift2 = round([C(2) + peak2(2) - sz/2, C(1) + peak2(1) - sz/2]);
%             
%             figure; imshow(region, [])
%             viscircles(centers, radii,'EdgeColor','b');      
%         end
        
    end
    
    % Hough circle detection
    [centers, radii] = imfindcircles(opened, [1 20]);
    
    viscircles(centers, radii,'EdgeColor','b');
    
    % Store the centroids on the image
    peaks = zeros(size(img, 1), size(img, 2));
    localized = zeros(size(original, 1), size(original, 2));
    peakCoords = [];
    angles = [];
    midpoints = [];
    
    L = 10;
    % Find the centroids using least squares
    for i = 1 : size(centers, 1)
        
        boxpeaks = zeros(L*2 + 1, L*2 + 1);
        c = centers(i, :);
        box = imgF(round(c(2))-L:round(c(2))+L,...
                   round(c(1))-L:round(c(1))+L);
        boxC = size(box, 1);
        
        % Box size
        [n, m] = size(box);
        % X, Y co-ordinates
        [X, Y] = meshgrid(1:n, 1:m);
        % Tolerance
        options = optimset('TolX', 1e-10);
        % guess [normalization, xc, yc, sigx, sigy]
        guess = [1, n/2, n/2, 1, 1];
        LB = [-inf, 1, 1, -inf, -inf];
        UB = [inf, n, n, inf, inf];
        % least-squares fitting for a 2D Gaussian
        params = lsqnonlin(@(P) objfun(P, X, Y, box), guess, LB, UB, options);
        % get the Gaussian fits       
        fitted = gauss2D(params, X, Y);
        fitted = reshape(fitted, [n, m]);
        
        % Shift the absolute co-ordinates
        coords = [c(1) - boxC / 2 + params(2), ...
                  c(2) - boxC / 2 + params(3)];
        
        peakCoords = [peakCoords; coords];  
        peaks(round(coords(2)), round(coords(1))) = 1;        
        boxpeaks(round(params(3)), round(params(2))) = 1;
%         figure; surf(X, Y, fitted)    
%         figure; imshowpair(boxpeaks, box)
    end    

    % Find peaks that correspond to the same DH-PSF
    % I.e. find two points that have the smallest distance between them
    % This requires the points to be adequately spaced
    numPts = size(peakCoords, 1)

    % each row of peakCoords stores the x,y coordinates
    % Find the minimum distance between to points
%     figure; imshow(imgF, [])
    for i = 1 : numPts - 1
        bestDist = 10000; 
        for j = i + 1 : numPts
            peak1 = [peakCoords(i, 1), peakCoords(i, 2)];
            peak2 = [peakCoords(j, 1), peakCoords(j, 2)];
            
            dist = sqrt((peak1(1) - peak2(1)) ^ 2 + ...
                        (peak1(2) - peak2(2)) ^ 2 );
                    
            % smallest dist indicates the closest pair
            if dist < bestDist
                bestDist = dist;
                Pair = [peak1, peak2];
            end
        end
        
        % Check if the pairs are correct and measure
        % the angle between them in order to estimate depth
%         hold on
%         plot([Pair(1), Pair(3)], ...
%              [Pair(2), Pair(4)], 'r')       
        peakCoords(peakCoords == Pair(3)) = inf;
        peakCoords(peakCoords == Pair(4)) = inf;
        
        % Measure angles w.r.t. the horizontal axis
        % dot(v1,v2) = mod(v1)mod(v2)cos(angle)
        X = Pair(3) - Pair(1);
        Y = Pair(4) - Pair(2);
        
        % Normalized vector between the lobes
        pairVect = [X, Y];
        pairLen = sqrt( X^2 + Y^2 );
        pairVect = pairVect ./ pairLen;
        
        % Unit vector in horizontal direction
        unitVect = [1, 0];        
        
        % Find the angle
        angle = acosd(abs(dot(pairVect, unitVect)));
        % Find the midpoint between peaks
        midpt = [Pair(1) + Pair(3), Pair(2) + Pair(4)] ./ 2;
        
        midpoints = [midpoints ; midpt];
        angles = [angles ; angle];
    end
    
    % Unique angles
    uniqMidpts = unique(midpoints, 'rows');
    uniqAngles = unique(angles, 'rows')
    
    % Get the depths
    depths = [];
    
    for i = 1 : size(uniqAngles, 1)
        angle = uniqAngles(i);
        midpt = round([uniqMidpts(i, 1), uniqMidpts(i, 2)]);
        depth = angle2depth(anglemap, depthmap, angle);
        depths = [depths ; depthmap(depth)];
        
        
        localized(midpt(2), midpt(1)) = 1;
    end
    
    img2D = localized;
    
%     figure; imshow(img, []);
%     figure; imshow(imgF, []);
%     figure; imshow(binary, []);
    
%     figure; imshow(imgF, []);
%     hold on
%     viscircles(centers, radii, 'EdgeColor', 'b');
%     
%     figure; imshow(peaks, []);
%     figure; imshow(localized, []);
%     figure; imshow(original - localized, []);
     
%     figure; imshowpair(original, peaks)
%     figure; imshowpair(original, houghpeaks)
%     figure; imshowpair(img, peaks)    
%     figure; imshowpair(peaks, houghpeaks)

end