function img = cortrans(img, PSFs)
% Perform Phase correlation to locate the rough peak positions
    
    imgF = img;
    
    corrstack = 0;
    for i = 1 : size(PSFs, 3)
    
        psf = PSFs(:, :, i);
        % Low pass filter in the frequency domain
        h = fspecial('gaussian', [6, 6], 1.5);
        H = fftshift(fft2(h, size(img, 1), size(img, 2)));
        
        % Find the Fourier transform of both images
        FTpsf = fftshift(fft2(psf));
        FTimg = fftshift(fft2(imgF));
        divcc = abs(FTimg) .* abs(FTpsf);
        divac = abs(FTimg) .* abs(FTimg);
        divac(divac == 0) = 0.000000000000001;
        divcc(divcc == 0) = 0.000000000000001;
        
        % Power spectrum
        CC = double(FTimg .* conj(FTpsf) .* H) ./ divcc;
        % Cross-Correlation
        cc = abs(ifftshift(ifft2(CC)));
        
        corrstack = corrstack + cc;
    end
    
    % Threshold the image
    thresh = multithresh(corrstack);
    binary = (imquantize(corrstack, thresh) - 1);

    % Use morphological operations to obtain connected regions, which
    % contain the peaks.
    diskelem = strel('disk', 6);
    closed = imdilate(binary, diskelem);
    figure; imshow(closed)
    
    % Segment these peaks 
    [L, num] = bwlabel(closed);
    rp = regionprops(L, 'Centroid');
    K = 25;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each peak pair find their midpoints and angles to indicate the
    % depth by fitting a 2D double Gaussian
    for i = 1 : num
        prop = rp(i);
        C = prop.Centroid;
        
        % Extract the region
        box = imgF(C(2)-K:C(2)+K, C(1)-K:C(1)+K);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Least square fitting routine for a double Gaussian
        % Arrays used for result storage
        peakCoords = []; c = size(psf) / 2;
        boxC = size(box, 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Initial guess for the centre locations using Hough circle
        % detector
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
                figure; imshow(box, [])
                        viscircles([centers(1,:); centers(k,:)], ...
                                   [radii(1); radii(k)], 'EdgeColor','b');
                break
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Initial guess for the peak heights using the peak values at the
        % centroid positions
        
        [n, m] = size(box); [X, Y] = meshgrid(1:n, 1:m);
        options = optimset('TolX', 1e-20);

        % guess [normalization, xc, yc, sigma,
        %        normalization, xc, yc, sigma]
        guess = [max(box(:)), xc1, yc1, 1, ...
                 max(box(:)), xc2, yc2, 1];
        LB = [max(box(:))/4, 1, 1, -20, ...
              max(box(:))/4, 1, 1, -20];
        UB = [max(box(:)), n, n, 20, ...
              max(box(:)), n, n, 20];

        % least square fit
        params = lsqnonlin(@(P) objfun(P, X, Y, box), guess, LB, UB, options);

        % Use the best fit parameters to compute the Gaussian
        fitted = gauss2D(params, X, Y);
        fitted = reshape(fitted, [n, m]);

        % Shift the absolute co-ordinates
        coords1 = [C(1) - boxC / 2 + params(2), ...
                   C(2) - boxC / 2 + params(3)];
        coords2 = [C(1) - boxC / 2 + params(6), ...
                   C(2) - boxC / 2 + params(7)];        
        
        % Plot peaks superimposed inside the box for testing
        figure; imshow(box, [])
                hold on; plot(coords1(1), coords1(2), '*b')
                hold on; plot(coords2(1), coords2(2), '*b')
        viscircles(centers(1:2,:), radii(1:2),'EdgeColor','b');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Take each peak pair and find the distance between them as well as
        % the angle they make with the horizontal axis to determine the
        % depth. The midpoint will be the localized co-ordinate in 2D.
        
        % Measure angles w.r.t. the horizontal axis
        % dot(v1,v2) = mod(v1)mod(v2)cos(angle)
        % X, Y are vectors pointing from one peak to another
        if coords1(1) > coords2(1) & coords1(2) > coords2(2)
            X = coords1(1) - coords2(1); Y = coords1(2) - coords2(2);
            
        elseif coords1(1) > coords2(1) & coords1(2) < coords2(2)
            X = coords1(1) - coords2(1); Y = coords2(2) - coords1(2);
            
        elseif coords1(1) < coords2(1) & coords1(2) > coords2(2)
            X = coords2(1) - coords1(1); Y = coords1(2) - coords2(2);
            
        else
            X = coords2(1) - coords1(1); Y = coords2(2) - coords1(2);
        end
        
        % Normalized vector between the lobes
        pairVect = [X, Y];
        pairLen = sqrt( X^2 + Y^2 );
        pairVect = pairVect ./ pairLen;

        % Unit vector in horizontal direction
        unitVect = refVect;
        
        % Use dot product to get the angle
        term = dot(pairVect, unitVect);
        ang = acosd(term);
        midpt = [coords1(1) + coords2(1), coords1(2) + coords2(2)] ./ 2;
        
        % Store the angle between the pair corresponding to each defocus
        % value. Also store the midpoint between them, which corresponds to
        % the localised 2D position of the particle.
        angles = [angles; ang];
        midPts = [midPts; midpt];
    
    end
end
    