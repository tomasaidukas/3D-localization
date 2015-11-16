function result = Gaussian(C, X, Y)
%------------------------------------------------------------%
% Function for a 2D Gaussian
% C1 Normalization term.
% C4 standard deviation.
% C2, C3 centre of the lobe.
%------------------------------------------------------------%

    result = C(1) .* exp(-...
            ( (X(:) - C(2)).^2 ./ (2 * C(4).^2) + ...
              (Y(:) - C(3)).^2 ./ (2 * C(5).^2) ));
end
