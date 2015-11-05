function result = DoubleGaussian(C, X, Y)
%------------------------------------------------------------%
% Function for a 2D double Gaussian.
% C1, C5 Normalization terms.
% C2, C3, C6, C7 peak position of the Gaussian peaks.
% C4, C8 standard deviations of each lobe.
%------------------------------------------------------------%

result = C(1) .* exp(-...
        ( (X(:) - C(2)).^2 ./ (2 * C(4).^2) + ...
          (Y(:) - C(3)).^2 ./ (2 * C(4).^2) )) + ...
        C(5) .* exp(-...
        ( (X(:) - C(6)).^2 ./ (2 * C(8).^2) + ...
          (Y(:) - C(7)).^2 ./ (2 * C(8).^2) ));
end