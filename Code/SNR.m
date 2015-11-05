function [SNR] = SNR(img, var)
%------------------------------------------------------------%
% Find the SNR of the image. var is the noise variance.
%------------------------------------------------------------%

SNR = (sum(img(:)) ./ numel(img)) ./ sqrt(var);

end