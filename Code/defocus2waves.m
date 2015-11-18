function [deviation] = defocus2waves(R, f, XYrange, W20);
%------------------------------------------------------------%
% Longitudinal length shift in terms of waves
%------------------------------------------------------------%

diameter = (R * 2);
fnumber = f / diameter;

deviation = 8 .* W20 .* fnumber.^2;

end