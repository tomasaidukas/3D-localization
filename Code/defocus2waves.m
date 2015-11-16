function [deviation] = defocus2waves(camera);
%------------------------------------------------------------%
% Longitudinal length shift in terms of waves
%------------------------------------------------------------%
[W20, maxDefocus, NoPts, XYrange, R, f] = camera{:};

diameter = (R * 2);
fnumber = f / diameter;

deviation = 8 .* W20 .* fnumber.^2;

end