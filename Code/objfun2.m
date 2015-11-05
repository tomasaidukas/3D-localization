function resids = objfun(params, X, Y, F)
%------------------------------------------------------------%
% Function used for least square fitting.
%------------------------------------------------------------%
ans = Gaussian(params, X, Y); 
resids = ans- F(:);

end