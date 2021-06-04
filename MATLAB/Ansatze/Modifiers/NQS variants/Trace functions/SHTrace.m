function [F] = SHTrace(Theta,B,HDim)
% Performs trace over number-like hidden unit h for an NQS with square bias
% terms. Nmax = HDim - 1, B is hidden unit square bias.
% N.B: assumes Theta is a (Nh x 1) vector and the same for B.
n = (-(HDim-1):2:(HDim - 1)) / 2; % 1 x HDim vector.
% In order to curb possible numerical instabilities, may have to divide
% through by exp((B(Nmax)^2)/2).
F = sum( exp((ones(numel(Theta),1)*n) .* ((Theta * ones(1,HDim)) + (B*n))),2)/HDim;% ...
     %./ exp(1.*(B*ones(1,HDim)) * (HDim-1)^2) ,2);
% Output F is same size as Theta / B.
end