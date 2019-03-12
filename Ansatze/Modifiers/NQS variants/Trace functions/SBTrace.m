function [F] = SBTrace(Theta,HDim,B)
% Performs trace over number-like hidden unit h for an NQS with square bias
% terms. Nmax = HDim - 1, B is hidden unit square bias.
% N.B: assumes Theta is a (Nh x 1) vector and the same for B.
n = 0:(HDim - 1); % 1 x HDim vector.
% In order to curb possible numerical instabilities, will divide through by
% exp((B(Nmax)^2)/2).
F = sum( exp((ones(numel(Theta),1)*n) .* ((Theta * ones(1,HDim)) + (B*n))) ...
    ./ exp( 0.5 * (B*ones(1,HDim)) * (HDim-1)^2) ,2);
% Output F is same size as Theta / B.
end