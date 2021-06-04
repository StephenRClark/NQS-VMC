function [F] = H2Trace(Theta,Phi,HDim)
% Performs trace over number-like hidden unit h for an NQS with
% hidden-squared terms. Nmax = HDim - 1, B is hidden unit square bias,
% Theta is linear in hidden units and Phi is quadratic in hidden units.
% N.B: assumes Theta is a (Nh x 1) vector and the same for Phi.
n = 0:(HDim - 1); % 1 x HDim vector.
F = sum( exp((Theta .* n) + (Phi.*(n.^2))),2)/HDim;
% Output F is same size as Theta / Phi.
end