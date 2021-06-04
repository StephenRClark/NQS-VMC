function [F] = H3Trace(Theta,Phi,Omega,HDim)
% Performs trace over number-like hidden unit h for an NQS with
% hidden-squared terms. Nmax = HDim - 1, Theta is linear in hidden units,
% Phi is quadratic in hidden units and Omega is cubic in hidden units.
% N.B: assumes Theta is a (Nh x 1) vector and the same for Phi and Omega.
n = 0:(HDim - 1); % 1 x HDim vector.
F = sum( exp((Theta .* n) + (Phi.*(n.^2)) + (Omega.*(n.^3))),2)/HDim;
% Output F is same size as Theta / Phi / Omega.
end