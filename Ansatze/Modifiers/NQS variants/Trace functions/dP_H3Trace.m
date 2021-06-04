function [dF] = dP_H3Trace(Theta,Phi,Omega,HDim)
% Performs logarithmic derivative of trace w.r.t. hidden-squared term Phi
% over number-like hidden unit h for an NQS with square bias terms. Nmax =
% HDim - 1, Theta is linear in hidden units, Phi is quadratic in hidden
% units and Omega is cubic in hidden units.
% N.B: assumes Theta is a (Nh x 1) vector and the same for Phi and Omega.
n = 0:(HDim - 1); M = exp((Theta .* n) + (Phi.*(n.^2)) + (Omega.*(n.^3)));
% First index is hidden index, second is hidden unit value from 0 to Nmax.
dF = sum(((ones(numel(Theta),1)*(n.^2)) .* M),2)./sum(M,2);
% Output dF is same size as Theta / Phi / Omega.
end