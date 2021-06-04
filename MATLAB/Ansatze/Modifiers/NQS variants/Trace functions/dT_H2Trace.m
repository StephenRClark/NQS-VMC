function [dF] = dT_H2Trace(Theta,Phi,HDim)
% Performs logarithmic derivative of trace w.r.t. hidden-linear term Theta
% over number-like hidden unit h for an NQS with square bias terms. Nmax =
% HDim - 1, Theta is linear in hidden units and Phi is quadratic in hidden
% units.
% N.B: assumes Theta is a (Nh x 1) vector and the same for Phi.
n = 0:(HDim - 1); M = exp((Theta .* n) + (Phi.*(n.^2)));
% First index is hidden index, second is hidden unit value from 0 to Nmax.
dF = sum(((ones(numel(Theta),1)*n) .* M),2)./sum(M,2);
% Output dF is same size as Theta / Phi.
end