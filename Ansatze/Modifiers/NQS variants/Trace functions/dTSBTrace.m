function [dF] = dTSBTrace(Theta,HDim,B)
% Performs logarithmic derivative of trace w.r.t Theta over number-like
% hidden unit h for an NQS with square bias terms. Nmax = HDim - 1, B is
% hidden unit square bias.
% N.B: assumes Theta is a (Nh x 1) vector and the same for B.
n = 0:(HDim - 1);
% In order to curb possible numerical instabilities, might divide through by
% exp((B(Nmax)^2)/2).
M = exp((ones(numel(Theta),1)*n) .* ((Theta * ones(1,HDim)) + (B*n))); % ...
    % ./ exp( 0.5 * (B*ones(1,HDim)) * (HDim-1)^2);
% First index is hidden index, second is hidden unit value from 0 to Nmax.
dF = sum(((ones(numel(Theta),1)*n) .* M),2)./sum(M,2);
% Output dF is same size as Theta / B.

end