function [dF] = dB_NHTrace(Theta,B,HDim)
% Performs logarithmic derivative of trace w.r.t B over number-like hidden
% unit h for an NQS with square bias terms. Nmax = HDim - 1, B is hidden
% unit square bias.
% N.B: assumes Theta is a (Nh x 1) vector and the same for B.
n = 0:(HDim - 1);
% In order to curb possible numerical instabilities, might divide through by
% exp((B(Nmax)^2)/2). This will crop up in the derivative as well.
M = exp((ones(numel(Theta),1)*n) .* ((Theta * ones(1,HDim)) + (B*n))); % ...
    % ./ exp( 0.5 * (B*ones(1,HDim)) * (HDim-1)^2);
% First index is hidden index, second is hidden unit value from 0 to Nmax.
dF = sum((ones(numel(Theta),1)*(n.^2) .* M),2)./sum(M,2);
    % - (0.5*(HDim-1)^2)*ones(numel(Theta),HDim)) 
% Output dF is same size as Theta / B.

end