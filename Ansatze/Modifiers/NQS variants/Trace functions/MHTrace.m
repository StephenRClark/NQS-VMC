function [F] = MHTrace(ThetaH,ThetaD,HDim)
% Performs trace over number-like hidden unit h for an NQS expressed with
% doublon and holon operators.
% N.B: assumes ThetaH / ThetaD is a (Nh x 1) vector.
n = 0:(HDim - 1); D = n - 1; D(1) = 0; H = [1 zeros(1,HDim-1)];
% Expressing doublons and holons in this form ensures expression is correct
% for all values of n, but does result in funky values for n > 2.
F = sum(exp(ThetaH.*H + ThetaD.*D),2)/HDim; 
% Output F will be a (Nh x 1) vector.
end