function [dTH, dTM] = dT_MHTrace(ThetaH,ThetaM,HDim)
% Performs logarithmic derivative of trace over hidden unit h for an NQS
% expressed with multiplon and holon operators.
% N.B: assumes ThetaH / ThetaM is a (Nh x 1) vector.
n = 0:(HDim - 1); M = n - 1; M(1) = 0; H = [1 zeros(1,HDim-1)]; Nh = numel(ThetaH);
% Expressing multiplons and holons in this form ensures expression is
% correct for all values of n. Alternative definitions may mean the
% operators cannot be interpreted as bona fide holons and multiplons.
F = exp(ThetaH.*H + ThetaM.*M);
dTH = sum(F .* (ones(Nh,1) .* H),2)./sum(F,2); 
dTM = sum(F .* (ones(Nh,1) .* M),2)./sum(F,2);
% Outputs dT will be (Nh x 1) vectors.
end