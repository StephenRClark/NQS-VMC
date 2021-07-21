function [dF] = dL_SqTrace(Theta,ThetaSq,HVals)
% Performs trace over number-like hidden unit h for an NQS with square bias
% terms. Nmax = HDim - 1, B is hidden unit square bias.
% N.B: assumes Theta / ThetaSq is a (Nh x 1) vector, HVals is a 1 x HDim vector.
HVSq = HVals.^2; FLin = Theta * HVals; FSq = ThetaSq * HVSq;
dF = sum((HVals.*exp(FLin + FSq)),2)./sum(exp(FLin + FSq),2);
% Output dF is same size as Theta / ThetaSq.
end