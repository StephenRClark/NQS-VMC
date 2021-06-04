function [Params] = paramListCPSTI(NQSObj)
Nv = NQSObj.Nv;
Alpha = NQSObj.Alpha;
Params = zeros(NQSObj.Np,1);
Params(1:2) = NQSObj.ati;
Params((1:(Alpha * 2)) + 2) = NQSObj.bti;
Params(2 + Alpha * 2 + (1:(4 * Nv * Alpha))) = NQSObj.Wti;
end
