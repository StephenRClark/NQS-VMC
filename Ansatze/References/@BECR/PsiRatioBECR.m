% --- General bosonic reference wave function amplitude ratio function ---

function [Ratio,OccP] = PsiRatioBECR(BECRObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed bosonic configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for Bose Einstein Condensate Reference:
% - BECR.Nb = number of bosons - assumed fixed for the most part.
% - BECR.SPO = (Nv x 1) vector - single particle boson orbital amplitudes.
% - BECR.Occ = (Nv x 1) vector - onsite boson occupation numbers.
% ---------------------------------
% Format for Update is a vector containing the new boson occupation numbers.
% ---------------------------------

OccP = BECRObj.Occ;
Ratio = prod(BECRObj.SPO(Diff.pos).'.^Diff.val);
for d = 1:Diff.num
    Ratio = Ratio * sqrt((OccP(Diff.pos(d))+(Diff.val(d)>0))^(-Diff.val(d)));
    OccP(Diff.pos(d)) = OccP(Diff.pos(d)) + Diff.val(d);
end
