% --- General bosonic reference wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioBECR(BECRObj,Diff)
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

Update.Occ = BECRObj.Occ; Update.Nb = BECRObj.Nb;
Ratio = prod(BECRObj.SPO(Diff.pos).'.^Diff.val);
for d = 1:Diff.num
    Update.Occ(Diff.pos(d)) = Update.Occ(Diff.pos(d)) + Diff.val(d);
    Update.Nb = Update.Nb + Diff.val(d);
    for m = 1:abs(Diff.val(d))
        Ratio = Ratio * sqrt(m+min([Update.Occ(Diff.pos(d)),BECRObj.Occ(Diff.pos(d))]))^(-sign(Diff.val(d)));
    end
end
for d = 1:abs(sum(Diff.val))
    Ratio = Ratio * sqrt(min([Update.Nb,BECRObj.Nb]))^(-sign(sum(Diff.val(d))));
end
