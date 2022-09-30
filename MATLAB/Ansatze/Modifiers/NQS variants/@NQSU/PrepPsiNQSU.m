% --- General NQS wave function preparation function ---

function [NQSObjNew] = PrepPsiNQSU(NQSObjNew,TestCfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" units.
% - NQSU.Nh = number of "hidden" units.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------
% Format for Update is a vector of new effective angles ThetaP and
% new unary vector UVecP.
% ---------------------------------

[Cfg_vec] = NQSObjNew.FullCfg(TestCfg); % Convert Cfg structure into a spin configuration vector.
Nmax = NQSObjNew.VDim-1; UVec = zeros(Nmax,NQSObjNew.Nv);
for v = 1:Nmax
    UVec(v,:) = (Cfg_vec.' == NQSObjNew.VList(v+1));
end
UVec = reshape(UVec,Nmax*NQSObjNew.Nv,1);
NQSObjNew.UVec = UVec;
NQSObjNew.Theta = NQSObjNew.b + NQSObjNew.W*UVec; % Compute the effective angle for this configuration state.