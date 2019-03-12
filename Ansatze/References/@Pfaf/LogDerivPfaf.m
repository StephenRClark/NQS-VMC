% --- General fermionic Pfaffian logarithmic derivative function ---

function dLogp = LogDerivPfaf(PfafObj,HilbertObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the fermionic Pfaffian ansatz, for a fermionic
% configuration specifed by the structure Cfg.
% ---------------------------------
% Format for fermionic Pfaffian Reference:
% - Pfaf.Nf = (1 x 2) vector - number of up/down fermions respectively.
% - Pfaf.PairMat = (2N x 2N) matrix - contains all pairing terms.
% - Pfaf.PfI = (Nf x Nf) matrix - inverse of reduced PfFull matrix.
% - Pfaf.PfG = (2N x Nf) matrix - matrix used for ratio calculations.
% - Pfaf.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
% - Pfaf.Np = number of variational parameters associated with Pfaf Reference.
% Pfaf properties used in variational version:
% - Pfaf.PfV = (2N x 2N) array - logs which variational parameters make up the elements of Pfaf.PairMat.
% - Pfaf.PfVR = (Nf x Nf) array - reduced matrix constructed from PfV.
% - Pfaf.PfVar = (Np x 1) vector - variational parameters in PfFull.
% ---------------------------------
% Format for dLogp vector is a (Np x 1) vector of parameter derivatives.
% ---------------------------------

% Make local copies to reduce notation in code below.
PfI = PfafObj.PfI; PfVR = PfafObj.PfVR; Np = PfafObj.Np;
dLogp = zeros(Np,1);

% LogDerivPfaf depends highly on PfV being correct - worth explicitly reconstructing.
Cfg_vec = reshape(HilbertObj.FullCfgRef(Cfg),numel(PfafObj.FermLoc),1);
PfafObj.PfVR = PfafObj.PfV(Cfg_vec==1,Cfg_vec==1);

for p = 1:Np
    % Construct dPfR/dp matrix by locating matrix elements featuring parameter p:
    dPfR = (PfVR == p) - (PfVR == -p);
    % Pfaffian logarithmic derivative:
    if isempty(dPfR) == 0
        dLogp(p) = (1/2) * trace(PfI * dPfR);
    end
end

% Nf x Nf elements appear in PfRed - dLogp should have derivatives for all
% the non-zero elements of PfRed

dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;