% --- General fermionic Pfaffian wave function preparation function ---

function [PfafObj] = PrepPsiPfaf(PfafObj,Cfg)
% This function initialises the fermionic Pfaffian ansatz structure 
% intermediate information (matrices) given an initial configuration.
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

% Make local copies to reduce notation in code below.
N = numel(PfafObj.FermLoc)/2; % Number of sites.
Nf = sum(PfafObj.Nf); % Number of fermions.

FermLoc = zeros(2*N,1);
for u = 1:numel(Cfg.up)
    FermLoc(Cfg.up(u)) = u; % Fermions numbering 1 to N_up
    % Location of up spin u is the location where FermLoc == u
end
for d = 1:numel(Cfg.dn)
    FermLoc(Cfg.dn(d) + N) = d + numel(Cfg.up); % Fermions numbering N_up+1 to N_f
    % Location of down spin d is the location where FermLoc == d, minus Nv
end

PfafObj.FermLoc = FermLoc;

PfInds = zeros(Nf,1);
for i = 1:Nf
    PfInds(i) = find(PfafObj.FermLoc == i);
end
PfRed = PfafObj.PairMat(PfInds,PfInds); 
PfafObj.PfVR = PfafObj.PfV(PfInds,PfInds,:);
PfafObj.PfG = PfafObj.PairMat(:,PfInds) * (PfRed^(-1));
PfafObj.PfI = PfRed^(-1);

% Sanity check for any unreasonably large or small elements
cap = 1e30; min = 1e-30;
PfafObj.PfG(abs(PfafObj.PfG)>cap)=0;
PfafObj.PfG(abs(PfafObj.PfG)<min)=0;
PfafObj.PfI(abs(PfafObj.PfI)>cap) = 0;
PfafObj.PfI(abs(PfafObj.PfI)<min) = 0;