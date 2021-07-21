% --- General fermionic determinant wave function preparation function ---

function [Ansatz] = PrepPsiDet_BfCr(Ansatz,Cfg)
% This function initialises the fermionic determinant ansatz structure 
% intermediate information (matrices) given an initial configuration.
% BfCr - backflow correlations, will alter determinant according to
% doublon-holon positions.
% ---------------------------------
% Format for fermionic determinant wavefunction addition:
% - Ansatz.Nf = number of fermions - set to Nv for spin models, specified at initialisation.
% - Ansatz.UFull = (2Nv x 2Nv) matrix - contains all available single particle orbitals.
% - Ansatz.UFe = (2Nv x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
% - Ansatz.WFe = (2Nv x Nf) matrix - elements are used for determinants in PsiRatio.
% - Ansatz.FermLoc = (2Nv x 1) vector - details locations of fermions by index for sign tracking purposes.
% Backflow correlation alterations:
% - Ansatz.UInv (Nf x 2Nv) replaces WFe as it is more convenient to work with.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = Ansatz.Nv; % Number of "visible" spins / sites.
Nf = Ansatz.Nf; % Number of fermions.
UMatUp = Ansatz.UFull(:,(sum(abs(Ansatz.UFull(1:Nv,:)),1)~=0)); % Collect all up-spin orbitals
UMatDn = Ansatz.UFull(:,(sum(abs(Ansatz.UFull((Nv+1):end,:)),1)~=0)); % Collect all down-spin orbitals
FermLoc = zeros(2*Nv,1);
UFe = zeros(2*Nv,Nf);
for u = 1:numel(Cfg.up)
    FermLoc(Cfg.up(u)) = u; % Fermions numbering 1 to N_up
    % Location of up spin u is the location where OpOrd == u
    UFe(:,u) = UMatUp(:,u); % Fill up single particle orbitals with up-spin fermions
end
for d = 1:numel(Cfg.dn)
    FermLoc(Cfg.dn(d) + Nv) = d + numel(Cfg.up); % Fermions numbering N_up+1 to N_f
    % Location of down spin d is the location where OpOrd == d, minus Nv
    UFe(:,d + numel(Cfg.up)) = UMatDn(:,d); % Fill up single particle orbitals with down-spin fermions
end
Ansatz.UFe = UFe;
UInv = UFe^(-1);
Ansatz.UInv = UInv;

Ansatz.FermLoc = FermLoc;

% Sanity check for any unreasonably large or small elements
cap = 1e30; min = 1e-30;
Ansatz.UFe(abs(Ansatz.UFe)>cap)=0;
Ansatz.UFe(abs(Ansatz.UFe)<min)=0;
Ansatz.UInv(abs(Ansatz.UInv)>cap) = 0;
Ansatz.UInv(abs(Ansatz.UInv)<min) = 0;