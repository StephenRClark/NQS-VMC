% --- General fermionic determinant wave function preparation function ---

function [SDetObj] = PrepPsiSDet(SDetObj,Cfg)
% This function initialises the fermionic determinant ansatz structure
% intermediate information (matrices) given an initial configuration.
% ---------------------------------
% Format for SDet Reference:
% - SDet.Nf = (1 x 2) vector - number of up/down fermions respectively.
% - SDet.Orbitals = (2N x 2N) matrix - contains all available single particle orbitals.
% - SDet.OrbMat = (2N x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
% - SDet.DetMat = (2N x Nf) matrix - elements are used for determinants in PsiRatio.
% - SDet.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
% - SDet.Np = number of variational parameters associated with SDet Reference.
% SDet properties used in variational version:
% - SDet.CArr = (2N x 2N x Np) array - connectivity array for the reference Hamiltonian.
% - SDet.WArr = (2N x 2N x Np) array - transformed connectivity array for the reference Hamiltonian.
% - SDet.EnFac = (2N x 2N) matrix - elements are used in LogDeriv function.
% - SDet.HVar = (Np x 1) vector - variational parameters in the reference Hamiltonian used.
% ---------------------------------

% Make local copies to reduce notation in code below.
N = SDetObj.N; % Number of sites.
SDetObj.Nf = [numel(Cfg.up) numel(Cfg.dn)]; % Assign Nf - may change later.
Nf = sum(SDetObj.Nf); % Number of fermions.

FermLoc = zeros(2*N,1);
OrbMat = zeros(2*N,Nf);
for u = 1:numel(Cfg.up)
    FermLoc(Cfg.up(u)) = u; % Fermions numbering 1 to N_up.
    OrbMat(:,u) = SDetObj.Orbitals(:,u);
    % Fill up single particle orbitals with up-spin fermions.
end
for d = 1:numel(Cfg.dn)
    FermLoc(Cfg.dn(d) + N) = d + numel(Cfg.up); % Fermions numbering N_up+1 to N_f.
    OrbMat(:,d + numel(Cfg.up)) = SDetObj.Orbitals(:,N+d);
    % Fill up single particle orbitals with down-spin fermions
end
SDetObj.OrbMat = OrbMat;
SDetObj.FermLoc = FermLoc;
OrbMatRed = zeros(Nf);
for l = 1:length(FermLoc)
    if FermLoc(l) ~= 0
        OrbMatRed(FermLoc(l),:) = SDetObj.OrbMat(l,:);
        % Construct reduced U with order determined by FermLoc.
    end
end
SDetObj.DetMat = OrbMat * (OrbMatRed^(-1));

% If variational, prepare SDet properties used in LogDeriv.
if SDetObj.VFlag == 1
    for h = 1:numel(SDetObj.HVar)
        % WArr(:,:,h) represents the connectivity matrix of parameter h
        % transformed into the basis that diagonalises the Hamiltonian.
        SDetObj.WArr(:,:,h) = (SDetObj.Orbitals*SDetObj.CArr(:,:,h)*SDetObj.Orbitals').*SDetObj.EnFac;
    end
end

% Sanity check for any unreasonably large or small elements
cap = 1e10; min = 1e-10;
SDetObj.OrbMat(abs(SDetObj.OrbMat)>cap)=0;
SDetObj.OrbMat(abs(SDetObj.OrbMat)<min)=0;
SDetObj.DetMat(abs(SDetObj.DetMat)>cap) = 0;
SDetObj.DetMat(abs(SDetObj.DetMat)<min) = 0;