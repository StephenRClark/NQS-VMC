% --- Fermionic Slater determinant reference state logarithmic derivative function ---

function dLogp = LogDerivSDet(SDetObj,HilbertObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the fermionic reference ansatz, for a fermion
% configuration specifed by the structure Cfg.
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
% Format for dLogp is a Np x 1 vector of parameter derivatives.
% ---------------------------------
% N.B. First parameter in HVar will automatically be excluded from the
% optimisation and act as the energy scale of the reference, as a
% preventative measure against runaway HVar scaling. It is recommended
% to set the first term in HVar as the hopping term.
% ---------------------------------

% Make local copies to reduce notation in code below.
N = HilbertObj.N; % Number of sites.
FermLoc = SDetObj.FermLoc; DetMat = SDetObj.DetMat;
Orbitals = SDetObj.Orbitals; N_up = numel(Cfg.up); N_dn = numel(Cfg.dn);
Nf = N_up + N_dn;

OrbIndsE = [(N_up+1):N, ((N_dn+1):N)+N];
OrbIndsF = [1:N_up, (1:N_dn)+N];

dLogp = zeros(SDetObj.Np,1); % Initialise full vector of derivatives.

for h = 1:SDetObj.Np
    if SDetObj.OptInds(p) ~= 0
        WRH = SDetObj.WArr(:,:,h); % Extract necessary connectivity matrix.
        WRH = Orbitals(:,OrbIndsE) * WRH(OrbIndsE,OrbIndsF) * (Orbitals(:,OrbIndsF)');
        for f = 1:Nf
            sites = 1:N;
            % Only DetMat matrix elements that represent hops from currently
            % occupied sites will contribute.
            if f <= N_up
                sites(Cfg.up) = [];
            else
                sites(Cfg.dn) = []; sites = sites + N;
            end
            dLogp(h) = dLogp(h) - sum(WRH((FermLoc==f),sites).*DetMat(sites,f).');
        end
    end
end

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
