% --- General Slater determinant fermionic reference initialisation function ---

function [SDetObj] = FixedInitPsiSDet(SDetObj,Params)
% This function constructs the Slater determinant Reference object. The
% input SDet object is assumed to have N, Nf and Np already assigned. The
% Params structure contains information controlling the form of the
% reference Hamiltonian generated, and subsequently the form of the single
% particle orbitals.
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
N = SDetObj.N; % Number of "visible" spins.

% Construct quadratic Hamiltonian from initial HVar provided in
% Params as well as their connectivity arrays CArr.
TMat = zeros(2*N);
for h = 1:numel(Params.HVar)
    TMat = TMat + Params.HVar(h)*Params.CArr(:,:,h);
end

% Assign hidden properties for variational use.
SDetObj.HVar = Params.HVar; SDetObj.CArr = Params.CArr;
SDetObj.WArr = SDetObj.CArr; % Placeholder for WArr, which is a transform of CArr.
SDetObj.FermLoc = zeros(2*N,1); % Placeholder to be populated later.

% Calculate eigenstates of newly constructed matrix for use
% with Slater determinant, and rearrange according to spin
% expectation then ascending energy. Additionally, require some
% energy prefactors for later use in LogDeriv.
[Orbitals,En] = eig(TMat); En = diag(En);
OrbMatUp = Orbitals(:,(sum(abs(Orbitals(1:N,:)),1)>=0.5)); % Collect all 'up-spin' orbitals.
OrbMatDn = Orbitals(:,(sum(abs(Orbitals((N+1):end,:)),1)>=0.5)); % Collect all 'down-spin' orbitals.
% Rearrange energies in same manner as Orbitals.
En = [En((sum(abs(Orbitals(1:N,:).^2),1)>1/2)); En((sum(abs(Orbitals((N+1):end,:)).^2,1)>1/2))];
% Recollect all Orbitals together as one matrix.
Orbitals = [OrbMatUp OrbMatDn]; SDetObj.Orbitals = Orbitals;
EnFac = zeros(2*N);
for n = 1:2*N
    for m = 1:2*N
        if abs(En(n) - En(m))>1e-15
            % Generate energy prefactors for LogDeriv function.
            EnFac(n,m) = 1/(En(n) - En(m));
        end
    end
end
SDetObj.EnFac = EnFac; 
SDetObj.OptInds = [ones(SDetObj.Np,1)]; 
% Set the first entry in HVar as the energy scale of the reference.