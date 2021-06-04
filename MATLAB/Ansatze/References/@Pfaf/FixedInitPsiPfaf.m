% --- General fermionic Pfaffian initialisation function ---

function [PfafObj] = FixedInitPsiPfaf(PfafObj,Params)
% This function constructs the fermionic Pfaffian Reference object. The
% input Pfaf object is assumed to have N, Nf and Np already assigned. The
% Params structure contains information controlling the form of the
% reference Hamiltonian generated, and subsequently the form of the
% fermionic pair amplitudes.
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

% Pfaf can be initialised in two ways - from the eigenstates of a quadratic
% Hamiltonian, and with random parameters with symmetries imposed by Graph.
% First option will initialise with FixedInitPsiPfaf, though if
% variational, will then take the individual terms and optimise without
% regard for initial Hamiltonian parameters.

% Make local copies to reduce notation in code below.
N = PfafObj.N; % Number of sites.
Nf = PfafObj.Nf; % 2 x 1 vector of fermion numbers.

% Construct quadratic Hamiltonian from initial HVar provided in
% Params as well as their connectivity arrays CArr.
TMat = zeros(2*N);
for h = 1:numel(Params.HVar)
    TMat = TMat + Params.HVar(h)*Params.CArr(:,:,h);
end

% Calculate eigenstates of newly constructed matrix for use
% with Slater determinant, and rearrange according to spin
% expectation then ascending energy.
[Orbitals,~] = eig(TMat);
OrbMatUp = Orbitals(:,(sum(abs(Orbitals(1:N,:)),1)>=0.5)); % Collect all 'up-spin' orbitals.
OrbMatDn = Orbitals(:,(sum(abs(Orbitals((N+1):end,:)),1)>=0.5)); % Collect all 'down-spin' orbitals.
% Recollect all Orbitals together as one matrix.
Orbitals = [OrbMatUp OrbMatDn];

% Construct full pairing term matrix from Orbitals.
PairMat = zeros(2*N);
for n = 1:(sum(Nf)/2)
    PairMat = PairMat + (Orbitals(:,n)*(Orbitals(:,n+N)')) - (Orbitals(:,n+N)*(Orbitals(:,n)'));
end
PairMat(abs(PairMat)<1e-12) = 0;
PfafObj.PairMat = PairMat;
% Find non-zero values and log their locations.
Np = 0; PfV = zeros(2*N); PfVar = zeros(N*(2*N-1),1); % Maximum number of potential parameters.
while sum(PairMat(:)~=0) ~= 0 % Check until all non-zero terms have been catalogued and removed.
    P = PairMat(find(abs(PairMat)>1e-15,1));
    Np = Np + 1; PfVar(Np) = Np;
    % Locate terms that are equal to P and -P, then mark by parameter
    % number in PfV.
    PfV(abs(PairMat - P)<1e-15) = Np;
    PfV(abs(PairMat + P)<1e-15) = -Np;
    % Set these entries in PairMat to zero.
    PairMat(abs(PairMat - P)<1e-15) = 0;
    PairMat(abs(PairMat + P)<1e-15) = 0;
end

PfVar = PfVar(1:Np); PfafObj.PfVar = PfVar; PfafObj.PfV = PfV; PfafObj.Np = Np;

% Assigning placeholders for configuration dependent information.
PfafObj.FermLoc = zeros(2*N,1); PfafObj.PfVR = zeros(sum(Nf));
PfafObj.PfI = zeros(sum(Nf)); PfafObj.PfG = zeros(2*N,sum(Nf));

PfafObj.OptInds = ones(PfafObj.Np,1);