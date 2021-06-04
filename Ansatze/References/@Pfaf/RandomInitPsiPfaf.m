% --- General fermionic Pfaffian initialisation function ---

function [PfafObj] = RandomInitPsiPfaf(PfafObj,Params)
% This function constructs the fermionic Pfaffian Reference object. The
% input Pfaf object is assumed to have N and Nf already assigned. The
% Params structure contains information controlling the form of the
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
% Second option will require Graph with desired symmetry, and will populate
% a matrix with small random parameters.

% Make local copies to reduce notation in code below.
N = PfafObj.N; % Number of sites.
Nf = PfafObj.Nf; % 2 x 1 vector of fermion numbers.
GraphObj = PfafObj.Graph; BondMap = GraphObj.BondMap;

% Construct antisymmetric matrix of parameter numbers PfV according to the
% symmetry of the Graph and the Pfaffian antisymmetry requirement.

PfV = zeros(2*N); Np = 0; % Count number of parameters needed to populate PfV.
if sum(Nf) == N && Nf(1) == Nf(2) % Half-filling case, will only populate terms with dissimilar spin.
    % Overall matrix structure will be PfV = [0 pfv; -pfv 0].
    for n = 1:N
        % Antisymmetry eliminates diagonal values of PfV.
        for m = (1:N) + N % Symmetry requirement means initially only need to populate on one side of the diagonal.
            if PfV(n,m) == 0 % If not populated with a variable index, needs assignment.
                Np = Np + 1;
                PfV(n,m) = Np; PfV(m,n) = -Np;
                for b = 1:numel(BondMap) % Assign variable indices for each translate.
                    Inds = [(BondMap{b}(1+mod(n-1,N))+N*(ceil(n/N)-1)),...
                        (BondMap{b}(1+mod(m-1,N))+N*(ceil(m/N)-1))];
                    if Inds(1) ~= 0 && Inds(2) ~= 0
                        if PfV(Inds(1),Inds(2)) == 0 && PfV(Inds(2),Inds(1)) == 0
                            PfV(Inds(1),Inds(2)) = Np; PfV(Inds(2),Inds(1)) = -Np;
                        end
                    end
                end
            end
        end
    end
else % Populate full 2N x 2N PfV matrix.
    for n = 1:(2*N-1)
        % Antisymmetry eliminates diagonal values of PfV.
        for m = (n+1):2*N % Symmetry requirement means initially only need to populate on one side of the diagonal.
            if PfV(n,m) == 0 % If not populated with a variable index, needs assignment.
                Np = Np + 1;
                PfV(n,m) = Np; PfV(m,n) = -Np;
                for b = 1:numel(BondMap) % Assign variable indices for each translate.
                    Inds = [(BondMap{b}(1+mod(n-1,N))+N*(ceil(n/N)-1)),...
                        (BondMap{b}(1+mod(m-1,N))+N*(ceil(m/N)-1))];
                    if Inds(1) ~= 0 && Inds(2) ~= 0
                        if PfV(Inds(1),Inds(2)) == 0 && PfV(Inds(2),Inds(1)) == 0
                            PfV(Inds(1),Inds(2)) = Np; PfV(Inds(2),Inds(1)) = -Np;
                        end
                    end
                end
            end
        end
    end
end

PfVar = zeros(Np,1); % Generate the random parameters that will populate PairMat.
for p = 1:Np
    PfVar(p) = Params.Pf * (1 - Params.nmag + 2*Params.nmag*rand) * exp(2i*pi*Params.nphs*rand);
end

PfafObj.PfVar = PfVar; PfafObj.PfV = PfV; PfafObj.Np = Np;
% Construct full pairing term matrix using PfVar and PfV.
PairMat = zeros(2*N);
PairMat((PfV>0)) = PfVar(PfV(PfV>0)); PairMat((PfV<0)) = - PfVar(-PfV(PfV<0));
PfafObj.PairMat = PairMat;

% Assigning placeholders for configuration dependent information.
PfafObj.FermLoc = zeros(2*N,1); PfafObj.PfVR = zeros(sum(Nf));
PfafObj.PfI = zeros(sum(Nf)); PfafObj.PfG = zeros(2*N,sum(Nf));

PfafObj.OptInds = ones(PfafObj.Np,1);
