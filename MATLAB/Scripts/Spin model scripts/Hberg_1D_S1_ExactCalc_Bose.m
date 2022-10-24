% ------------------------------------
% Exact code for a bosonic lattice for fixed particle number
% ------------------------------------

% System setup :
L = 14; % Number of lattice sites.
N = 14; % Number of atoms.
N_m = 2; % Maximum occupation on a lattice site (set N_m = 1 for hard-core limit).
not_fermion = 0; % Computing bosons only so no fermionic phases.
bc = 0; % (1 = O.B.C., 0 = P.B.C.)

% --- Construct bosonic operators :

% NOTE: This step can be very time-consuming - consider computing once and saving
% the output if multiple calculations for a given system size are needed.

% Compute the operators and basis :
tic
[basis,W] = num_basis(L,N,N_m); % Full number state basis in the chosen total particle number sector.
[hop,num,list] = bose_operators(L,basis,W,not_fermion); % Build the commonly used operators.
toc
% ---

basis = basis - 1;

% Pick system geometry (comment/uncomment):
ind = nn_oned(L,list,bc); % For chain.
% ind = nn_twod(L,list,bc); % For square lattice.
list(ind,:); % Check the pairs nearest-neighbour hopping terms retained.

JH = 1;

% BHM parameters :
J = JH; % Hopping amplitude.
U = JH; % On-site interaction strength.

% N.N Hopping :
HJ = sparse(1,1,0,W,W,W); % W = basis size for this sector.
for m=1:length(ind)
    HJ = HJ + J*hop{ind(m)}; % Add up all the relevant hopping terms.
end
% On-site Interaction :
HU = sparse(1,1,0,W,W,W);
for m=1:L
    n = 1 + mod(m,L);
    HU = HU + U*(num{m}-speye(W))*(num{n}-speye(W)); % Use speye(W) and not just 1 since the latter is not the identity but a matrix of ones!.
end
% Full Hamiltonian :
H = HJ + HU;

tic
opts.disp = 0;
[gs,gs_en] = eigs(H,1,'sr',opts); % Compute ground state only.
toc