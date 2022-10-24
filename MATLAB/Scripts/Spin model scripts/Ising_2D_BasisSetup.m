L = 16; N = L; S = 1/2; SiteDim = 2*S + 1;

% This section will generate all the configurations possible for the
% specified lattice.
IndAncilla = ones(1,N); Nconf = SiteDim^N; S0Inds = zeros(Nconf,1);
for n = 1:(N-1)
    IndAncilla(1:n) = IndAncilla(1:n) * SiteDim;
end
Basis = zeros(Nconf,N);
for b = 1:Nconf
    Conf = mod(floor((b-1)./IndAncilla),SiteDim) - S;
    Basis(b,:) = Conf*2; % For ease of use, changing from +/- 1/2 to +/- 1.
    if sum(Conf) == 0 % Recording sites with total spin zero.
        S0Inds(b) = 1;
    end
end

% Constructing the Hamiltonian matrix using Kronecker products.
sz = sparse([1 0; 0 -1]); SZ = cell(L,1);
sx = sparse([0 1; 1 0]); SX = cell(L,1);
HFull = sparse(Nconf,Nconf); J = -1; Bx = 2;
% Construct operator matrices.
for i=1:L
    Dim1 = SiteDim^(i-1);
    Dim2 = SiteDim^(L-i);
    SZ{i}=kron(kron(spdiags(ones(Dim1,1),0,Dim1,Dim1),sz),spdiags(ones(Dim2,1),0,Dim2,Dim2));
    SX{i}=kron(kron(spdiags(ones(Dim1,1),0,Dim1,Dim1),sx),spdiags(ones(Dim2,1),0,Dim2,Dim2));
end
for i = 1:4
    for j = 1:4
    ip = 1 + mod(i,4); jp = 1 + mod(j,4);
    s0 = i + (j-1)*4; sv = ip + (j-1)*4; sh = i + (jp-1)*4;
    HFull = HFull + (0.25*J*(SZ{s0}*(SZ{sh} + SZ{sv}))+Bx*0.5*SX{s0})/N;
    end
end

[GS,E0] = eigs(HFull,1,'sa'); % GS will be a Nconf x 1 vector, E0 is the energy.
% Note that this particular construction assumes a configuration ordering
% like I've generated above to get Basis, i.e. only use wavefunctions
% generated from this Basis to compare with this HFull, or your results
% will be completely borked.


% Using the basis generated here, write and use a PsiGenerate function for
% the CPS and compare its amplitudes with GS, or the fidelity F =
% abs(GS'*(your wavefunction))^2.

% This will give the expectation value of the transverse spin projection
% per site - replace GS with your wavefunction from PsiGenerate to compare.
SxExpVal = 0;
for i = 1:N
    SxExpVal = SxExpVal + S*GS'*SX{i}*GS/N; 
end
% For setting up the Hamiltonian, use the Hilbert, Graph and Operator
% objects here:

HilbertObj = Spin(N,S,[]);

GraphObj = HypCub(L,'PBC',1,1);

SzSzOp = OperatorDg(HilbertObj,GraphObj,@SzSz_CfgVal);
SxOp = Operator1S(HilbertObj,GraphObj,@Sx_OpMatEls);

HamiltonianObj = Hamiltonian({SzSzOp; SxOp},[J; Bx]);