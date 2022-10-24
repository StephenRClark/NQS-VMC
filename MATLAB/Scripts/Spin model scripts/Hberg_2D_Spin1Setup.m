L = 3; N = L^2; Dim = [L L]; S = 1; SiteDim = 2*S + 1;

% This section will generate all the configurations possible for the
% specified lattice.
IndAncilla = ones(1,N); Nconf = SiteDim^N; S0Inds = zeros(Nconf,1);
for n = 1:(N-1)
    IndAncilla(1:n) = IndAncilla(1:n) * SiteDim;
end
Basis = zeros(Nconf,N); 
for b = 1:Nconf
    Conf = mod(floor((b-1)./IndAncilla),SiteDim) - S;
    Basis(b,:) = Conf;
    if sum(Conf) == 0 % Recording sites with total spin zero.
        S0Inds(b) = 1;
    end
end

% Constructing the Hamiltonian matrix using Kronecker products.
sz = sparse([1 0 0; 0 0 0; 0 0 -1]); SZ = cell(L,L);
sp = sparse([0 0 0; sqrt(2) 0 0 ; 0 sqrt(2) 0]); SP = cell(L,L);
sm = sparse([0 sqrt(2) 0 ;0 0 sqrt(2); 0 0 0]); SM = cell(L,L);
HFull = sparse(Nconf,Nconf); J = 1;
% Construct operator matrices.
for i=1:Dim(1)
    for j=1:Dim(2)
        m = (i-1)*Dim(1)+j;
        Dim1 = SiteDim^(m-1);
        Dim2 = SiteDim^(Dim(1)*Dim(2)-m);
        SZ{i,j}=kron(kron(spdiags(ones(Dim1,1),0,Dim1,Dim1),sz),spdiags(ones(Dim2,1),0,Dim2,Dim2));
        SP{i,j}=kron(kron(spdiags(ones(Dim1,1),0,Dim1,Dim1),sp),spdiags(ones(Dim2,1),0,Dim2,Dim2));
        SM{i,j}=kron(kron(spdiags(ones(Dim1,1),0,Dim1,Dim1),sm),spdiags(ones(Dim2,1),0,Dim2,Dim2));
    end
end
for i = 1:Dim(1)
    for j = 1:Dim(2)
        ip = 1 + mod(i,Dim(1)); jp = 1 + mod(j,Dim(2));
        HFull = HFull + J*(SZ{i,j}*(SZ{ip,j}+SZ{i,jp})+0.5*(SP{i,j}*(SM{i,jp}+SM{ip,j})+SM{i,j}*(SP{i,jp}+SP{ip,j})))/N;
    end
end
HRed = HFull(S0Inds==1,S0Inds==1); % Reduced Hamiltonian matrix that only concerns S = 0 configurations.

[GS,E0] = eigs(HFull,1,'sa'); % GS will be a Nconf x 1 vector, E0 is the energy.

% Using the basis generated here, write and use a PsiGenerate function for
% the CPS and compare its amplitudes with GS.

% For setting up the Hamiltonian, use the Hilbert, Graph and Operator
% objects here:

HilbertObj = Spin(N,S,0);

GraphObj = HypCub(Dim,'PBC',eye(2),1);

OperatorObj = Operator2S(HilbertObj,GraphObj,@SiSj_S1_OpMatEls);