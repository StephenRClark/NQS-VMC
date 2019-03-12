% Script for constructing and optimising Ansatz objects with Jastrow
% Modifiers and Fermi sea References for the Fermi Hubbard model on a
% 1D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
N = 20; Dim = N;

%% Specify Graph object.
LVecs = 1;
GraphObj = HypCub(Dim,'PBC',LVecs);

%% Specify Hilbert object.
N_up = N/2; N_dn = N/2;
HilbertObj = Ferm(N,N_up,N_dn);

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT); HypCub(Dim,'PBC',dRU)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [0 1/4 1/3 1/2 2/3 1 3/2 2 5/2 3 7/2 4 9/2 5 6 8];
HamiltonianObj = Hubbard(HilbertObj,HGraphs,t,U(1)); 

%% Initialise Sampler.
SamplerObj = Sampler(HilbertObj,GraphObj,HamiltonianObj,{});

SamplerObj = SetNsamp(SamplerObj,20000);

%% Initialise Stochastic Reconfig Optimiser.
Npass = 500; Ncore = 8;
SR = StochasticReconfig(Npass,Ncore);

SR = SetExtraSamples(SR,20000);

SR = SetSRTolerances(SR,3,1e-6);

%% Initialise starting Ansatz.

% Specify desired Reference and Modifier(s), flagging with 1 if
% optimisation is desired.
Ref = {'SDet',0}; Mod = {'JastTI',1};

% Construct the single particle Hamiltonian necessary for the SDet Reference.
CArr = zeros(2*N);
Bonds = GraphObj.Bonds;
for n = 1:N
    for b = 1:size(Bonds,2)
        CArr(n,Bonds(n,b)) = 1; CArr(Bonds(n,b),n) = 1;
        CArr(n+N,Bonds(n,b)+N) = 1; CArr(Bonds(n,b)+N,n+N) = 1;
    end
end
Params.CArr = CArr; Params.HVar = -1;
% Add the fields necessary for the Jastrow Modifier to be initialised.
Params.Js = 0.01; Params.nmag = 0.5;

% Initialise Ansatz object.
AnsatzObj = Ansatz(Ref,Mod,HilbertObj,GraphObj,Params);

for u = 1:numel(U)
    if u > 1
        % Create new Hamiltonian and cycle out old one.
        NewH = Hubbard(HilbertObj,HGraphs,t,U(u));
        SamplerObj = SetHamiltonian(SamplerObj,NewH);
        % Create new Jastrow Modifier and cycle out old one.
        NewJast = JastTI(HilbertObj,GraphObj,Params,1);
        AnsatzObj = ModReplace(AnsatzObj,NewJast,1);
        % Reset energy tolerance in SR to avoid getting stuck early in
        % optimisation.
        SR = SetSRTolerances(SR,3+2*U(u),1e-6);
    end
    tic;
    [AnsatzObj,EnIter] = SR.Optimise(SamplerObj,AnsatzObj);
    RunTime = toc;
    save(['FHM U ' num2str(U(u)) ' N ' num2str(N) ' SDet-JastTI Logs.mat'],'AnsatzObj','EnIter','RunTime');    
end