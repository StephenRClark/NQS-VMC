% Script for constructing and optimising Ansatz objects with Jastrow
% Modifiers and Bose condensate References for the Bose Hubbard model on a
% 1D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
N = 20; Dim = N;

%% Specify Graph object.
LVecs = 1;
GraphObj = HypCub(Dim,'PBC',LVecs);

%% Specify Hilbert object.
Nb = N; Nmax = 4;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT); HypCub(Dim,'PBC',dRU)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [2 2.1 2.2 2.25 2.3 2.35 2.375 2.4 2.425 2.45 2.475 2.5 2.525 2.55 2.575 2.6 2.7 2.8 2.9 3];
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
Ref = {'BECR',0}; Mod = {'JastTI',1};

% Construct the single particle Hamiltonian necessary for the BECR Reference.
SPH = zeros(N);
Bonds = GraphObj.Bonds;
for n = 1:N
    for b = 1:size(Bonds,2)
        SPH(n,Bonds(n,b)) = -1; SPH(Bonds(n,b),n) = -1;
    end
end
Params.SPH = SPH;
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
    save(['BHM U ' num2str(U(u)) ' N ' num2str(N) ' BECR-JastTI Logs.mat'],'AnsatzObj','EnIter','RunTime');    
end