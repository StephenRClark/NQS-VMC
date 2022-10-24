% Script for constructing and optimising Ansatz objects for the spin-1/2
% transverse field Ising model in 2D using NQS.

%% Specify the basic information for the lattice.
L = 4; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

%% Specify Hilbert object.
S = 1/2; Sector = []; % Allow the total Sz of the system to change freely.
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.
dRJ = [1 0; 0 1]; dRG = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRJ,1); HypCub(Dim,'PBC',dRG,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
J = -1; G = [0 1/2 1 3/2 2 3 4 5 6 8 10 12 16 20 25];

% Specify the Operators that will be used in the Hamiltonian.
SzSzOp = OperatorDg(HilbertObj,HGraphs(1),@SzSz_CfgVal);
% SzSzOp describes the two-site Sz-Sz operator in the Ising model. It is
% diagonal in the configuration basis, hence the Dg subclass.
GSzOp = Operator1S(HilbertObj,HGraphs(2),@Sx_OpMatEls);
% GSzOp describes the transverse field Sx operator in the Ising model. It
% is not diagonal in the configuration basis and acts on single sites,
% hence the 1S subclass.
Operators = {SzSzOp; GSzOp}; HParams = [J, J*G(1)];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj = SetNsamp(SamplerObj,3200);

%% Initialise Stochastic Reconfig Optimiser.
Npass = 250; Ncore = 4;

SR = StochasticReconfig(Npass,Ncore);

SR = SetExtraSamples(SR,4800);

SR = SetSRTolerances(SR,5+G(1),1e-6);

SR = SetRegularisation(SR,1e4,1e-3,0.9);

%% Initialise starting Ansatz.

% Initialise Reference and Modifier(s) to slot into Ansatz.

% Chosen Reference: Plus (fixed)
% Necessary fields: N/A

Ref = Plus();

% Chosen Modifier: NQSTI
% Necessary fields: Hilbert, Graph, Params.(a, b, c, nmag, nphs, Nh/Alpha),
% VFlag (set to 1)

Alpha = 1;
ModParams.Nh = Alpha * N; ModParams.a = 0.01; ModParams.b = 0.01; ModParams.c = -0.01;
ModParams.nmag = 0.05; ModParams.nphs = 0;

NQSObj = NQSTI(HilbertObj,GraphObj,ModParams,1);

% Initialise Ansatz object. Mods need to be passed as a cell list.
Mod = {NQSObj};

AnsatzObj = Ansatz(Ref,Mod,HilbertObj);

AnsStr = ['Plus-NQSTI Alpha ' num2str(Alpha)]; % Ansatz name string for easy identification later.

for g = 1:numel(G)
    if g > 1
        % Alter Hamiltonian parameter corresponding to interaction.
        HamiltonianObj.HParams(2) = J*G(g);
        SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
        % Cycle out optimised NQS Modifier with unoptimised one to avoid
        % biasing.
        AnsatzObj = ModReplace(AnsatzObj,NQSObj,1);
        % Reset energy tolerance in SR to avoid getting stuck early in
        % optimisation.
        SR = SetSRTolerances(SR,3+G(g),1e-6);
    end
    tic;
    % Perform optimisation of AnsatzObj with SR.
    [AnsatzObj,EnIter] = SR.Optimise(SamplerObj,AnsatzObj);
    RunTime = toc;
    save(['Ising 2D G ' num2str(G(g)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','EnIter','RunTime');
    % Save AnsatzObj and run details for later analysis.
end