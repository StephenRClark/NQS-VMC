% Script for constructing and optimising Ansatz objects for the spin-1/2
% XXZ model in 2D using NQS.

%% Specify the basic information for the lattice.
L = 4; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; Bound = [1 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,Bound,LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

%% Specify Hilbert object.
S = 1/2; Sector = 0; % Set the total Sz to be fixed at zero.
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.
% dRXX = [1 0; 0 1]; dRZZ = [1 0; 0 1];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
% HGraphs = [HypCub(Dim,Bound,dRXX,1); HypCub(Dim,Bound,dRZZ,1)];

HGraph = HypCub(Dim,Bound,LVecs,1);
for i = 1:L
    for j = 1:L
        HGraph.SLInds(i+(j-1)*L) = 1+mod(i+j,2);
    end
end

% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
J = 1;
% Jx = 1; Jz = 1; % [0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5...
% 0.6 0.7 0.75 0.8 0.85 0.9 0.925 0.95 0.97 0.98 0.99 1];

SiSjOp = Operator2S(HilbertObj,HGraph,@SiSj_GT_OpMatEls);

% Specify the Operators that will be used in the Hamiltonian.
% SzSzOp = OperatorDg(HilbertObj,HGraphs(1),@SzSz_CfgVal);
% SzSzOp describes the two-site Sz-Sz operator in the XXZ model. It is
% diagonal in the configuration basis, hence the Dg subclass.
% SpSmOp = Operator2S(HilbertObj,HGraphs(2),@SpmSmp_GT_OpMatEls);
% GSzOp describes the transverse field S+S- + h.c. operator in the XXZ
% model. It is not diagonal in the configuration basis and acts on pairs of
% sites, hence the 2S subclass.
% Operators = {SzSzOp; SpSmOp}; HParams = [Jz(1), Jx];

Operators = {SiSjOp}; HParams = [J];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj = SetNsamp(SamplerObj,4000);

%% Initialise Stochastic Reconfig Optimiser.
Npass = 1000; Ncore = 1;

SR = StochasticReconfig(Npass,Ncore);

SR = SetExtraSamples(SR,4000);

SR = SetSRTolerances(SR,3,1e-6);

SR = SetRegularisation(SR,1e3,1e-3,exp(-6*log(10)/500));

SR = SetLearnRate(SR,1);

%% Initialise starting Ansatz.

% Initialise Reference and Modifier(s) to slot into Ansatz.

% Chosen Reference: Plus (fixed)
% Necessary fields: N/A

Ref = Plus();

% Chosen Modifier: NQSTI
% Necessary fields: Hilbert, Graph, Params.(a, b, c, nmag, nphs, Nh/Alpha),
% VFlag (set to 1)

Alpha = 4;
ModParams.Nh = Alpha * N; ModParams.a = 0.01; ModParams.b = 0.01; ModParams.c = -0.01;
ModParams.W = -0.01; ModParams.nmag = 0.05; ModParams.nphs = 0;

NQSObj = NQSTI(HilbertObj,GraphObj,ModParams,1);

% Initialise Ansatz object. Mods need to be passed as a cell list.
Mod = {NQSObj};

AnsatzObj = Ansatz(Ref,Mod,HilbertObj);

AnsStr = ['Plus-NQSTI Alpha ' num2str(Alpha)]; % Ansatz name string for easy identification later.

for z = 1 % :numel(Jz)
    % Alter Hamiltonian parameter corresponding to interaction.
    % HamiltonianObj.HParams(1) = Jz(z);
    % SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    % Cycle out optimised NQS Modifier with unoptimised one to avoid
    % biasing.
    AnsatzObj = ModReplace(AnsatzObj,NQSObj,1);
    tic;
    % Perform optimisation of AnsatzObj with SR.
    [AnsatzObj,EnIter] = SR.Optimise(SamplerObj,AnsatzObj);
    RunTime = toc;
    save('Heisenberg 2D N 16 example optimisation.mat','AnsatzObj','EnIter','RunTime');
    % save(['XXZ 2D Jz ' num2str(Jz(z)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],...
    %    'AnsatzObj','EnIter','RunTime');
    % Save AnsatzObj and run details for later analysis.
end