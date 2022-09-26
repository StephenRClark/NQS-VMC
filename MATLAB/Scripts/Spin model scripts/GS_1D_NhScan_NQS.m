% Script for constructing and optimising Ansatz objects for a spin-1/2
% graph state Hamiltonian in 1D.

%% Specify the basic information for the lattice.
N = 12; Dim = N; Nt = N/3;

%% Specify Graph object.
LVecs = 1; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,1,LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

% Set site groupings in GraphObj.ExtraLabels:
TripletList = zeros(Nt,3);
LeftQuadList = zeros(Nt,4);
RightQuadList = zeros(Nt,4);
for t = 1:Nt
    TripletList(t,:) = [1 2 3] + (t-1)*3;
    LeftQuadList(t,:) = mod([-1 0 1 2] + (t-1)*3,N) + 1;
    RightQuadList(t,:) = mod([0 1 2 3] + (t-1)*3,N) + 1;
end
GraphObj.ExtraLabels = {TripletList; LeftQuadList; RightQuadList};

%% Specify Hilbert object.
S = 1/2; Sector = []; % Set the total Sz to be fixed at zero.
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.

% Specify the Operators that will be used in the Hamiltonian.
XYXOp = OperatorMS(HilbertObj,GraphObj,@SMS_XYX_OpMatEls,1);
% Triplet term - XYX.

LXZZXOp = OperatorMS(HilbertObj,GraphObj,@SMS_XZZX_OpMatEls,2);
% Left quad term - XZZX.

RXZZXOp = OperatorMS(HilbertObj,GraphObj,@SMS_XZZX_OpMatEls,3);
% Right quad term - XZZX.

Operators = {XYXOp; LXZZXOp; RXZZXOp}; HParams = [-1, -1, -1];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj = SetNsamp(SamplerObj,3200);

%% Initialise Stochastic Reconfig Optimiser.
Npass = 250; Ncore = 4;

SR = StochasticReconfig(Npass,Ncore);

SR = SetExtraSamples(SR,4800);

SR = SetSRTolerances(SR,3,1e-6);

SR = SetRegularisation(SR,1e4,1e-3,0.9);

%% Initialise starting Ansatz.

% Initialise Reference and Modifier(s) to slot into Ansatz.

% Chosen Reference: Plus (fixed)
% Necessary fields: N/A

Ref = Plus();

% Chosen Modifier: NQSTI
% Necessary fields: Hilbert, Graph, Params.(a, b, c, nmag, nphs, Nh/Alpha),
% VFlag (set to 1)

ModParams.Nh = 1; ModParams.a = 0.01; ModParams.b = 0.01; ModParams.c = -0.01;
ModParams.nmag = 0.05; ModParams.nphs = 0;

NQSObj = NQS(HilbertObj,GraphObj,ModParams,1);

% Initialise Ansatz object. Mods need to be passed as a cell list.
Mod = {NQSObj};

AnsatzObj = Ansatz(Ref,Mod,HilbertObj);

AnsStr = ['Plus-NQS Nh ' num2str(Nh)]; % Ansatz name string for easy identification later.

for z = 1:numel(Jz)
    if z > 1
        % Alter Hamiltonian parameter corresponding to interaction.
        HamiltonianObj.HParams(1) = Jz(z);
        SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
        % Cycle out optimised NQS Modifier with unoptimised one to avoid
        % biasing.
        AnsatzObj = ModReplace(AnsatzObj,NQSObj,1);
    end
    tic;
    % Perform optimisation of AnsatzObj with SR.
    [AnsatzObj,EnIter] = SR.Optimise(SamplerObj,AnsatzObj);
    RunTime = toc;
    save(['XXZ 1D Jz ' num2str(Jz(z)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','EnIter','RunTime');
    % Save AnsatzObj and run details for later analysis.
end