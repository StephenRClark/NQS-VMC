% Script for constructing and optimising Ansatz objects with Jastrow
% Modifiers and Bose condensate References for the Bose Hubbard model on a
% 2D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
N = 10; Dim = N;

%% Specify Graph object.
LVecs = 1; SFlag = 1;
% SFlag determines whether to generate all possible separations in
% GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours as well as
% for multiple nearest neighbour separations.

%% Specify Hilbert object.
Nb = N; Nmax = 3; NStr = ' N ';
HilbertObj = Bose(N,Nb,Nmax);
% Specify Hilbert for NQS hidden units.
HilbertNQS = [Bose(N,Nb,1); Bose(N,Nb,2); Bose(N,Nb,3); Bose(N,Nb,4)];

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [2 6 10];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t; U(1)/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj = SetNsamp(SamplerObj,8000);

%% Initialise Stochastic Reconfig Optimiser.
Npass1 = 300; Npass2 = 200; Ncore = 4;

SR1 = StochasticReconfig(Npass1,Ncore);
SR2 = StochasticReconfig(Npass2,Ncore);

SR1 = SetExtraSamples(SR1,8000);
SR2 = SetExtraSamples(SR2,8000);

SR1 = SetSRTolerances(SR1,5+2*U(1),1e-6);
SR2 = SetSRTolerances(SR2,5+2*U(1),1e-6);

SR1 = SetRegularisation(SR1,1e4,1e1,0.92);
SR2 = SetRegularisation(SR2,1e1,1e-2,0.96);

SR2 = SetLearnRate(SR2,0.5);

%% Initialise starting Ansatz.

% Initialise Reference and Modifier(s) to slot into Ansatz.

% Chosen Reference: BECR Necessary fields: Hilbert, Graph, Params.SPH

% Can construct a single particle Hamiltonian (hopping only) using
% Graph2Array.
CArr = Graph2Array(GraphObj,0); RefParams.SPH = -CArr;

Ref = BECR(HilbertObj,GraphObj,RefParams);

% Chosen Modifier: NQSNHTI Necessary fields: Hilbert, Graph, Params.(a, b,
% c, nmag, nphs, A, B, Nh/Alpha), VFlag (set to 1)

Alpha = 1;
ModParams.Nh = Alpha * N; ModParams.a = 0; ModParams.b = -0.01; ModParams.W = -0.01;
ModParams.nmag = 0.05; ModParams.nphs = 0; ModParams.A = -0.01; ModParams.B = 0;

for d = 1:numel(HilbertNQS)
    
    NQSObj = NQSNHTI(HilbertNQS(d),GraphObj,ModParams,1);
    
    NQSObj.OptInds([1 3 4]) = 0; % Fixes a, b, B.
    
    % Initialise Ansatz object. Mods need to be passed as a cell list.
    Mod = {NQSObj};
    
    AnsatzObj = Ansatz(Ref,Mod,HilbertObj);
    
    AnsStr = ['BECR-NQSNHTI-HD' num2str(NQSObj.HDim) ' Alpha ' num2str(Alpha)];
    % Ansatz name string for easy identification later.
    
    mkdir([AnsStr NStr num2str(N)]);
    
    for u = 1:numel(U)
        if u > 1
            % Alter Hamiltonian parameter corresponding to interaction.
            HamiltonianObj.HParams(2) = U(u)/2;
            SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
            % Cycle out old NQS Modifier with unoptimised one.
            AnsatzObj = ModReplace(AnsatzObj,NQSObj,1);
            % Reset energy tolerance in SR to avoid getting stuck early in
            % optimisation.
            SR1 = SetSRTolerances(SR1,5+2*U(u),1e-6);
            SR2 = SetSRTolerances(SR2,5+2*U(u),1e-6);
        end
        tic;
        % Perform optimisation of AnsatzObj with SR.
        EnIter = cell(2,1);
        [AnsatzObj,EnIter{1}] = SR1.Optimise(SamplerObj,AnsatzObj);
        [AnsatzObj,EnIter{2}] = SR2.Optimise(SamplerObj,AnsatzObj);
        RunTime = toc; EnIter = cell2mat(EnIter);
        save([AnsStr NStr num2str(N) '/BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
            'AnsatzObj','EnIter','RunTime');
        % Save AnsatzObj and run details for later analysis.
    end
    
end

%% Fine tune optimisation.

HamiltonianObj.HParams(1) = U(1)/2;

SamplerObjFT = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObjFT = SetNsamp(SamplerObjFT,16000);

SRFT1 = StochasticReconfig(Npass1,Ncore);
SRFT2 = StochasticReconfig(Npass2,Ncore);

SRFT1 = SetExtraSamples(SRFT1,16000);
SRFT2 = SetExtraSamples(SRFT2,16000);

SRFT1 = SetSRTolerances(SRFT1,5+2*U(1),1e-6);
SRFT2 = SetSRTolerances(SRFT2,5+2*U(1),1e-6);

SRFT1 = SetRegularisation(SRFT1,1e4,1e-3,0.92);
SRFT2 = SetRegularisation(SRFT2,1e-3,1e-3,0.96);

SRFT1 = SetLearnRate(SRFT1,0.1);
SRFT2 = SetLearnRate(SRFT2,0.05);

for d = 1:numel(HilbertNQS)
    AnsStr = ['BECR-NQSNHTI-HD' num2str(NQSObj.HDim) ' Alpha ' num2str(Alpha)];
    for u = 1:numel(U)
        load([AnsStr NStr num2str(N) '/BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');
        if u > 1
            % Alter Hamiltonian parameter corresponding to interaction.
            HamiltonianObj.HParams(2) = U(u)/2;
            SamplerObjFT = SetHamiltonian(SamplerObjFT,HamiltonianObj);
            % Reset energy tolerance in SR to avoid getting stuck early in
            % optimisation.
            SRFT1 = SetSRTolerances(SRFT1,5+2*U(u),1e-6);
            SRFT2 = SetSRTolerances(SRFT2,5+2*U(u),1e-6);
        end
        tic;
        % Perform optimisation of AnsatzObj with SR.
        EnIter = cell(2,1);
        [AnsatzObj,EnIter{1}] = SRFT1.Optimise(SamplerObjFT,AnsatzObj);
        [AnsatzObj,EnIter{2}] = SRFT2.Optimise(SamplerObjFT,AnsatzObj);
        RunTime = toc; EnIter = cell2mat(EnIter);
        save([AnsStr NStr num2str(N) '/BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' FT Logs.mat'],...
            'AnsatzObj','EnIter','RunTime');
        % Save AnsatzObj and run details for later analysis.
    end
end