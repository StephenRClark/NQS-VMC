% Script for constructing and optimising Ansatz objects with Jastrow
% Modifiers and Bose condensate References for the Bose Hubbard model on a
% 2D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
L = 6; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 2D L x L lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours in x and
% y, as well as for multiple nearest neighbour separations.

%% Specify Hilbert object.
% Nb = N-1; NStr = ' N-1 ';
Nb = N; NStr = ' N ';
% Nb = N+1; NStr = ' N+1 ';
Nmax = 4;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = [1 0; 0 1]; dRU = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t; U(1)/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

EvalHamiltonian = Hamiltonian(Operators,HParams);

%% Initialise Sampler.

SamplerObj1 = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj1 = SetNsamp(SamplerObj1,6400); SamplerObj1 = SetNblock(SamplerObj1,L);

SamplerObj2 = SetNsamp(SamplerObj1,9600); SamplerObj2 = SetNblock(SamplerObj2,N);

%% Initialise Sampler.

% Construct the Operator objects for sampling and place in a cell array. In
% this case, operators are variance in site occupation, two-site
% density-density correlations, static structure factor, doublon-holon
% correlations and occupation number fractions/probabilities.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@VarN_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@NiNj_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@DbHl_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@OccFrac_Bose_CfgVal),...
    Operator2S(HilbertObj,HGraphs(1),@BiBj_OpMatEls)};

EvalSampler = Sampler(HilbertObj,HamiltonianObj,SampOperators);

EvalSampler = SetNsamp(EvalSampler,40000);

%% Initialise Stochastic Reconfig Optimiser.
Npass1 = 300; Npass2 = 100; Ncore = 8; LR1 = 0.005*N; LR2 = 0.001*N;

SR1 = StochasticReconfig(Npass1,Ncore);
SR2 = StochasticReconfig(Npass2,Ncore);

SR1 = SetExtraSamples(SR1,6400);
SR2 = SetExtraSamples(SR2,6400);

SR1 = SetSRTolerances(SR1,5+2*U(1),1e-6);
SR2 = SetSRTolerances(SR2,5+2*U(1),1e-6);

SR1 = SetRegularisation(SR1,1e4,1e-3,0.9);
SR2 = SetRegularisation(SR2,1e-3,1e-3,1);

SR1 = SetLearnRate(SR1,LR1);
SR2 = SetLearnRate(SR2,LR2);


%% Initialise starting Ansatz.

AnsStr0 = 'BECR-NQSSHTI-HD5-bB Alpha 1'; % Ansatz name string for easy identification later.
AnsStr = 'BECR-NQSSHTI-HD5-bB Alpha 1 FT';

for u = 1:numel(U)
    load(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr0 ' Logs.mat'],'AnsatzObj');
    if u > 1
        % Alter Hamiltonian parameter corresponding to interaction.
        HamiltonianObj.HParams(1) = t/U(u);
        EvalHamiltonian.HParams(2) = U(u)/2;
        SamplerObj1 = SetHamiltonian(SamplerObj1,HamiltonianObj);
        SamplerObj2 = SetHamiltonian(SamplerObj2,HamiltonianObj);
        EvalSampler = SetHamiltonian(EvalSampler,EvalHamiltonian);
        % Reset energy tolerance in SR to avoid getting stuck early in
        % optimisation.
        SR1 = SetSRTolerances(SR1,5+2*U(u),1e-6);
        SR2 = SetSRTolerances(SR2,5+2*U(u),1e-6);
    end
    tic;
    % Perform optimisation of AnsatzObj with SR.
    EnIter = cell(3,1);
    [AnsatzObj,EnIter{1}] = SR1.Optimise(SamplerObj1,AnsatzObj);
    [AnsatzObj,EnIter{2}] = SR2.Optimise(SamplerObj2,AnsatzObj);
    [AnsatzObj,EnIter{3},Params] = SR2.ParamAvgOptimise(SamplerObj2,AnsatzObj);
    EnIter = cell2mat(EnIter);
    RunTime = toc;    
    save(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','EnIter','RunTime');
    [EnAvg,EnSamp,EvalAvg,~] = MultiChainSample(EvalSampler,AnsatzObj,Ncore);
    % Store sampled quantities.
    EneGS = EnAvg/N; VarN = EvalAvg{1}; NiNj = reshape(EvalAvg{2},L,L);
    DbHl = reshape(EvalAvg{3},L,L); OcFr = reshape(EvalAvg{4},1,Nmax+1);
    VarE = mean((EnSamp(:)/N - EnAvg/N).^2); BiBj = EvalAvg{5}; 
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
    save(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','VarN','NiNj','VarE','DbHl','OcFr','BiBj','Params');
    % Save AnsatzObj and run details for later analysis.
end