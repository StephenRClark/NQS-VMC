% Script for constructing and optimising Ansatz objects with NQS
% Modifiers and Fermi sea References for the Fermi Hubbard model on a
% 1D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
N = 20; Dim = N; 

%% Specify Graph object.
Bound = 1; LVecs = 1; SFlag = 1; 
% Bound = 1 specifies periodic boundary conditions.
% LVecs correspond to primitive lattice vectors for the lattice - +1 here.
% SFlag = 1 tells the Graph to generate all mappings through all distinct
% combinations of provided primitive lattice vectors.
GraphObj = HypCub(Dim,Bound,LVecs,SFlag);

%% Specify Hilbert object.
Nferm = N; % Set number of fermions - leave as empty if not fixed.
SzTotal = 0; % Set total Sz projection, leave as empty if not fixed.
HilbertObj = Ferm(N,Nferm,SzTotal);

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% dRT sets separation (1) of site pairs for use with hopping operator.
% dRU = 0 specifies single site density correlations for interaction operator.

% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,Bound,dRT,1); HypCub(Dim,Bound,dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [2 4 6 8 10]; HParams = [t; U(1)/2];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@CpCm_OpMatEls);
% HopOp describes the fermionic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Ferm_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2).
Operators = {HopOp; IntOp};

HamiltonianObj = Hamiltonian(Operators,HParams);

SysStr = 'FHM 1D';

%% Initialise Sampler.

% Setting up three Sampler objects of differing sample numbers.
SamplerObj1 = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj1 = SetNsamp(SamplerObj1,2000);
SamplerObj1 = SetNblock(SamplerObj1,N); SamplerObj1 = SetNequil(SamplerObj1,1000);

SamplerObj2 = SetNsamp(SamplerObj1,3000);
SamplerObj2 = SetNblock(SamplerObj2,N); SamplerObj2 = SetNequil(SamplerObj2,1000);

SamplerObj3 = SetNsamp(SamplerObj2,4000);
SamplerObj3 = SetNblock(SamplerObj3,N); SamplerObj3 = SetNequil(SamplerObj3,1000);

% Set up evaluation Sampler for after optimisation.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@DbOc_Ferm_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@NiNj_Ferm_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@Rgdy_Ferm_CfgVal)};

EvalSampler = Sampler(HilbertObj,HamiltonianObj,SampOperators);

EvalSampler = SetNsamp(EvalSampler,10000);

%% Initialise Stochastic Reconfig Optimiser.
% Setting up a three stage optimisation process with separate SR objects.
Ncore = 1; % Single threading version - set to desired number of cores.

Npass1 = 300; Npass2 = 400; Npass3 = 300;
dEV1 = 0.2/U(1); dEV2 = 0.1/U(1); dEV3 = 0.05/U(1);
% Set Npass, number of optimisation steps.
SR1 = StochasticReconfig(Npass1,Ncore);
SR2 = StochasticReconfig(Npass2,Ncore);
SR3 = StochasticReconfig(Npass3,Ncore);
% Set ExtraSamples if bad energy sample occurs.
SR1 = SetExtraSamples(SR1,2000);
SR2 = SetExtraSamples(SR2,3000);
SR3 = SetExtraSamples(SR3,4000);
% Set bounds to detect unphysical energy.
SR1 = SetEnergyTolerances(SR1,dEV1,0.05);
SR2 = SetEnergyTolerances(SR2,dEV2,0.05);
SR3 = SetEnergyTolerances(SR3,dEV3,0.05);
% Optimisation stability hyperparameters for SR1.
RegMax1 = [1e5 1e5 1e5 1e5 1e5];
RegMin1 = [1e3 1e3 1e3 1e3 1e3];
LearnRate1 = [1 1 1 1 1];
% Optimisation stability hyperparameters for SR2.
RegMax2 = [1e3 1e3 1e3 1e3 1e3];
RegMin2 = [1e1 1e1 1e1 1e1 1e1];
LearnRate2 = [0.1 0.1 0.1 0.1 0.1];
% Optimisation stability hyperparameters for SR3.
RegMax3 = [1e1 1e1 1e1 1e1 1e1];
RegMin3 = [1e-1 1e-1 1e-1 1e-1 1e-1];
LearnRate3 = [0.01 0.01 0.01 0.01 0.01];
% Regularisation decay factors and number of runs at regularisation floor.
Nfloor = 50; 
RegFac1 = exp(-2*log(10)/(Npass1-Nfloor)); 
RegFac2 = exp(-2*log(10)/(Npass2-Nfloor)); 
RegFac3 = exp(-2*log(10)/(Npass3-Nfloor));

%% Initialise starting Ansatz.

% Construct the BECR (Bose Einstein Condensate Reference) Reference object.
% Necessary fields: Hilbert, Graph, Params.[CArr, HVar]
CArr = Graph2Array(GraphObj,1); RefParams.CArr = CArr; RefParams.HVar = -1;
% Single particle Hamiltonian constructed inside SDet reference using CArr
% and HVar fields.
Ref = SDet(HilbertObj,GraphObj,RefParams);

% Construct the NQSTI Modifier object.
% NQSTI - translation invariant NQS
% Necessary fields: Hilbert, Graph, Params.[a, b, W, Alpha or Nh], VFlag (set to 1)
ModParams.a = 0.1; ModParams.b = 0.1; ModParams.W = 0.1; ModParams.Alpha = 1;
% Random noise for starting parameters in magnitude and phase.
ModParams.nmag = 0.1; ModParams.nphs = 0;
NQSObj = NQSTI(HilbertObj,GraphObj,ModParams,1);
% Initialise Ansatz object. Mods need to be passed as a cell list.
Mod = {NQSObj};

% Initialise Ansatz object.
AnsatzObj = Ansatz(Ref,Mod,HilbertObj);

% Set Ansatz identifier - convention is Ref-Mod1-Mod2...
AnsStr = 'SDet-NQSTI Alpha 1';

% Set values of U to sample with urange.
urange = 1:numel(U);

% Set up directory for results to be saved into.
DirStr = [SysStr ' N ' num2str(N) ' ' AnsStr]; mkdir(DirStr);

for u = urange
    % Replace Modifier after each optimisation to avoid biasing.
    AnsatzObj = AnsatzObj.ModReplace(NQSObj,1);
    % SaveStr is the filename the Ansatz and sampled observables will be
    % saved to, within directory specified by DirStr.
    SaveStr = [SysStr ' U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'];
    % Adjust Hamiltonian objects and reload into Samplers.
    HamiltonianObj.HParams(2) = U(u)/2;
    EvalHamiltonian.HParams(2) = U(u)/2;
    SamplerObj1 = SetHamiltonian(SamplerObj1,HamiltonianObj);
    SamplerObj2 = SetHamiltonian(SamplerObj2,HamiltonianObj);
    SamplerObj3 = SetHamiltonian(SamplerObj3,HamiltonianObj);
    EvalSampler = SetHamiltonian(EvalSampler,EvalHamiltonian);
    % Set energy tolerances and regularisation for SR objects.
    SR1 = SetSRTolerances(SR1,2+4*U(u),1e-6); SR1 = SetEnergyTolerances(SR1,dEV1,0.05);
    SR2 = SetSRTolerances(SR2,2+4*U(u),1e-6); SR2 = SetEnergyTolerances(SR2,dEV2,0.05);
    SR3 = SetSRTolerances(SR3,2+4*U(u),1e-6); SR3 = SetEnergyTolerances(SR3,dEV3,0.05);
    SR1 = SetLearnRate(SR1,LearnRate1(u)); SR1 = SetRegularisation(SR1,RegMax1(u),RegMin1(u),RegFac1);
    SR2 = SetLearnRate(SR2,LearnRate2(u)); SR2 = SetRegularisation(SR2,RegMax2(u),RegMin2(u),RegFac2);
    SR3 = SetLearnRate(SR3,LearnRate3(u)); SR3 = SetRegularisation(SR3,RegMax3(u),RegMin3(u),RegFac3);
    % Optimisation using each of the three SR objects.
    tic; EnIter = cell(3,1);    
    [AnsatzObj,EnIter{1}] = SR1.Optimise(SamplerObj1,AnsatzObj);
    save([DirStr '/' SaveStr],'AnsatzObj','EnIter');
    [AnsatzObj,EnIter{2}] = SR2.Optimise(SamplerObj2,AnsatzObj);
    save([DirStr '/' SaveStr],'AnsatzObj','EnIter');
    [AnsatzObj,EnIter{3}] = SR3.Optimise(SamplerObj3,AnsatzObj);
    EnIter = cell2mat(EnIter); Params = AnsatzObj.ParamList; RunTime = toc;
    save([DirStr '/' SaveStr],'AnsatzObj','EnIter','RunTime','Params');
    % Perform a sampling run to evaluate observables.
    tic;
    [EnAvg,EnSamp,EvalAvg,~] = MultiChainSample(EvalSampler,AnsatzObj,Ncore);
    EvalTime = toc;   
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
    % Save observable values in file containing Ansatz and run logs.
    EneGS = EnAvg/N; VarE = mean((EnSamp(:)/N - EnAvg/N).^2); 
    DbOc = reshape(EvalAvg{1},N,1); NiNj = reshape(EvalAvg{2},N,1); 
    Rgdy = reshape(EvalAvg{3},N,1); RunDate = date;
    save([SysStr num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','Params','EvalTime','RunDate','VarE',...
        'DbOc','NiNj','Rgdy');
end