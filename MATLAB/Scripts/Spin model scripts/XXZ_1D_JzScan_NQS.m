% Script for constructing and optimising Ansatz objects with NQS% Modifiers
% for the XXZ model on a 1D lattice with periodic boundary conditions.

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
S = 1/2; 
Sector = 0; % Sector fixes the total Sz projection, set to empty if variable.
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.

% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
Jx = 1; JzV = [0 0.25 0.5 1 2 4];

% Specify the Operators that will be used in the Hamiltonian.
SzSzOp = OperatorDg(HilbertObj,GraphObj,@SzSz_CfgVal);
% SzSzOp describes the two-site Sz-Sz operator in the XXZ model. It is
% diagonal in the configuration basis, hence the Dg subclass.
SpSmOp = Operator2S(HilbertObj,GraphObj,@SpmSmp_OpMatEls);
% GSzOp describes the transverse field S+S- + h.c operator in the XXZ
% model. It is not diagonal in the configuration basis and acts on pairs of
% sites, hence the 2S subclass.
Operators = {SzSzOp; SpSmOp}; HParams = [JzV(1), Jx/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

SysStr = 'XXZ 1D';

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
SampOperators = {OperatorDg(HilbertObj,GraphObj,@SzSz_CfgVal),...
    Operator2S(HilbertObj,GraphObj,@SiSj_OpMatEls)};

EvalSampler = Sampler(HilbertObj,HamiltonianObj,SampOperators);

EvalSampler = SetNsamp(EvalSampler,10000);

%% Initialise Stochastic Reconfig Optimiser.
% Setting up a three stage optimisation process with separate SR objects.
Ncore = 1; % Single threading version - set to desired number of cores.

Npass1 = 300; Npass2 = 400; Npass3 = 300;
dEV1 = 0.2; dEV2 = 0.1; dEV3 = 0.05;
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
RegMax1 = [1e5 1e5 1e5 1e5 1e5 1e5];
RegMin1 = [1e3 1e3 1e3 1e3 1e3 1e5];
LearnRate1 = [1 1 1 1 1 1];
% Optimisation stability hyperparameters for SR2.
RegMax2 = [1e3 1e3 1e3 1e3 1e3 1e3];
RegMin2 = [1e1 1e1 1e1 1e1 1e1 1e1];
LearnRate2 = [0.1 0.1 0.1 0.1 0.1 0.1];
% Optimisation stability hyperparameters for SR3.
RegMax3 = [1e1 1e1 1e1 1e1 1e1 1e1];
RegMin3 = [1e-1 1e-1 1e-1 1e-1 1e-1 1e-1];
LearnRate3 = [0.01 0.01 0.01 0.01 0.01 0.01];
% Regularisation decay factors and number of runs at regularisation floor.
Nfloor = 50; 
RegFac1 = exp(-2*log(10)/(Npass1-Nfloor)); 
RegFac2 = exp(-2*log(10)/(Npass2-Nfloor)); 
RegFac3 = exp(-2*log(10)/(Npass3-Nfloor));

%% Initialise starting Ansatz.

% Construct the Plus Reference object (equal superposition).
% Necessary fields: N/A
Ref = Plus();

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
AnsStr = 'Plus-NQSTI Alpha 1';

% Set values of Jz to sample with jrange.
jrange = 1:numel(JzV);

% Set up directory for results to be saved into.
DirStr = [SysStr ' N ' num2str(N) ' ' AnsStr]; mkdir(DirStr);

for j = jrange
    % Replace Modifier before each optimisation to avoid biasing.
    AnsatzObj = AnsatzObj.ModReplace(NQSObj,1);
    % SaveStr is the filename the Ansatz and sampled observables will be
    % saved to, within directory specified by DirStr.
    SaveStr = [SysStr ' Jz ' num2str(JzV(j)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'];
    % Adjust Hamiltonian objects and reload into Samplers.
    HamiltonianObj.HParams(1) = JzV(j);
    SamplerObj1 = SetHamiltonian(SamplerObj1,HamiltonianObj);
    SamplerObj2 = SetHamiltonian(SamplerObj2,HamiltonianObj);
    SamplerObj3 = SetHamiltonian(SamplerObj3,HamiltonianObj);
    EvalSampler = SetHamiltonian(EvalSampler,EvalHamiltonian);
    % Set energy tolerances and regularisation for SR objects.
    SR1 = SetSRTolerances(SR1,2+JzV(j),1e-6); SR1 = SetEnergyTolerances(SR1,dEV1,0.05);
    SR2 = SetSRTolerances(SR2,2+JzV(j),1e-6); SR2 = SetEnergyTolerances(SR2,dEV2,0.05);
    SR3 = SetSRTolerances(SR3,2+JzV(j),1e-6); SR3 = SetEnergyTolerances(SR3,dEV3,0.05);
    SR1 = SetLearnRate(SR1,LearnRate1(j)); SR1 = SetRegularisation(SR1,RegMax1(j),RegMin1(j),RegFac1);
    SR2 = SetLearnRate(SR2,LearnRate2(j)); SR2 = SetRegularisation(SR2,RegMax2(j),RegMin2(j),RegFac2);
    SR3 = SetLearnRate(SR3,LearnRate3(j)); SR3 = SetRegularisation(SR3,RegMax3(j),RegMin3(j),RegFac3);
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
    disp(['Sampling for Jz = ' num2str(JzV(j)) ' complete.']);
    % Save observable values in file containing Ansatz and run logs.
    EneGS = EnAvg/N; VarE = mean((EnSamp(:)/N - EnAvg/N).^2); 
    SzSz = reshape(EvalAvg{1},N,1); SiSj = reshape(EvalAvg{2},N,1); 
    save([SysStr num2str(JzV(j)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','Params','EvalTime','RunDate','VarE',...
        'SzSz','SiSj');
end