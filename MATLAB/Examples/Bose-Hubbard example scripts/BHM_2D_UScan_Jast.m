% Script for constructing and optimising Ansatz objects with Jastrow
% Modifiers and Bose condensate References for the Bose Hubbard model on a
% 1D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
L = 10; Dim = [L L]; N = prod(Dim);

%% Specify Graph object.
Bound = [1 1]; LVecs = eye(2); SFlag = 1; 
% Bound = 1 specifies periodic boundary conditions.
% LVecs correspond to primitive lattice vectors for the lattice - +1 here.
% SFlag = 1 tells the Graph to generate all mappings through all distinct
% combinations of provided primitive lattice vectors.
GraphObj = HypCub(Dim,Bound,LVecs,SFlag);

%% Specify Hilbert object.
Nb = N; % Nb specifies number of bosons, 
Nmax = 4; % Nmax specifies maximum on-site occupation.
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = eye(2); dRU = 0;
% dRT sets separation (1) of site pairs for use with hopping operator.
% dRU = 0 specifies single site n^2 correlations for interaction operator.

% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,Bound,dRT,1); HypCub(Dim,Bound,dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [16 32]; HParams = [t; 1/2];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2).
Operators = {HopOp; IntOp};

HamiltonianObj = Hamiltonian(Operators,HParams);

SysStr = 'BHM 2D'; % Identifier for the Hamiltonian / system.

%% Initialise Sampler.

% Setting up three Sampler objects of differing sample numbers.
Samp1 = Sampler(HilbertObj,HamiltonianObj,{});

Samp1 = SetNsamp(Samp1,4000);
Samp1 = SetNblock(Samp1,N); Samp1 = SetNequil(Samp1,2500);

Samp2 = SetNsamp(Samp1,6000);
Samp2 = SetNblock(Samp2,N); Samp2 = SetNequil(Samp2,2500);

Samp3 = SetNsamp(Samp2,8000);
Samp3 = SetNblock(Samp3,N); Samp3 = SetNequil(Samp3,5000);

SamplerArray = [Samp1; Samp2; Samp3];

% Set up evaluation Sampler for after optimisation.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@VarN_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@NiNj_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@DbHl_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@OccFrac_Bose_CfgVal)...
    Operator2S(HilbertObj,HGraphs(1),@BiBj_OpMatEls)};

EvalSampler = Sampler(HilbertObj,HamiltonianObj,SampOperators);

EvalSampler = SetNsamp(EvalSampler,40000);

%% Initialise Stochastic Reconfig Optimiser.
% Setting up a three stage optimisation process with separate SR objects.
Ncore = 4; % Single threading version - set to desired number of cores.

Nopt = 3; Npass = [300; 400; 300];
dEV = [0.2; 0.1; 0.05]; dER = [0.2; 0.1; 0.05];
ExtraSamples = [2000; 3000; 4000];
% Set Npass, number of optimisation steps.
SR1 = StochasticReconfig(Npass(1),Ncore);
SR2 = StochasticReconfig(Npass(2),Ncore);
SR3 = StochasticReconfig(Npass(3),Ncore);
SRArray = [SR1; SR2; SR3];

for n = 1:Nopt
    % Set ExtraSamples if bad energy sample occurs.
    SRArray(n) = SRArray(n).SetExtraSamples(ExtraSamples(n));
    % Set bounds to detect unphysical energy.
    SRArray(n) = SRArray(n).SetEnergyTolerances(dEV(n)/U(1),dER(n));
end

% Optimisation stability hyperparameters for SRArray.
RegMax = [1e5 1e5; 1e3 1e3; 1e1 1e1];
RegMin = [1e3 1e3; 1e1 1e1; 1e-1 1e-1];
LearnRate = [10 10; 1 1; 0.1 0.1];
% Regularisation decay factors and number of runs at regularisation floor.
Nfloor = 50; 
RegFac = exp((log(RegMin)-log(RegMax))./(Npass-Nfloor)); 

%% Initialise starting Ansatz.

% Construct the BECR (Bose Einstein Condensate Reference) Reference object.
% Necessary fields: Hilbert, Graph, Params.SPH
CArr = Graph2Array(GraphObj,0); RefParams.SPH = -CArr;
% Construct single particle hopping Hamiltonian from connectivity array 
% generated from Graph.
Ref = BECR(HilbertObj,GraphObj,RefParams);

% Construct the Jastrow Modifier object.
% Necessary fields: Hilbert, Graph, Params.Js, VFlag (set to 1)
ModParams.Js = 0.1; ModParams.nmag = 0.1; ModParams.nphs = 0;
JastObj = Jast(HilbertObj,GraphObj,ModParams,1);
% Initialise Ansatz object. Mods need to be passed as a cell list.
Mod = {JastObj};

% Initialise Ansatz object.
AnsatzObj = Ansatz(Ref,Mod,HilbertObj);

% Set Ansatz identifier - convention is Ref-Mod1-Mod2...
AnsStr = 'BECR-Jast';

% Set values of U to sample with urange.
urange = 1:numel(U);

% Set up directory for results to be saved into.
DirStr = [SysStr ' N ' num2str(N) ' ' AnsStr]; mkdir(DirStr);

for u = urange
    % Replace Modifier before each optimisation to avoid biasing.
    AnsatzObj = AnsatzObj.ModReplace(JastObj,1);
    % SaveStr is the filename the Ansatz and sampled observables will be
    % saved to, within directory specified by DirStr.
    SaveStr = [SysStr ' U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'];
    EnIter = cell(Nopt,1);  
    % Adjust Hamiltonian objects and reload into Samplers.
    HamiltonianObj.HParams(1) = t/U(u);
    EvalHamiltonian.HParams(2) = U(u)/2;
    EvalSampler = SetHamiltonian(EvalSampler,EvalHamiltonian);
    for n = 1:Nopt
        SamplerArray(n) = SetHamiltonian(SamplerArray(n),HamiltonianObj);
        % Set energy tolerances and regularisation for SR objects.
        SRArray(n) = SRArray(n).SetSRTolerances(5/U(u)+4,1e-6); 
        SRArray(n) = SRArray(n).SetEnergyTolerances(dEV(n)/U(u),dER(n));
        SRArray(n) = SRArray(n).SetLearnRate(LearnRate(n,u));
        SRArray(n) = SRArray(n).SetRegularisation(RegMax(n,u),RegMin(n,u),RegFac(n,u));
    end
    % Optimisation using each of the SR objects.
    tic;
    for n = 1:Nopt
        [AnsatzObj,EnIter{n}] = SRArray(n).Optimise(SamplerArray(n),AnsatzObj);
        save([DirStr '/' SaveStr],'AnsatzObj','EnIter','n','SRArray','SamplerArray');
    end
    EnIter = cell2mat(EnIter); Params = AnsatzObj.ParamList; RunTime = toc;
    save([DirStr '/' SaveStr],'AnsatzObj','EnIter','RunTime','Params');
    % Perform a sampling run to evaluate observables.
    tic;
    [EnAvg,EnSamp,EvalAvg,~] = MultiChainSample(EvalSampler,AnsatzObj,Ncore);
    EvalTime = toc; EneGS = EnAvg/N;VarE = mean(EnSamp(:).^2 - mean(EnSamp(:))^2)/N;
    VarN = EvalAvg{1}; NiNj = reshape(EvalAvg{2},L,L); DbHl = reshape(EvalAvg{3},L,L); 
    OcFr = reshape(EvalAvg{4},1,Nmax+1); BiBj = EvalAvg{5}; RunDate = datetime("today");
    % Save AnsatzObj and run details for later analysis.
    save([DirStr '/' SaveStr],'AnsatzObj','RunTime','EnIter','EvalTime','RunDate',...
        'Params','EneGS','VarE','NiNj','VarN','DbHl','OcFr','BiBj'); 
end