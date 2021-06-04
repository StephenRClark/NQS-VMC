% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators in one spatial dimension.

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
Nb = N; % Nb specifies number of bosons, 
Nmax = 4; % Nmax specifies maximum on-site occupation.
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% dRT sets separation (1) of site pairs for use with hopping operator.
% dRU = 0 specifies single site n^2 correlations for interaction operator.

% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,Bound,dRT,1); HypCub(Dim,Bound,dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [2 4 6 8 10]; HParams = [t; U(1)/2];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2).
Operators = {HopOp; IntOp};

HamiltonianObj = Hamiltonian(Operators,HParams);

SysStr = 'BHM 1D U '; % Identifier for the Hamiltonian / system.

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

SamplerObj = Sampler(HilbertObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,10000); Ncore = 1;

% Ansatz identifier string - used to select file containing relevant Ansatz.
AnsStr = 'BECR-Jast';

urange = 1:numel(U); % Specify which values of U to sample for.

for u = urange
    HamiltonianObj.HParams(2) = U(u)/2;
    SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    load([SysStr ' U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    tic;
    [EnAvg,EnSamp,EvalAvg,~] = EvalSample(SamplerObj,AnsatzObj);
    EvalTime = toc;   
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
    % Save observable values in file containing Ansatz and run logs.
    EneGS = EnAvg/N; VarE = mean((EnSamp(:)/N - EnAvg/N).^2); VarN = EvalAvg{1}; 
    NiNj = reshape(EvalAvg{2},N,1); DbHl = reshape(EvalAvg{3},N,1); 
    OcFr = reshape(EvalAvg{4},Nmax+1,1); BiBj = EvalAvg{5}; RunDate = date;
    save([SysStr num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','Params','EvalTime','RunDate','VarE',...
        'VarN','NiNj','DbHl','OcFr','BiBj');
end