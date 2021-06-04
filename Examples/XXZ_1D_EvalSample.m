% Script for sampling pre-optimised Ansatz objects with the spin-1/2 XXZ
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

% Construct the Operator objects for sampling and place in a cell array. In
% this case, operators are two-site Sz-Sz correlations and two-site full
% spin-spin correlations.
SampOperators = {OperatorDg(HilbertObj,GraphObj,@SzSz_CfgVal),...
    Operator2S(HilbertObj,GraphObj,@SiSj_OpMatEls)};

SamplerObj = Sampler(HilbertObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,10000); Ncore = 1;

% Ansatz identifier string - used to select file containing relevant Ansatz.
AnsStr = 'Plus-NQSTI Alpha 1';

jrange = 1:numel(JzV); % Specify which values of Jz to sample for.

for j = jrange
    HamiltonianObj.HParams(1) = JzV(j);
    SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    load([SysStr ' Jz ' num2str(JzV(j)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    tic;
    [EnAvg,EnSamp,EvalAvg,~] = MultiChainSample(EvalSampler,AnsatzObj,Ncore);
    EvalTime = toc;   
    disp(['Sampling for Jz = ' num2str(U(j)) ' complete.']);
    % Save observable values in file containing Ansatz and run logs.
    EneGS = EnAvg/N; VarE = mean((EnSamp(:)/N - EnAvg/N).^2); 
    SzSz = reshape(EvalAvg{1},N,1); SiSj = reshape(EvalAvg{2},N,1); 
    save([SysStr num2str(JzV(j)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','Params','EvalTime','RunDate','VarE',...
        'SzSz','SiSj');
end