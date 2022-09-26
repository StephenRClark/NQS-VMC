% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

%% Add the folders containing class definitions to MatLab's path.
addpath(genpath('Variational Monte Carlo'));

%% Specify the basic information for the lattice.
N = 40; Dim = N;

%% Specify Graph object.
LVecs = 1; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.
for n = 1:N
    GraphObj.SLInds(n) = 1+mod(n-1,2);
end

%% Specify Hilbert object.
S = 1/2; Sector = 0;
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.
dR = [1 0; 0 1]; J = 1;
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.

% Specify the Operators that will be used in the Hamiltonian.
SiSjOp = Operator2S(HilbertObj,GraphObj,@SiSj_GT_OpMatEls);

Operators = {SiSjOp}; HParams = J;

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.

SpSmQ = Operator2S(HilbertObj,GraphObj,@SpSmQ_GT_OpMatEls);

OpOrderPH1 = {'L',1,[1; 1; 3]; 'R',1,[2; 2; 3]}; % Operator ordering input.
% First row specifies left operator, 1st from configuration, first index is
% spatial separation, 3rd is vector output from matrix elements. Will
% output an Nx1xN array on its own.

% Second row specifies right operator, 1st from configuration, second index
% is spatial separation, 3rd is vector output from matrix elements. Will
% output an 1xNxN array on its own.

OpOrderPH2 = {'L',1,[1; 1; 3]; 'R',2,[2; 2; 3]}; HOrder = {[],1}; % Operator ordering input.
% HOrder specifies H as 1st from configuration.

ExOp = OperatorC(HilbertObj,{SpSmQ SpSmQ},OpOrderPH1); % Output is NxNxN.
HExOp = OperatorH(HilbertObj,HamiltonianObj,{SpSmQ SpSmQ},OpOrderPH2,HOrder);

SampOperators = {ExOp; HExOp};

SamplerObj = Sampler(HilbertObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,100); Ncore = 16;

AnsStr = 'Plus-NQSTI Alpha 4'; % Ansatz identifier string.

load(['Heisenberg 1D N ' num2str(N) ' ' AnsStr ' Logs.mat']);
% Initiate Monte Carlo sampling to evaluate energy and other
% observables.
[EnAvg,EnSamp,EvalAvg,~] = EvalSample(SamplerObj,AnsatzObj);
% Store sampled quantities.
EneGS = EnAvg/N; VarE = mean((EnSamp(:)/N - EnAvg/N).^2); % Local per-site energy variance.
ExOpEval = reshape(EvalAvg{1},N,N,N); HExOpEval = reshape(EvalAvg{2},N,N,N);
    
save(['BHM 1D N ' num2str(N) ' U Scan ' AnsStr ' excitation observable values.mat'],...
    'EneGS','VarE','ExOpEval','HExOpEval');

exit