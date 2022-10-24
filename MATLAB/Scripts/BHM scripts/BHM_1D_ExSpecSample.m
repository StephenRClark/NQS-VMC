% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

%% Specify the basic information for the lattice.
N = 10; Dim = N;

%% Specify Graph object.
LVecs = 1; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

%% Specify Hilbert object.
Nb = N; Nmax = 3;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = 1; dRU =  0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 25/4 20/3 15/2 35/4 10 12 16 20 30 50 100];
% Reentrance U vector: [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 25/4 20/3 15/2 35/4 10 12 16 20 30 50 100]
% Exact test U vector: [0 1 2 3 4 5 5.5 6 6.5 7 8 9 10 12 16 24]

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

BpQ = Operator1S(HilbertObj,GraphObj,@BpQ_OpMatEls);
BmQ = Operator1S(HilbertObj,GraphObj,@BmQ_OpMatEls);

OpOrderA1 = {'L',1,[1; 3]; 'R',1,[2; 3]};
OpOrderA2 = {'L',1,[1; 3]; 'R',2,[2; 3]};

BpBmQ = Operator2S(HilbertObj,GraphObj,@BpBmQ_OpMatEls);

OpOrder1 = {'L',1,[1; 1; 3]; 'R',1,[2; 2; 3]}; % Operator ordering input.
% First row specifies left operator, 1st from configuration, first index is
% spatial separation, 3rd is vector output from matrix elements. Will
% output an Nx1xN array on its own.

% Second row specifies right operator, 1st from configuration, second index
% is spatial separation, 3rd is vector output from matrix elements. Will
% output an 1xNxN array on its own.

OpOrder2 = {'L',1,[1; 1; 3]; 'R',2,[2; 2; 3]}; HOrder = {[],1}; % Operator ordering input.
% HOrder specifies H as 1st from configuration on 'right' (ket-acting) side.

ExOp = OperatorC(HilbertObj,{BpBmQ BpBmQ},OpOrder1); % Output is NxNxN.

HExOp = OperatorH(HilbertObj,HamiltonianObj,{BpBmQ BpBmQ},OpOrder2,HOrder);

ExOpP = OperatorC(HilbertObj,{BpQ BpQ},OpOrderA1);
ExOpM = OperatorC(HilbertObj,{BmQ BmQ},OpOrderA1);

HExOpP = OperatorH(HilbertObj,HamiltonianObj,{BpQ BpQ},OpOrderA2,HOrder);
HExOpM = OperatorH(HilbertObj,HamiltonianObj,{BmQ BmQ},OpOrderA2,HOrder);

SampOperators = {ExOp; HExOp; ExOpP; HExOpP; ExOpM; HExOpM};

SamplerObj = Sampler(HilbertObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,40000); Ncore = 1;

% Initialise storage for the sampled quantities.
EneGS = zeros(numel(U),1); VarE = zeros(numel(U),1);
ExOpEval = cell(numel(U),1); HExOpEval = cell(numel(U),1);

AnsStr = 'BECR-Jast FT'; % Ansatz identifier string.

for u = 1:numel(U)
    if u > 1
        % Alter Hamiltonian parameter.
        HamiltonianObj.HParams(2) = U(u)/2;
        SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    end
    load(['BHM 1D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    % Initiate Monte Carlo sampling to evaluate energy and other
    % observables.
    [EnAvg,EnSamp,EvalAvg,~] = EvalSample(SamplerObj,AnsatzObj);
    % Store sampled quantities.
    EneGS(u) = EnAvg/N; VarE(u) = mean((EnSamp(:)/N - EnAvg/N).^2); % Local per-site energy variance.
    ExOpEval{u} = reshape(EvalAvg{1},N,N,N);
    HExOpEval{u} = reshape(EvalAvg{2},N,N,N);
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
end

save(['BHM 1D N ' num2str(N) ' U Scan ' AnsStr ' excitation observable values.mat'],...
    'EneGS','VarE','ExOpEval','HExOpEval');

exit