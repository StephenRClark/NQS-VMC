% Script for constructing and optimising Ansatz objects with number-hidden
% NQS Modifiers and Bose condensate References for the Bose Hubbard model
% on a 1D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
L = 4; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 2D L x L lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours in x and
% y, as well as for multiple nearest neighbour separations.

%% Specify Hilbert object.
Nb = N; NStr = ' N '; Nmax = 4;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = [1 0; 0 1]; dRU = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = 1;

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t; U/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

EvalHamiltonian = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj1 = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj1 = SetNsamp(SamplerObj1,2500);
SamplerObj1 = SetNblock(SamplerObj1,N/2); SamplerObj1 = SetNequil(SamplerObj1,1000);

SamplerObj2 = SetNsamp(SamplerObj1,5000);
SamplerObj2 = SetNblock(SamplerObj2,N); SamplerObj2 = SetNequil(SamplerObj2,2500);

SamplerObj3 = SetNsamp(SamplerObj1,10000);
SamplerObj3 = SetNblock(SamplerObj3,N); SamplerObj3 = SetNequil(SamplerObj3,5000);
% Emulating usage in Optimiser passes.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@SR_EnLogDerivCorr),...
    OperatorDg(HilbertObj,HGraphs(2),@SR_LogDerivCorr)};

SamplerObj1 = SetOperator(SamplerObj1,SampOperators{1},1); SamplerObj1 = SetOperator(SamplerObj1,SampOperators{2},2);
SamplerObj2 = SetOperator(SamplerObj2,SampOperators{1},1); SamplerObj2 = SetOperator(SamplerObj2,SampOperators{2},2);
SamplerObj3 = SetOperator(SamplerObj3,SampOperators{1},1); SamplerObj3 = SetOperator(SamplerObj3,SampOperators{2},2);

%% Initialise starting Ansatz.

% Chosen Modifier: NQSSHTI
% Necessary fields: Hilbert, Graph, Params.(a, b, c, nmag, nphs, A, B,
% Nh/Alpha), VFlag (set to 1)
CArr = Graph2Array(GraphObj,0); RefParams.SPH = -CArr;
Ref = BECR(HilbertObj,GraphObj,RefParams);

Alpha = 1;
ModParams.Nh = Alpha * N; ModParams.a = 0; ModParams.b = 0.1; ModParams.W = -0.01;
ModParams.nmag = 0; ModParams.nphs = 0; ModParams.A = -0.1; ModParams.B = -0.1;
NQSObj = NQSSHTI(HilbertObj,GraphObj,ModParams,1);

AnsatzObj = Ansatz(Ref,{NQSObj},HilbertObj);
for f = 1:5
tic;
[EnAvg1,dLogpAvg1,EvalAvg1,MRate1] = SamplerObj1.MCMCSample(AnsatzObj);
SampTime1 = toc; disp(['First sampling stage time: ' num2str(SampTime1)]);
tic;
[EnAvg2,dLogpAvg2,EvalAvg2,MRate2] = SamplerObj2.MCMCSample(AnsatzObj);
SampTime2 = toc; disp(['Second sampling stage time: ' num2str(SampTime2)]);
tic;
[EnAvg3,dLogpAvg3,EvalAvg3,MRate3] = SamplerObj3.MCMCSample(AnsatzObj);
SampTime3 = toc; disp(['Third sampling stage time: ' num2str(SampTime3)]);

save(['BHM 2D N ' num2str(N) ' sampling performance tests ' num2str(f) '.mat'],...
    'EnAvg1','SampTime1','SamplerObj1','MRate1',...
    'EnAvg2','SampTime2','SamplerObj2','MRate2',...
    'EnAvg3','SampTime3','SamplerObj3','MRate3');
end