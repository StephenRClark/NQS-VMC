%% Specify the basic information for the lattice.
L = 10; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1; Bound = [1 1];
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,Bound,LVecs,SFlag);
% GraphObj describes a 2D L x L lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours in x and
% y, as well as for multiple nearest neighbour separations.

%% Specify Hilbert object.
Nb = N; NStr = ' N '; Nmax = 4;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = [1 0; 0 1]; dRU = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,Bound,dRT,1); HypCub(Dim,Bound,dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [1 14 15 16 17 18 23 32];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t/(U(1)); 1/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

EvalHamiltonian = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj1 = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj1 = SetNsamp(SamplerObj1,8000);
SamplerObj1 = SetNblock(SamplerObj1,N); SamplerObj1 = SetNequil(SamplerObj1,5000);

SamplerObj2 = SetNsamp(SamplerObj1,9600);
SamplerObj2 = SetNblock(SamplerObj2,N); SamplerObj2 = SetNequil(SamplerObj2,5000);

% Set up evaluation Sampler.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@VarN_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@NiNj_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@DbHl_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@OccFrac_Bose_CfgVal),...
    Operator2S(HilbertObj,HGraphs(1),@BiBj_OpMatEls),...
    OperatorDg(HilbertObj,HGraphs(1),@DiDj_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@HiHj_Bose_CfgVal)};

EvalSampler = Sampler(HilbertObj,HamiltonianObj,SampOperators);

EvalSampler = SetNsamp(EvalSampler,40000);

% Binning analysis parameters for energy errorbar.
Nbin = 40; Lbin = 1000; EneBin = zeros(Nbin,1);

%% Initialise Stochastic Reconfig Optimiser.
Ncore = 16;

Npass1 = 600; Npass2 = 400; Nfloor = 50; dEV1 = 0.05/U(1); dEV2 = 0.01/U(1);

SR1 = StochasticReconfig(Npass1,Ncore); SR1 = SetExtraSamples(SR1,4000); 
SR2 = StochasticReconfig(Npass2,Ncore); SR2 = SetExtraSamples(SR2,4000); 

RegMax1 = [1e3 1e3 1e3 1e3 1e3 1e3 1e3 1e3]; RegMin1 = [1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1 1e-1]; 
LearnRate1 = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; RegFac1 = exp((log(RegMin1)-log(RegMax1))/(Npass1-Nfloor)); 

RegMax2 = [1e1 1e1 1e1 1e1 1e1 1e1 1e1 1e1]; RegMin2 = [1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3 1e-3]; 
LearnRate2 = [0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01]; RegFac2 = exp((log(RegMin2)-log(RegMax2))/(Npass2-Nfloor));  