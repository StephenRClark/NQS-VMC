function BHM_2D_NQS_Op1(L,U)

if ischar(L)
    L = round(str2num(L));
end
if ischar(U)
    U = round(str2num(U));
end

%% Specify the basic information for the lattice.
N = L^2; Dim = [L L];

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
t = -1;

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t/U; 1/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.
SamplerObj1 = Sampler(HilbertObj,HamiltonianObj,{});

SamplerObj1 = SetNsamp(SamplerObj1,400);
SamplerObj1 = SetNblock(SamplerObj1,L); SamplerObj1 = SetNequil(SamplerObj1,1000);

SamplerObj2 = SetNsamp(SamplerObj1,500);
SamplerObj2 = SetNblock(SamplerObj2,2*L); SamplerObj2 = SetNequil(SamplerObj2,2500);

SamplerObj3 = SetNsamp(SamplerObj1,600);
SamplerObj3 = SetNblock(SamplerObj3,N); SamplerObj3 = SetNequil(SamplerObj3,2500);

%% Initialise Stochastic Reconfig Optimiser.
Npass1 = 40; Npass2 = 160; Npass3 = 50; Ncore = 1;
LR1 = 0.01*N; LR2 = 0.005*N; LR3 = 0.001*N;
dEV1 = 0.5/U; dEV2 = 0.2/U; dEV3 = 0.1/U;

SR1 = StochasticReconfig(Npass1,Ncore);
SR2 = StochasticReconfig(Npass2,Ncore);
SR3 = StochasticReconfig(Npass3,Ncore);

SR1 = SetExtraSamples(SR1,0);
SR2 = SetExtraSamples(SR2,0);
SR3 = SetExtraSamples(SR3,0);

SR1 = SetSRTolerances(SR1,2+4/U(1),1e-6);
SR2 = SetSRTolerances(SR2,2+4/U(1),1e-6);
SR3 = SetSRTolerances(SR3,2+4/U(1),1e-6);

SR1 = SetRegularisation(SR1,1e5,1e2,0.8);
SR2 = SetRegularisation(SR2,1e2,1e-2,0.95);
SR3 = SetRegularisation(SR3,1e-2,1e-3,0.95);

SR1 = SetLearnRate(SR1,LR1);
SR2 = SetLearnRate(SR2,LR2);
SR3 = SetLearnRate(SR3,LR3);

SR1 = SetEnergyTolerances(SR1,dEV1,0.1);
SR2 = SetEnergyTolerances(SR2,dEV2,0.1);
SR3 = SetEnergyTolerances(SR3,dEV3,0.1);

%% Initialise starting Ansatz.

% Initialise Reference and Modifier(s) to slot into Ansatz.

% Chosen Reference: BECR
% Necessary fields: Hilbert, Graph, Params.SPH

% Can construct a single particle Hamiltonian (hopping only) using Graph2Array.
CArr = Graph2Array(GraphObj,0); RefParams.SPH = -CArr;

Ref = BECR(HilbertObj,GraphObj,RefParams);

% Chosen Modifier: NQSNHTI
% Necessary fields: Hilbert, Graph, Params.(a, b, c, nmag, nphs, A, B,
% Nh/Alpha), VFlag (set to 1)

Alpha = 1;
ModParams.Nh = Alpha * N; ModParams.a = 0; ModParams.b = 0; ModParams.W = -0.01;
ModParams.nmag = 0.1; ModParams.nphs = 0; ModParams.A = -0.1; ModParams.B = 0;

NQSObj = NQSSHTI(HilbertObj,GraphObj,ModParams,1);

ExInds = [1 3 4]; NQSObj.OptInds(ExInds) = 0;

dP = [zeros(4,1); -10*ModParams.W; zeros(Alpha*N-1,1)]; NQSObj = PsiUpdate(NQSObj,dP);

% Initialise Ansatz object. Mods need to be passed as a cell list.
Mod = {NQSObj};

AnsatzObj = Ansatz(Ref,Mod,HilbertObj);

AnsStr = 'BECR-NQSSHTI-HD5 Alpha 1 Op3'; % Ansatz name string for easy identification later.

tic;
% Perform optimisation of AnsatzObj with SR.
EnIter = cell(3,1);
[AnsatzObj,EnIter{1}] = SR1.Optimise(SamplerObj1,AnsatzObj);
[AnsatzObj,EnIter{2}] = SR2.Optimise(SamplerObj2,AnsatzObj);
[AnsatzObj,EnIter{3}] = SR3.Optimise(SamplerObj3,AnsatzObj);
RunTime = toc; EnIter = cell2mat(EnIter); Params = AnsatzObj.ParamList;
disp(['Total run time: ' num2str(RunTime) ' seconds.']);
if isdeployed
    savepath = strcat(pwd,filesep,'Output');
    if isdir(savepath) == 0
        mkdir(savepath);
    end
    savepath = strcat(savepath,filesep);
else
    savepath = [];
end
save([savepath 'BHM 2D U ' num2str(U) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
    'AnsatzObj','EnIter','RunTime','Params','SamplerObj1','SamplerObj2','SamplerObj3','SR1','SR2','SR3');
% Save AnsatzObj and run details for later analysis.
end