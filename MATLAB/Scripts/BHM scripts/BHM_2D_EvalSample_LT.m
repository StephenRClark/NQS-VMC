% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

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
Nmax = 4;
NStr = ' N '; Nb = N;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = [1 0; 0 1]; dRU = [0 0]; Bound = [1 1];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,Bound,dRT,1); HypCub(Dim,Bound,dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [0];

% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t; U(1)/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

SysStr = 'BHM 2D U ';

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

Nbin = 40; Lbin = 1000; EneBin = zeros(Nbin,1);

SamplerObj = SetNsamp(SamplerObj,40000); Ncore = 4;

% Initialise storage for the sampled quantities.
EneGSU = zeros(numel(U),1); VarEU = zeros(numel(U),1);
VarNU = zeros(numel(U),1); NiNjU = cell(numel(U),1);
DbHlU = cell(numel(U),1); OcFrU = zeros(numel(U),Nmax+1);
BiBjU = cell(numel(U),1);

% AnsStr = 'BECR-Gutz';
% AnsStr = 'BECR-NQSSHTI-HD5 Alpha 2 SR';
% AnsStr = 'BECR-NQS-HD5-ABSH-DHSR Alpha 2';
% AnsStr = 'BECR-NQSNHTI-HD5-NNMB Alpha 2';
% AnsStr = 'BECR-NQSMHTI-HD5-LR Alpha 1';
% AnsStr = 'BECR-Jast-NNMB FT';
% AnsStr = 'BECR-NQSTI-I Alpha 1';
% AnsStr = 'BECR-NQSOHTI-VD5 Alpha 1';
AnsStr = 'BECR';

RefParams.SPH = -Graph2Array(GraphObj,0);
Ref = BECR(HilbertObj,GraphObj,RefParams);
AnsatzObj = Ansatz(Ref,{None()},HilbertObj);

urange = 1;

for u = urange
    % Alter Hamiltonian parameter corresponding to interaction.
    HamiltonianObj.HParams(2) = U(u)/2;
    SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    % load([SysStr num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat']);
    % Initiate Monte Carlo sampling to evaluate energy and other
    % observables.
    tic;
    [EnAvg,EnSamp,EvalAvg,~] = EvalSample(SamplerObj,AnsatzObj);
    EvalTime = toc;   
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
    % Save observable values in file containing Ansatz and run logs.
    EneGS = EnAvg/N; VarN = EvalAvg{1}; EnSamp = EnSamp(:)/N;
    for n = 1:Nbin
        EneBin(n) = mean(EnSamp((1:Lbin)+(n-1)*Lbin));
    end
    VarE = mean((EneBin-mean(EneBin)).^2);
    NiNj = reshape(EvalAvg{2},L,L); DbHl = reshape(EvalAvg{3},L,L); 
    OcFr = reshape(EvalAvg{4},Nmax+1,1); BiBj = EvalAvg{5}; RunDate = date;
    save([SysStr num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','EvalTime','RunDate','VarE',...
        'VarN','NiNj','DbHl','OcFr','BiBj');
end