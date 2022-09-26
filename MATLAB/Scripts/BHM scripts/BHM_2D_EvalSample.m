% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

%% Specify the basic information for the lattice.
L = 6; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 2D L x L lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours in x and
% y, as well as for multiple nearest neighbour separations.

NStr = {' N-1 ',' N ',' N+1 '};
Nvec = [N-1 N N+1];

for nb = 2 % 1:numel(Nvec)
%% Specify Hilbert object.
Nb = Nvec(nb); Nmax = 4;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = [1 0; 0 1]; dRU = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

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

SamplerObj = SetNsamp(SamplerObj,40000); Ncore = 8;

% Initialise storage for the sampled quantities.
EneGSU = zeros(numel(U),1); VarEU = zeros(numel(U),1);
VarNU = zeros(numel(U),1); NiNjU = cell(numel(U),1);
DbHlU = cell(numel(U),1); OcFrU = zeros(numel(U),Nmax+1);
BiBjU = cell(numel(U),1);

% AnsStr = 'BECR-Jast-NNMB SeqOpt';
% AnsStr = 'BECR-NQSSHTI-HD5-AWGS Alpha 1'; % Ansatz identifier string.
AnsStr = 'BECR-NQSMHTI-HD5-LR Alpha 1';

for u = 17 % 1:numel(U)
    if u > 1
        % Alter Hamiltonian parameter corresponding to interaction.
        HamiltonianObj.HParams(2) = U(u)/2;
        SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    end
    load(['BHM 2D U ' num2str(U(u)) NStr{nb} num2str(N) ' ' AnsStr ' Logs.mat']);
    % Initiate Monte Carlo sampling to evaluate energy and other
    % observables.
    [EnAvg,EnSamp,EvalAvg,~] = MultiChainSample(SamplerObj,AnsatzObj,Ncore);
    % Store sampled quantities.
    EneGSU(u) = EnAvg/N; VarNU(u) = EvalAvg{1}; NiNjU{u} = reshape(EvalAvg{2},L,L);
    DbHlU{u} = reshape(EvalAvg{3},L,L); OcFrU(u,:) = reshape(EvalAvg{4},1,Nmax+1);
    VarEU(u) = mean((EnSamp(:)/N - EnAvg/N).^2); BiBjU{u} = EvalAvg{5}; % Local per-site energy variance.
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
    EneGS = EneGSU(u); VarN = VarNU(u); VarE = VarEU(u); NiNj = NiNjU{u};
    DbHl = DbHlU{u}; OcFr = OcFrU(u,:).'; BiBj = BiBjU{u};
    save(['BHM 2D U ' num2str(U(u)) NStr{nb} num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','EneGS','VarN','NiNj','VarE','DbHl','OcFr','BiBj');
end

% for u = 1:numel(U)
%     load(['BHM 2D U ' num2str(U(u)) NStr{nb} num2str(N) ' ' AnsStr ' Logs.mat']);
%     EneGSU(u) = EneGS; VarNU(u) = VarN; VarEU(u) = VarE; NiNjU{u} = NiNj;
%     DbHlU{u} = DbHl; OcFrU(u,:) = OcFr.'; BiBjU{u} = BiBj;
% end

% figure(1); hold on; plot(U,real(EneGSU));

% save(['BHM 2D U Scan' NStr{nb} num2str(N) ' ' AnsStr ' MMC expectation values.mat'],...
%     'EneGSU','VarNU','NiNjU','VarEU','DbHlU','OcFrU','BiBjU');

end