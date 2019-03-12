% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

%% Specify the basic information for the lattice.
N = 20; Dim = N;

%% Specify Graph object.
LVecs = 1;
GraphObj = HypCub(Dim,'PBC',LVecs);

%% Specify Hilbert object.
Nb = N; Nmax = 10;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT); HypCub(Dim,'PBC',dRU)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [2 2.1 2.2 2.25 2.3 2.35 2.375 2.4 2.425 2.45 2.475 2.5 2.525 2.55 2.575 2.6 2.7 2.8 2.9 3];
HamiltonianObj = Hubbard(HilbertObj,HGraphs,t,U(1)); 

%% Initialise Sampler.

% Construct the Operator objects for sampling and place in a cell array. In
% this case, operators are variance in site occupation, two-site density-density
% correlations and static structure factor.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@VarN_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(1),@NiNj_Bose_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@StatSF_Bose_CfgVal)};

SamplerObj = Sampler(HilbertObj,GraphObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,10000);

% Initialise storage for the sampled quantities.
EneGS = zeros(numel(U),1); VarE = zeros(numel(U),1);
VarN = zeros(numel(U),1); NiNj = zeros(numel(U),N);
StatSF = zeros(numel(U),N);

AnsStr = 'BECR-NQSTI-HD5 Alpha 1'; % Ansatz identifier string.

for u = 1:numel(U)
    if u > 1
        % Create new Hamiltonian and cycle out old one.
        NewH = Hubbard(HilbertObj,HGraphs,t,U(u));
        SamplerObj = SetHamiltonian(SamplerObj,NewH);
    end
    tic;
    load(['BHM U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj','EnIter','RunTime');
    [EnAvg,EnSamp,EvalAvg,~] = EvalSample(SamplerObj,AnsatzObj);
    % Store sampled quantities.
    EneGS(u) = EnAvg/N; VarN(u) = EvalAvg{1}; NiNj(u,:) = EvalAvg{2};
    StatSF(u,:) = EvalAvg{3}.';
    VarE(u) = mean(((EnSamp - EnAvg)/EnAvg).^2); % Normalised local energy variance.
    disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
end

figure(1); plot(U,real(EneGS));

save(['BHM N ' num2str(N) ' U Scan ' AnsStr ' expectation values.mat'],'EneGS','VarN','NiNj','VarE','StatSF');