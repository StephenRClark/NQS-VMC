% Script for sampling pre-optimised Ansatz objects with the XXZ model and
% any specified Operators.

%% Specify the basic information for the lattice.
L = 4; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

%% Specify Hilbert object.
S = 1/2; Sector = 0; % Set the total Sz to be fixed at zero.
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.
dRXX = [1 0; 0 1]; dRZZ = [1 0; 0 1];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRXX,1); HypCub(Dim,'PBC',dRZZ,1)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
Jx = 1; Jz = [0 0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5...
    0.6 0.7 0.75 0.8 0.85 0.9 0.925 0.95 0.97 0.98 0.99 1];

% Specify the Operators that will be used in the Hamiltonian.
SzSzOp = OperatorDg(HilbertObj,HGraphs(1),@SzSz_CfgVal);
% SzSzOp describes the two-site Sz-Sz operator in the XXZ model. It is
% diagonal in the configuration basis, hence the Dg subclass.
SpSmOp = Operator2S(HilbertObj,HGraphs(2),@SpmSmp_OpMatEls);
% GSzOp describes the transverse field S+S- + h.c. operator in the XXZ
% model. It is not diagonal in the configuration basis and acts on pairs of
% sites, hence the 2S subclass.
Operators = {SzSzOp; SpSmOp}; HParams = [Jz(1), Jx/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.

% Construct the Operator objects for sampling and place in a cell array.

% Relevant operators for the XXZ model are SzSz and S.S operators. The
% first has been made already. Note that the SzSz and S.S operator outputs
% will be given over a range of separations (i.e. 1st entry is nearest
% neighbour, second next-nearest, etc) through the use of EvalSample.

SSOp = Operator2S(HilbertObj,GraphObj,@SiSj_OpMatEls);

SampOperators = {SzSzOp,SSOp};

SamplerObj = Sampler(HilbertObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,10000); 

% Initialise storage for the sampled quantities.
EneGS = zeros(numel(G),1); VarE = zeros(numel(G),1);
SSProfiles = cell(numel(G),1); SzSzProfiles = cell(numel(G),1);

AnsStr = 'Plus-NQSTI Alpha 1'; % Ansatz identifier string.

for z = 1:numel(Jz)
    if z > 1
        % Alter Hamiltonian parameter.
        HamiltonianObj.HParams(2) = Jz(z);
        SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    end
    load(['XXZ 1D Jz ' num2str(U(z)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    % Initiate Monte Carlo sampling to evaluate energy and other
    % observables.
    [EnAvg,EnSamp,EvalAvg,EvalSamp] = EvalSample(SamplerObj,AnsatzObj);
    % Store sampled quantities.
    EneGS(z) = EnAvg/N; VarE(z) = mean((EnSamp(:)/N - EnAvg/N).^2); % Local per-site energy variance.
    SSProfiles{z} = reshape(EvalAvg{3},L,L); SzSzProfiles{z} = reshape(EvalAvg{4},L,L);
    disp(['Sampling for Jz = ' num2str(Jz(z)) ' complete.']);
end

figure(1); hold on; plot(U,real(EneGS));

save(['XXZ 2D N ' num2str(N) ' Jz Scan ' AnsStr ' expectation values.mat'],...
    'EneGS','EnSamp','EvalSamp','SSProfiles','SzSzProfiles');