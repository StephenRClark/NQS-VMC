% Script for sampling pre-optimised Ansatz objects with the Ising model and
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
S = 1/2; Sector = []; % Allow the total Sz of the system to change freely.
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.
dRJ = [1 0; 0 1]; dRG = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRJ,1); HypCub(Dim,'PBC',dRG,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
J = -1; G = [0 1/2 1 3/2 2 3 4 5 6 8 10 12 16 20 25];

% Specify the Operators that will be used in the Hamiltonian.
SzSzOp = OperatorDg(HilbertObj,HGraphs(1),@SzSz_CfgVal);
% SzSzOp describes the two-site Sz-Sz operator in the Ising model. It is
% diagonal in the configuration basis, hence the Dg subclass.
GSzOp = Operator1S(HilbertObj,HGraphs(2),@Sx_OpMatEls);
% GSzOp describes the transverse field Sx operator in the Ising model. It
% is not diagonal in the configuration basis and acts on single sites,
% hence the 1S subclass.
Operators = {SzSzOp; GSzOp}; HParams = [J, J*G(1)];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.

% Construct the Operator objects for sampling and place in a cell array.

% Relevant operators for the Ising model are total Sz, Sz profile, Sx and
% SzSz operators. The latter two have been made already. Note that the SzSz
% operator output will be given over a range of separations (i.e. 1st entry
% is nearest neighbour, second next-nearest, etc) through the use of
% EvalSample.

SzOp = OperatorDg(HilbertObj,GraphObj,@Sz_CfgVal);

SzProfOp = OperatorDg(HilbertObj,GraphObj,@DenProf_CfgVal);

SampOperators = {SzOp, SxOp, SzProfOp, SzSzOp};

SamplerObj = Sampler(HilbertObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,10000); 

% Initialise storage for the sampled quantities.
EneGS = zeros(numel(G),1); VarE = zeros(numel(G),1);
SzTotal = zeros(numel(G),1); SxTotal = zeros(numel(G),1);
SzProfiles = cell(numel(G),1); SzSzProfiles = cell(numel(G),1);

AnsStr = 'Plus-NQSTI Alpha 1'; % Ansatz identifier string.

for g = 1:numel(G)
    if g > 1
        % Alter Hamiltonian parameter.
        HamiltonianObj.HParams(2) = J*G(g);
        SamplerObj = SetHamiltonian(SamplerObj,HamiltonianObj);
    end
    load(['Ising 2D G ' num2str(U(g)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    % Initiate Monte Carlo sampling to evaluate energy and other
    % observables.
    [EnAvg,EnSamp,EvalAvg,EvalSamp] = EvalSample(SamplerObj,AnsatzObj);
    % Store sampled quantities.
    EneGS(g) = EnAvg/N; VarE(g) = mean((EnSamp(:)/N - EnAvg/N).^2); % Local per-site energy variance.
    SzTotal(g) = EvalAvg{1}; SxTotal(g) = EvalAvg{2};
    SzProfiles{g} = reshape(EvalAvg{3},L,L); SzSzProfiles{g} = reshape(EvalAvg{4},L,L);
    disp(['Sampling for G = ' num2str(G(g)) ' complete.']);
end

figure(1); hold on; plot(U,real(EneGS));

save(['Ising 2D N ' num2str(N) ' G Scan ' AnsStr ' expectation values.mat'],...
    'EneGS','EnSamp','EvalSamp','SzTotal','SxTotal','SzProfiles','SzSzProfiles');