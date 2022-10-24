%% Set some options for the demonstration

show_figures = true;
short_on_time = false;

%% Add relevant folders to path

addpath(genpath("../../.")); % Adds all folders of MATLAB section.

%% Generate desired geometry with Graph

% Chosen lattice geometry: 4x4 grid, nearest neighbour connectivity

Dim = [4 4]; Bound = [1 1]; % (1 x Ndim) vectors for dimensions and boundary conditions.
LVecs = [1 0; 0 1]; % Primitive lattice vectors, nearest neighbours in x and y.
SFlag = 1; % Tells Graph object to generate full map of site connections in Graph.BondMap.
N = prod(Dim); % Total number of sites.

GraphObj = HypCub(Dim,Bound,LVecs,SFlag); % Create HypCub (hypercube) instance.

% Can show connections with Graph2Array. Option 0 as second argument
% specifies an ordinary lattice suited for spins / bosons.
GraphMat = Graph2Array(GraphObj,0);
if show_figures
    figure; imagesc(GraphMat); xlabel('Site index'); ylabel('Site index');
    title('Connectivity matrix of chosen Graph');
end

%% Parameterise configuration generation with Hilbert

% Chosen system: spin-1/2s, fixed total Sz projection of 0.
S = 1/2; Sector = 0; % Specify spin magnitude S and desired spin Sector.

HilbertObj = Spin(N,S,Sector); % Create Spin Hilbert instance.

% Generate and example configuration:
Cfg_example = HilbertObj.RandomCfg();
% FullCfg generates a vector of spin projections from site 1 to site N.
Cfg_vec = HilbertObj.FullCfg(Cfg_example);
% Verify total spin projection is zero:
disp(['Sum of total spin projections: ' num2str(sum(Cfg_vec))]);
if show_figures
    figure; imagesc(reshape(Cfg_vec,Dim(1),Dim(2))); xlabel('x'); ylabel('y');
    title('Example spin configuration');
end

%% Create Hamiltonian using Operators

% Chosen Hamiltonian: 2D nearest neighbour spin-1/2 Heisenberg interactions
% Heisenberg Hamiltonian has a single parameter J:
J = 1;
% Hamiltonian parameters passed in HParams.
HParams = J;

% Heisenberg interaction is full spin-spin dot product across site pairs.
% The Heisenberg model has a sign structure that can be accounted for with
% a gauge transform, built into one of the Operator functions.

% Create an Operator2S (2 site) with SiSj_GT_OpMatEls as matrix element
% generating function. Other matrix element functions can be found in
% 'Operators/Operator matrix elements', or for Operators diagonal in the
% configuration basis, look in 'Operators/Operator configuration values'.
MatElFunc = @SiSj_GT_OpMatEls;

% Need to add custom sublattice indexing to GraphObj for gauge transform -
% can do so by modifying Graph.SLInds.
for i = 1:Dim(1)
    for j = 1:Dim(1)
        GraphObj.SLInds(i+(j-1)*Dim(1)) = 1 + mod(i+j,2);
    end
end

% Need to pass Hilbert, Graph with desired connectivity and matrix element 
% function handle for an Operator instance.
SiSjOp = Operator2S(HilbertObj,GraphObj,MatElFunc);

% Operators create a set of Diff structs and matrix element values when
% passed a configuration.
[DiffOp,OpMatEls] = SiSjOp.CorrMatEls(Cfg_example);
% Diffs are used to determine ratios of wavefunction amplitudes.

% Create Hamiltonian instance by passing cell array of operators, and
% HParams.
OperatorArray = {SiSjOp};
HamObj = Hamiltonian(OperatorArray,HParams);

%% Set up variational wavefunction in Ansatz

% Chosen wavefunction: NQS with one hidden layer of hidden unit density 4,
% on an equal superposition of basis states.

% Equal superposition of basis states is the Plus reference.
Ref = Plus();

% NQS belongs to the Modifier class. Modifiers require Hilbert, Graph, and
% a struct of Modifier parameters determining how the starting NQS
% parameters are initialised. The last argument VFlag is a flag to allow /
% prevent optimisation of variational parameters. More fine-tuned
% optimisation flags on a per-parameter basis can be found in
% Modifier.OptInds.
VFlag = 1; 

% Noise terms determine degree of randomness. Random values from 0 to nmag
% are added, then a random phase between -nphs to nphs is applied as an
% exp(1i*phase) multiplier. Want small randomness added, and only real
% parameters.
ModParams.nmag = 0.01; ModParams.nphs = 0;

% NQS requires additional parameters set: Alpha (hidden unit density), a,
% b, W (starting average values of these parameters).
ModParams.a = 0.01; ModParams.b = 0.01; ModParams.W = 0.01; ModParams.Alpha = 4;

% Generate NQS modifier instance. NQS modifier will retain the symmetries
% present in the Graph, in this case it is invariant under translations in
% x and y. For no translational invariance, pass a Graph with one lattice
% vector of zeros, describing disconnected sites.
NQSObj = NQS(HilbertObj,GraphObj,ModParams,VFlag);

% Full Ansatz is generated by passing the Reference, a cell array of
% Modifiers, and the Hilbert for the system.
ModArray = {NQSObj};
AnsatzObj = Ansatz(Ref,ModArray,HilbertObj);

%% Set Monte Carlo parameters with Sampler

% Sampler needs Hilbert, Hamiltonian, and can be optionally passed a cell
% array of Operators to sample at appropriate points in the Markov chain.
% No Operators require sampling during optimisation, passing empty cells.
OptSampler = Sampler(HilbertObj,HamObj,{});

% Chosen parameters for OptSampler: 5000 samples per sampling phase, with 
% samples spaced by N to decorrelate and 2500 burn-in / equilibration
% Markov steps.
Nsamp_opt = 5000; Nblock_opt = N; Nequil_opt = 2500;

OptSampler = OptSampler.SetNsamp(Nsamp_opt);
OptSampler = OptSampler.SetNblock(Nblock_opt);
OptSampler = OptSampler.SetNequil(Nequil_opt);

% Set up a separate EvalSampler for evaluation of observables after the
% optimisation. Need SampleOperators set up.

SzSzOp = OperatorDg(HilbertObj,GraphObj,@SzSz_CfgVal); % Two site Sz correlation, diagonal in Sz basis.
SzOp = OperatorDg(HilbertObj,GraphObj,@Sz_CfgVal); % Single site Sz average, diagonal in Sz basis.
SxOp = Operator1S(HilbertObj,GraphObj,@Sx_OpMatEls); % Single site Sx average, non-diagonal in Sz basis.
SampleOperators = {SzSzOp;SzOp;SxOp};

EvalSampler = Sampler(HilbertObj,HamObj,SampleOperators);

% Chosen parameters for EvalSampler: 40000 samples per sampling phase, with 
% samples spaced by N to decorrelate and 5000 burn-in / equilibration
% Markov steps.
Nsamp_samp = 40000; Nblock_samp = N; Nequil_samp = 5000;

EvalSampler = EvalSampler.SetNsamp(Nsamp_samp);
EvalSampler = EvalSampler.SetNblock(Nblock_samp);
EvalSampler = EvalSampler.SetNequil(Nequil_samp);

%% Set optimisation parameters with StochasticReconfig

% Chosen optimisation method: stochastic reconfiguration, for 500
% optimisation steps on four cores. Multi-stage optimisations can be set up
% using multiple StochasticReconfig objects.
Npass = 500; Ncore = 4;

SRObj = StochasticReconfig(Npass,Ncore);

% Stochastic reconfiguration has several hyperparameters that influence the
% stability of the optimisation. The chosen parameters here are a starting
% regularisation of 10^4, minimum regularisation of 10^-3, regularisation
% decay factor of 0.96 and a learning rate of 1. 
RegMax = 10^4; RegMin = 10^-3; RegMult = 0.96; LearnRate = 1;

SRObj = SRObj.SetRegularisation(RegMax,RegMin,RegMult);
SRObj = SRObj.SetLearnRate(LearnRate);

%% Pass Ansatz object for optimisation

if ~short_on_time
    % Passing an Ansatz object to an Optimiser with a Sampler will give an
    % updated Ansatz object and a vector of sampled energies by iteration.
    [AnsatzOpt, EnIter] = SRObj.Optimise(OptSampler,AnsatzObj); 
    save('Heisenberg 2D example optimisation.mat','AnsatzObj','AnsatzOpt','EnIter');
else
    load('Heisenberg 2D example optimisation.mat');
end
if show_figures
    figure; plot(EnIter); xlabel('Iteration'); ylabel('Sampled energy');
end

if ~short_on_time
    % After an optimsation, pass the Ansatz object to a Sampler with
    % assigned Operators to evaluate observables. Using MultiChainSample
    % will utilise multiple cores through parpool.
    % Sampling methods will return mean energy of the samples, energy per
    % sample, mean values of sampled Operators, and Operator values per
    % sample. The per-sample variables are omitted here to avoid clutter.
    [EnAvg,~,EvalAvg,~] = EvalSampler.MultiChainSample(AnsatzOpt,Ncore);
    % EnAvg is total system energy, but normally energy density is more
    % interesting.
    EnDen = EnAvg/N; 
    % Unpack the other observables into more presentable forms.
    SzSzAvg = reshape(EvalAvg{1},N,N);
    SzAvg = EvalAvg{2}; SxAvg = EvalAvg{3};
    save('Heisenberg 2D example evaluation.mat','AnsatzObj','EnDen','SzSzAvg','SzAvg','SxAvg');
else
    load('Heisenberg 2D example evaluation.mat')
end
disp(['Average energy per site: ' num2str(EnDen)]);
disp(['Average Sz projection: ' num2str(SzAvg)]);
disp(['Average Sx projection: ' num2str(SxAvg)]);
if show_figures
    figure; imagesc(SzSzAvg); xlabel('Site index'); ylabel('Site index');
end