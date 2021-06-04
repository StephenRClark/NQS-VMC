% --- Markov Chain Monte Carlo sampling function ---

function [EnAvg,EnSamp,EvalAvg,EvalSamp] = EvalSample(SamplerObj,AnsatzObj)

% Object-orientation version, which feeds Hamiltonian and SampOperators
% into a Sampler object. Ansatz contains all the information relevant to
% the variational wavefunction, and both contain references to the Hilbert
% space of configurations that they occupy.

% EvalSample differs from MCMCSample when it comes to sampling operators -
% EvalSample will produce an Nsamp cell list of values for error analysis
% purposes. Additionally, the EvalSamp functions will invoke GraphSample
% for each operator.

% EvalSample is not meant to be used during optmisation, and does not
% output dLogpAvg. MCMCSample is to be used during optimisation.

Nrun = SamplerObj.Nsamp * SamplerObj.Nblock; % Total number of MCMC steps.

Cfg = AnsatzObj.Hilbert.RandomCfg; % Generate an initial configuration Cfg for Markov chain.
AnsatzObj = AnsatzObj.PrepPsi(Cfg); % Prepare the ansatz by computing intermediate information for the initial configuration.
HilbertObj = AnsatzObj.Hilbert; PropMove = HilbertObj.PropMoveFunc;
% Initialise averages of sampled quantities:
EnAvg = 0; MRate = 0; % Average energy and move acceptance rate accumulators.
NumEvals = numel(SamplerObj.Operators); % Number of additional quantities for averaging.
EvalAvg = cell(NumEvals,1); % Create a cell array for storing the additional averaged quantities.
EnSamp = zeros(1,SamplerObj.Nsamp); % Lists local energy sample values.
EvalSamp = cell(NumEvals,SamplerObj.Nsamp); % Cell array for storing samples of the evaluated quantities.
% Run the evaluation functions once, multiplying the output by zero, to
% initialise each element of the cell array.

for k=1:NumEvals
    EvalAvg{k} = 0*SamplerObj.Operators{k}.GraphSample(Cfg,0,0,AnsatzObj);
end

% Run the Markov chain steps:
for q=1:SamplerObj.Nequil
    % Perform Markov chain Monte Carlo step:
    [Diff,CfgP] = PropMove(Cfg); % Propose a move to a new configuration CfgP.
    % Compute ansatz wave function amplitude ratio Psi(CfgP)/Psi(Cfg):
    [PsiRatio,CfgP_update] = AnsatzObj.PsiRatio(Diff); % Intermediate ansatz information CfgP_update is also returned for updates.
    PChange = min(abs(PsiRatio)^2,1); % Compute the Metropolis acceptance probability min(|Psi(CfgP)|^2/|Psi(Cfg)|^2,1).
    if rand<=PChange
        % Move is accepted:
        Cfg = CfgP; % Update current configuration.
        AnsatzObj = AnsatzObj.PsiCfgUpdate(CfgP_update); % Apply the update of the configuration state information in the ansatz.
    end
    % First Nequil runs are to reach "burn in" equilibrium for the chain.
end
for q=1:SamplerObj.Nsamp
    for n = 1:SamplerObj.Nblock
        % Perform Markov chain Monte Carlo step:
        [Diff,CfgP] = PropMove(Cfg); % Propose a move to a new configuration CfgP.
        % Compute ansatz wave function amplitude ratio Psi(CfgP)/Psi(Cfg):
        [PsiRatio,CfgP_update] = AnsatzObj.PsiRatio(Diff); % Intermediate ansatz information CfgP_update is also returned for updates.
        PChange = min(abs(PsiRatio)^2,1); % Compute the Metropolis acceptance probability min(|Psi(CfgP)|^2/|Psi(Cfg)|^2,1).
        if rand<=PChange
            % Move is accepted:
            Cfg = CfgP; % Update current configuration.
            MRate = MRate + (1/Nrun); % Accumulator for end move acceptance rate.
            AnsatzObj = AnsatzObj.PsiCfgUpdate(CfgP_update); % Apply the update of the configuration state information in the ansatz.
        end
    end
    % --- Perform evaluation of standard quantities:
    % (1) Start with the average energy via the local energy estimator:
    [EnLoc] = SamplerObj.Hamiltonian.EnergySample(Cfg,AnsatzObj);
    EnAvg = EnAvg + (1/SamplerObj.Nsamp)*EnLoc; % The average of the local energy (once all samples are made).
    EnSamp(q) = EnLoc;
    % --- Now compute other estimators as specified by the evaluation function:
    for k=1:NumEvals
        SampVal = SamplerObj.Operators{k}.GraphSample(Cfg,EnLoc,0,AnsatzObj);
        EvalAvg{k} = EvalAvg{k} + (1/SamplerObj.Nsamp)*SampVal;
        EvalSamp{k,n} = SampVal;
    end
end
end
