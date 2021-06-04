% --- Markov Chain Monte Carlo sampling function ---

function [EnAvg,dLogpAvg,EvalAvg,MRate] = MCMCSample(SamplerObj,AnsatzObj)

% Object-orientation version, which feeds Hamiltonian and SampOperators
% into a Sampler object. Ansatz contains all the information relevant to
% the variational wavefunction, and both contain references to the Hilbert
% space of configurations that they occupy.

Nrun = SamplerObj.Nsamp * SamplerObj.Nblock; % Total number of MCMC steps.

Cfg = AnsatzObj.Hilbert.RandomCfg; % Generate an initial configuration Cfg for Markov chain.
AnsatzObj = AnsatzObj.PrepPsi(Cfg); % Prepare the ansatz by computing intermediate information for the initial configuration.
PropMove = AnsatzObj.Hilbert.PropMoveFunc;
% Initialise averages of sampled quantities:
EnAvg = 0; MRate = 0; % Average energy and move acceptance rate accumulators.
dLogpAvg = zeros(AnsatzObj.NpTotal,1); % Average log-derivative vector.
NumEvals = numel(SamplerObj.Operators); % Number of additional quantities for averaging.
EvalAvg = cell(NumEvals,1); % Create a cell array for storing the additional averaged quantities.
% Run the evaluation functions once, multiplying the output by zero, to
% initialise each element of the cell array.
for k=1:NumEvals
    EvalAvg{k} = 0*SamplerObj.Operators{k}.LocalSample(Cfg,0,0,AnsatzObj);
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
end
for q = 1:SamplerObj.Nsamp
    % First Nequil runs are to reach "burn in" equilibrium for the chain.
    for n = 1:SamplerObj.Nblock
        % Perform Markov chain Monte Carlo step:
        [Diff,CfgP] = PropMove(Cfg); % Propose a move to a new configuration CfgP.
        % Compute ansatz wave function amplitude ratio Psi(CfgP)/Psi(Cfg):
        [PsiRatio,CfgP_update] = AnsatzObj.PsiRatio(Diff); % Intermediate ansatz information CfgP_update is also returned for updates.
        PChange = min(abs(PsiRatio)^2,1); % Compute the Metropolis acceptance probability min(|Psi(CfgP)|^2/|Psi(Cfg)|^2,1).
        if rand<=PChange
            % Move is accepted:
            Cfg = CfgP; % Update current configuration.
            MRate = MRate + (1/Nrun);
            AnsatzObj = AnsatzObj.PsiCfgUpdate(CfgP_update); % Apply the update of the configuration state information in the ansatz.
        end
    end
    % --- Perform evaluation of standard quantities:
    % (1) Start with the average energy via the local energy estimator:
    [EnLoc] = SamplerObj.Hamiltonian.EnergySample(Cfg,AnsatzObj);
    EnAvg = EnAvg + (1/SamplerObj.Nsamp)*EnLoc; % The average of the local energy (once all samples are made).    
    % (2) Next compute the average of the logarithmic derivative of the ansatz:
    dLogp = AnsatzObj.LogDeriv(Cfg);
    dLogpAvg = dLogpAvg + (1/SamplerObj.Nsamp)*dLogp; % The average of the logarithmic derivative vector (once all samples are made).    
    % --- Now compute other estimators as specified by the evaluation function:
    for k=1:NumEvals
        EvalAvg{k} = EvalAvg{k} + (1/SamplerObj.Nsamp)*SamplerObj.Operators{k}.LocalSample(Cfg,EnLoc,dLogp,AnsatzObj);
    end
end
end
