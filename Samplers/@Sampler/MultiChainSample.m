% --- Markov Chain Monte Carlo sampling function ---

function [EnAvg,EnSamp,EvalAvg,EvalSamp] = MultiChainSample(SamplerObj,AnsatzObj,Ncore)

% Object-orientation version, which feeds Hamiltonian and SampOperators
% into a Sampler object. Ansatz contains all the information relevant to
% the variational wavefunction, and both contain references to the Hilbert
% space of configurations that they occupy.

% MultiChainSample is the multithreaded counterpart of EvalSample, intended
% for use with already optimised Ansatz objects rather than during
% optimisation (which implements a one-to-one core to chain ratio in
% multithreaded versions).

if Ncore == 1
    error('For single thread sampling, use EvalSample.')
end

Nsamp0 = SamplerObj.Nsamp; Nsamp = round(Nsamp0/Ncore,-1); % Round to nearest ten.

SamplerObj = SetNsamp(SamplerObj,Nsamp);

AnsatzStack = AnsatzObj; SamplerStack = SamplerObj;
for n = 1:Ncore
    AnsatzStack(n) = AnsatzObj; SamplerStack(n) = SamplerObj;
end

% Initialise averages of sampled quantities:
EnAvgStack = zeros(1,Ncore); NumEvals = numel(SamplerObj.Operators); % Number of additional quantities for averaging.
EvalAvgStack = cell(NumEvals,Ncore); % Create a cell array for storing the additional averaged quantities.
EnSamp = zeros(1,SamplerObj.Nsamp,Ncore); % Lists local energy sample values for each chain.
EvalSamp = cell(NumEvals,SamplerObj.Nsamp,Ncore); % Cell array for storing samples of the evaluated quantities.

MMCpool = parpool(Ncore); % Assigning object handle for later deletion / closure.

parfor c = 1:Ncore
    [EnAvgStack(c),EnSamp(:,:,c),EvalAvgStack(:,c),EvalSamp(:,:,c)] = ...
        EvalSample(SamplerStack(c),AnsatzStack(c));
end

delete(MMCpool);

EnAvg = mean(EnAvgStack); EvalAvg = cell(NumEvals,1);
for n = 1:NumEvals
    EvalAvg{n} = EvalAvgStack{n,1}/numel(EnAvgStack);
    for c = 2:numel(EnAvgStack)
        EvalAvg{n} = EvalAvg{n} + EvalAvgStack{n,c}/numel(EnAvgStack);
    end
end
end