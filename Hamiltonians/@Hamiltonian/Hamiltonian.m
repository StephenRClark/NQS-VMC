classdef Hamiltonian
    % Hamiltonian - an object containing the details of the Hamiltonian
    % being applied to the physical system.
    %   Hamiltonian objects are fed into Sampler objects for sampling with
    %   an Ansatz.
    
    properties (Abstract) 
        Hilbert % Contains configuration information for the Markov chain.
        Operator % List of Operator objects that compose the Hamiltonian.
        HParams % List of parameters that characterise the Hamiltonian.
        TFuncs % Cell list of time-dependent functions that modify the HParams. Can be left empty.
        % Note that each Operator needs its own Graph that specifies their
        % connectivity on the lattice.
    end
    
    properties (Hidden, Abstract)
        HParams0 % Original input Hamiltonian parameters - used as a reference when time evolving.
    end
    
    methods (Abstract)
        % Generate the matrix elements and the differences created by the
        % action of the Hamiltonian on a configuration Cfg.
        [Diff,HMatEls] = HamMatEls(obj,Cfg);
        
        % Sample the individual operators composing the Hamiltonian and
        % output the local energy EnLoc.
        [EnLoc] = EnergySample(obj,Cfg,Ansatz);
        
        % Modify HParams according to any TFuncs specified at
        % initialisation. Only called during time-dependent optimisations.
        [obj] = TimeEvolveH(obj,time);
    end    
end

