classdef Optimiser
    % Optimiser - given an Ansatz and Sampler, the Optimiser will change
    % the parameters of Ansatz via Monte Carlo sampling mediated by the
    % Sampler according to the action of the Hamiltonian in Sampler.
    %    Optimiser is overarching class - current subclasses include
    %    TimeVariational and StochasticReconfig.
    
    properties
        Npass % Number of optimisation passes made by Optimiser.
        Ncore % Number of threads / Markov chains the sampling is to be split across.
    end
    
    properties (Hidden)
        ExtraSamp = 5000; % Extra samples in case of an error.
        dERTol = 0.05; % Maximum change in energy ratio tolerance - used as an error checker.
        dEVTol = 0.25; % Maximum change in energy value tolerance - used as an error checker.
        dPTol = 1e-15; % Minimum parameter change magnitude tolerance - anything below is ignored.
        ParamTol = 1e-15; % Minimum parameter magnitude tolerance - anything below is ignored.
    end
    
    methods
        % General constructor for the class - simply assigns properties inputted.
        function obj = Optimiser(Npass,Ncore)
            obj.Npass = Npass; obj.Ncore = Ncore;
        end
        
        % Reassign dEVTol and dERTol if defaults are unsuitable.
        function [obj] = SetEnergyTolerances(obj,dEV_p,dER_p)
            obj.dEVTol = dEV_p; obj.dERTol = dER_p;
        end
        
        % Reassign extra samples if default is unsuitable.
        function [obj] = SetExtraSamples(obj,ExtraSamp_p)
            obj.ExtraSamp = ExtraSamp_p;
        end
    end
    
    
    methods (Abstract)
        [Ansatz,EvalIter] = Optimise(obj,Sampler,Ansatz)
    end
    
end

