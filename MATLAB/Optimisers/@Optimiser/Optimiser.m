classdef Optimiser
    % Optimiser - given an Ansatz and Sampler, the Optimiser will change
    % the parameters of Ansatz via Monte Carlo sampling mediated by the
    % Sampler according to the action of the Hamiltonian in Sampler.
    %    Optimiser is overarching class - current subclasses include
    %    TimeVariational and StochasticReconfig.
    
    properties (SetAccess = protected)
        Npass % Number of optimisation passes made by Optimiser.
        Ncore % Number of threads / Markov chains the sampling is to be split across.
    end
    
    properties (Hidden, SetAccess = protected)
        ExtraSamp = 5000; % Extra samples in case of an error.
        dERTol = 0.05; % Maximum change in energy ratio tolerance - used as an error checker.
        dEVTol = 0.25; % Maximum change in energy value tolerance - used as an error checker.
        dPTol = 1e-15; % Minimum parameter change magnitude tolerance - anything below is ignored.
        ParamTol = 1e-15; % Minimum parameter magnitude tolerance - anything below is ignored.
        OFrac = exp(-1); % Default parameter fraction for random batch optimisation.
    end
    
    methods
        % General constructor for the class - simply assigns properties
        % inputted.
        function obj = Optimiser(Npass,Ncore)
            obj.Npass = Npass; obj.Ncore = Ncore;
        end
        
        % SetEnergyTolerances: Reassign dEVTol and dERTol if defaults are
        % unsuitable.
        function [obj] = SetEnergyTolerances(obj,dEV_p,dER_p)
            obj.dEVTol = dEV_p; obj.dERTol = dER_p;
        end
        
        % SetParamTolerances: Reassign dPTol and ParamTol if defaults are
        % unsuitable.
        function [obj] = SetParamTolerances(obj,dP_p,P_p)
            obj.dPTol = dP_p; obj.ParamTol = P_p;
        end
        
        % SetExtraSamples: Reassign extra samples if default is unsuitable.
        function [obj] = SetExtraSamples(obj,ExtraSamp_p)
            obj.ExtraSamp = ExtraSamp_p;
        end
        
        % SetBatchFraction: Reassign random batch optimisation fraction if
        % default is unsuitable.
        function [obj] = SetBatchFraction(obj,OFrac_p)
            if OFrac_p <= 0 || OFrac_p > 1
                error('Random batch fraction should be value between 0 and 1 (1 inclusive)');
            else
                obj.OFrac = OFrac_p;
            end
        end
    end    
    
    methods (Abstract)
        [Ansatz,EvalIter] = Optimise(obj,Sampler,Ansatz)
    end
    
end