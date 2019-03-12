classdef StochasticReconfig < Optimiser
    % StochasticReconfig - an Optimiser subclass that is used to minimise
    % the variational energy of a given Ansatz for the Hamiltonian provided
    % in a Sampler object through the Stochastic Reconfiguration method.
    %   Optimiser is overarching class. StochasticReconfig features single
    %   and multi-thread options.
    
    properties 
        lambda0 = 10000; % Regularisation parameter - initial regularisation strength.
        lambda_min = 0.001; % Regularisation parameter - minimum regularisation strength.
        b = 0.95; % Regularisation parameter - regularisation decay by iteration.        
    end
    
    properties (Hidden) % Generally do not need to modify these, as these are just numerical checks.
        ETol = 5; % Maximum energy magnitude tolerance - physics dependent, used as an error checker.
        LRate = 1; % Learning rate - tunes the parameter change magnitude.
        PSave % Error correction parameter - number of steps to step back in event of error.
        PShift % Error correction parameter - shift regularisation by b^(-PShift) in the event of error.
        STol = 1e-30; % Tolerance in overlap matrix element magnitude - anything below is ignored.
        MRTol = 1e-5; % Move rate tolerance - if Metropolis move acceptance rate falls below, flag as potential error.
        ESens = 1e-12; % Minimum energy change tolerance - used to check if solution is stuck in a local minimum.
    end
    
    methods
        % General constructor for subclass - additional arguments in ETols and SRParams.
        function obj = StochasticReconfig(Npass,Ncore)
            obj@Optimiser(Npass,Ncore);
            obj.PSave = round(obj.Npass/20); obj.PShift = 2*obj.PSave;
        end
        
        % Reassign regularisation parameters if defaults are unsuitable.
        function [obj] = SetRegularisation(obj,lambda0_p,lambda_min_p,b_p)
            obj.lambda0 = lambda0_p; obj.lambda_min = lambda_min_p; obj.b = b_p;
        end
        
        % Reassign learning rate if default is unsuitable.
        function [obj] = SetLearnRate(obj,LRate_p)
            obj.LRate = LRate_p;
        end
        
        % Reassign energy and move rate tolerances if defaults are
        % unsuitable.
        function [obj] = SetSRTolerances(obj,ETol_p,MRTol_p)
            obj.ETol = ETol_p; obj.MRTol = MRTol_p;            
        end
        
        % Reassign reroll parameters if defaults are unsuitable.
        function [obj] = SetRollback(obj,PSave_p,PShift_p)
            obj.PSave = PSave_p; obj.PShift = PShift_p;
        end
    end
    
    methods
        function [Ansatz,EvalIter] = Optimise(obj,Sampler,Ansatz)
            [Ansatz,EvalIter] = SROptimise(obj,Sampler,Ansatz);            
        end
    end
    
end

