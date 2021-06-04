classdef (Abstract) Modifier
    % Modifier - an object containing parameters to modify configuration
    % amplitudes for variational Monte Carlo. Plugs into Ansatz object.
    %   Modifier is the overarching class - subclasses will define specific
    %   functions for the methods listed below.    
    
    properties (Abstract)
        VFlag % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties (Abstract, Hidden)
        OptInds % Individual parameter flags for variational purposes.
    end
    
    properties (Abstract, SetAccess = protected)
        Np % Number of parameters.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Abstract, Hidden, SetAccess = protected)
        FullCfg % Function used by Modifer to interface with Cfg structs.
    end
    
    methods (Abstract)
        % PsiUpdate: Update Modifier variational parameters according to changes dP.
        [obj] = PsiUpdate(obj,dP);
        
        % PsiCfgUpdate: Update Modifier configuration information inside Update.
        [obj] = PsiCfgUpdate(obj,Update);
        
        % PrepPsi: Initialise Modifier configuration information given a starting Cfg.
        [obj] = PrepPsi(obj,Cfg);
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by Diff.
        [Ratio,Update] = PsiRatio(obj,Diff);
        
        % LogDeriv: Logarithmic derivative for the variational parameters in Modifier.
        [dLogp] = LogDeriv(obj,Cfg);
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        [Params] = ParamList(obj);
    end
    
    methods
        % RndBatchSelect: Randomly select some proportion of parameters to
        % optimise and disable optimisation of the rest. Used for random
        % batch optimisation schemes.
        function [obj] = RndBatchSelect(obj,OFrac)
            NpR = round(min(max(1,OFrac*obj.Np),obj.Np)); % Ensure at least 1 or at most Np are selected.
            OptIndsP = zeros(obj.Np,1); OptIndsP(randperm(obj.Np,NpR)) = 1;
            obj.OptInds = OptIndsP;
        end
        
    end

end

