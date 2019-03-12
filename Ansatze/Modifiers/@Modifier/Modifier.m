classdef (Abstract) Modifier
    % Modifier - an object containing parameters to modify configuration
    % amplitudes for variational Monte Carlo. Plugs into Ansatz object.
    %   Modifier is the overarching class - subclasses will define specific
    %   functions for the methods listed below.    
    
    properties (Abstract)
        VFlag % Flag for whether to vary the parameters specified in this modifier.
        Np % Number of parameters.
        OptInds % Individual parameter flags for variational purposes.
    end
    
    methods (Abstract)
        % Update Modifier variational parameters according to changes dP.
        [Modifier] = PsiUpdate(Modifier,dP);
        
        % Update Modifier configuration information inside Update.
        [Modifier] = PsiCfgUpdate(Modifier,Update);
        
        % Initialise Modifier configuration information given a starting Cfg.
        [Modifier] = PrepPsi(Modifier,Cfg);
        
    end
    
    methods (Abstract, Static)
        % Ratio of amplitudes for two configurations separated by Diff.
        [Ratio,Update] = PsiRatio(Modifier,Diff);
        
        % Logarithmic derivative for the variational parameters in Modifier.
        [dLogp] = LogDeriv(Modifier,Cfg);
    end
    
    
end

