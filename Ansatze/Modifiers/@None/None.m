classdef None < Modifier
    % None - a Modifier object that leaves configuration amplitudes
    % untouched. The Modifier equivalent of the Plus subclass of
    % Reference.
    %   Modifier is overarching class.    
    
    properties 
        VFlag = 0; % Nothing to vary here.
        Np = 0; % No parameters.
        OptInds = 0; % No parameters to flag.
    end
    
    methods
        % Constructor for None Modifier object:
        function [obj] = None()
            obj.VFlag = 0; obj.Np = 0; obj.OptInds = 0;
        end
        
        % Update Modifier variational parameters according to changes dP.
        function [obj] = PsiUpdate(obj,~)
            % Don't need to do anything here, really. Just need the function to exist.
        end
        
        % Update Modifier configuration information inside Update.
        function [obj] = PsiCfgUpdate(obj,~)
            % Don't need to do anything here, really. Just need the function to exist.
        end
        % Initialise Modifier configuration information given a starting Cfg.
        function [obj] = PrepPsi(obj,~)
            % Don't need to do anything here, really. Just need the function to exist.
        end
    end
    
    methods (Static)
        % Ratio of amplitudes for two configurations separated by Diff.
        function [Ratio,Update] = PsiRatio(~,~)
            Ratio = 1; Update = 0; % Need these outputs to be defined.
        end
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(~,~)
            dLogp = 0; % Need this output to be defined. 
            % Function won't ever be called, but needed for compatability.
        end
    end
    
    
end

