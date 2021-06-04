classdef None < Modifier
    % None - a Modifier object that leaves configuration amplitudes
    % untouched. The Modifier equivalent of the Plus subclass of
    % Reference.
    %   Modifier is overarching class.
    
    properties
        VFlag = 0; % Nothing to vary here.
    end
    
    properties (SetAccess = protected)
        Np = 0; % No parameters.
        Graph = []; % No need to worry about the lattice.
    end
    
    properties (Hidden)
        OptInds = 0; % No parameters to flag.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = []; % No need to interface with Cfg structs.
    end
    
    methods
        % Constructor for None Modifier object:
        function [obj] = None()
            obj.VFlag = 0; obj.Np = 0; obj.OptInds = 0;
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function [obj] = PsiUpdate(obj,~)
            % Don't need to do anything here, really. Just need the
            % function to exist.
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function [obj] = PsiCfgUpdate(obj,~)
            % Don't need to do anything here, really. Just need the
            % function to exist.
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function [obj] = PrepPsi(obj,~)
            % Don't need to do anything here, really. Just need the
            % function to exist.
        end
    
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(~,~)
            Ratio = 1; Update = 0; % Need these outputs to be defined.
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(~,~)
            dLogp = 0; % Need this output to be defined.
            % Function won't ever be called, but needed for compatability.
        end
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(~)
            Params = [];
        end
    end    
    
end