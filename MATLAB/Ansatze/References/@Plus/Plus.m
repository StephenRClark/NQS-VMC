classdef Plus < Reference
    % Plus - a reference state that is an equal superposition of all
    % possible states in the Hilbert space. 
    %   Reference is overarching class. Plus cannot be made variational.
    
    % ---------------------------------
    % As an effective 'blank slate' to start with, Plus has no real
    % structure or format, as it does not alter any amplitudes.
    % ---------------------------------
    
    properties         
        VFlag = 0; % Set to 0 as Plus cannot be made variational.        
    end
    
    properties (SetAccess = protected)
        Type = 'Plus'; % Identifier for the reference state.
        Np = 0; % No variational parameters associated with Plus.        
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = []; % No need to carry out any Cfg conversions.
    end
    
    methods 
        function obj = Plus()
            obj.VFlag = 0; obj.Type = 'Plus';
        end
        
        % Update Ansatz configuration information according to Update.
        function [obj] = PsiCfgUpdate(obj,~)
            0;
            % Don't need to do anything. Later version should include some
            % way to avoid the function call, but initial version will have
            % this empty function called each time rather than 'if' checks.
        end
        
        % Initialise Ansatz configuration values given a starting Cfg.
        function [obj] = PrepPsi(obj,~)
            0;
            % Don't need to do anything. Later version should include some
            % way to avoid the function call, but initial version will have
            % this empty function called each time rather than 'if' checks.
        end
   
        % Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,~)
           Ratio = 1; Update = [];
           % Don't need to do anything. Later version should include some
           % way to avoid the function call, but initial version will have
           % this empty function called each time rather than 'if' checks.
        end

        % PsiGenerate is just a flat state of equal amplitudes:
        function [Psi] = PsiGenerate(obj,Basis)
            N_cfgs = size(Basis,1); Psi = ones(N_cfgs)/sqrt(N_cfgs);
        end
    end
    
end

