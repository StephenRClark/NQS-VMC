classdef Reference
    % Reference - a reference state that gives base amplitudes within the
    % Hilbert space to be modified with projectors.
    %   Main types of references - NoRF (none) / Bose (Bose condensate) /
    %   SDet (Slater determinant) / Pfaf (pairing amplitude Pfaffian). Only
    %   SDet and Pfaf offer variational versions.
    
    properties (Abstract)
        Type % Identifier for the reference state.
        VFlag % Set to 1 if the reference is meant to be variational.
        Np % All References need Np for NpTotal tracking purposes.
    end
    
    methods (Abstract)
        % Ratio between two configurations differing by Diff.
        [Ratio,Update] = PsiRatio(obj,Diff)
        
        % Update Reference configuration information according to Update.
        [Reference] = PsiCfgUpdate(obj,Update)
        
        % Initialise Ansatz configuration values given a starting Cfg.
        [Reference] = PrepPsi(obj,Cfg)
    end
    
end

