classdef Reference
    % Reference - a reference state that gives base amplitudes within the
    % Hilbert space to be modified with projectors.
    %   Main types of references - NoRF (none) / Bose (Bose condensate) /
    %   SDet (Slater determinant) / Pfaf (pairing amplitude Pfaffian). Only
    %   SDet and Pfaf offer variational versions.
    
    properties (Abstract)
        VFlag % Set to 1 if the reference is meant to be variational.
    end
    
    properties (Abstract, SetAccess = protected)
        Type % Identifies Hilbert compatibility.
        Np % All References need Np for NpTotal tracking purposes.
    end
    
    properties (Abstract, Hidden, SetAccess = protected)
        FullCfg % Function used by Reference to interface with Cfg structs.
    end
    
    methods (Abstract)
        % PsiCfgUpdate: Update Reference configuration information inside Update.
        [obj] = PsiCfgUpdate(obj,Update);
        
        % PrepPsi: Initialise Reference configuration information given a starting Cfg.
        [obj] = PrepPsi(obj,Cfg);
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by Diff.
        [Ratio,Update] = PsiRatio(obj,Diff);
        
        % PsiUpdate, LogDeriv and RndBatchSelect are included case-by-case,
        % since not all References are allowed to be variational.
    end
    
end

