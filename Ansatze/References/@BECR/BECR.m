classdef BECR < Reference
    % BECR - a Bose condensate reference state that places all the
    % particles in the lowest eigenstate of the provided single-particle
    % Hamiltonian.
    %   Reference is overarching class. Currently no variational Bose
    %   variants are available.
    
    % ---------------------------------
    % Format for Bose Einstein Condensate Reference:
    % - BECR.Nb = number of bosons - assumed fixed for the most part.
    % - BECR.SPO = (Nv x 1) vector - single particle boson orbital amplitudes.
    % - BECR.Occ = (Nv x 1) vector - onsite boson occupation numbers.
    % ---------------------------------
    % Format for Update is a vector containing the new boson occupation numbers.
    % ---------------------------------    
    
    properties
        Type = 'BECR'; % Identifier for the reference state.
        VFlag = 0; % Variational optimisation of Bose is disabled.
        Nb = 1; % Number of bosons, inherited from input Hilbert.
        Np = 0; % Currently no variational parameters associated with BECR.
    end
    
    properties (Hidden)
        Occ = 1; % Vector of boson number on each site.
        SPO = 1; % Single particle orbital.
    end
    
    methods
        % Constructor for BECR Reference object:
        function obj = BECR(Hilbert,Params)
            obj.Nb = Hilbert.Nb; SPH = Params.SPH;
            % Calculate lowest energy single particle orbital from single
            % particle Hamiltonian in Params.
            [SPO,~] = eigs(SPH,1,'sa');
            obj.SPO = SPO;
            % Placeholder for Occ, which is initialised later.
            obj.Occ = (Hilbert.Nb / Hilbert.N) * ones(1,Hilbert.N);
        end
        
        % Initialise Reference configuration values given a starting Cfg.
        function [obj] = PrepPsi(obj,~,Cfg)
            obj = PrepPsiBECR(obj,Cfg);
        end
        
        % Update Reference configuration information according to Update.
        function [obj] = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdateBECR(obj,Update);
        end
    end
    
    methods (Static)
        % Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = PsiRatioBECR(obj,Diff);
        end
    end
    
end

