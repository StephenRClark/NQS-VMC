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
        VFlag = 0; % Variational optimisation of Bose is disabled.
    end
    
    properties (SetAccess = protected)
        Type = 'Bose'; % Identifier for the reference state.
        Nb = 1; % Number of bosons, inherited from input Hilbert.
        Np = 0; % Currently no variational parameters associated with BECR.
        Graph % Details connectivity of sites - used to include symmetries.
    end
    
    properties (Hidden, SetAccess = protected)
        Occ = 1; % Vector of boson number on each site.
        SPO = 1; % Single particle orbital.
        FullCfg = @FullBoseCfg; % Function used by Reference to interface with Cfg structs.
    end
    
    methods
        % Constructor for BECR Reference object:
        function obj = BECR(Hilbert,Graph,Params)
            obj.Nb = Hilbert.Nb; SPH = Params.SPH;
            % Calculate lowest energy single particle orbital from single
            % particle Hamiltonian in Params.
            [SPO,~] = eigs(SPH,1,'sa');
            obj.SPO = SPO;
            % Placeholder for Occ, which is initialised later.
            obj.Occ = (Hilbert.Nb / Hilbert.N) * ones(1,Hilbert.N);
            obj.Graph = Graph;
        end
        
        % PrepPsi: Initialise Reference configuration values given a
        % starting Cfg.
        function [obj] = PrepPsi(obj,Cfg)
            obj = PrepPsiBECR(obj,Cfg);
        end
        
        % PsiCfgUpdate: Update Reference configuration information
        % according to Update.
        function [obj] = PsiCfgUpdate(obj,Update)
            obj.Occ = Update.Occ; obj.Nb = Update.Nb;
        end
        
        % PsiGenerate: Generate full normalised amplitudes for a given set
        % of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateBECR(obj,Basis);
        end
        
        % PsiRatio: Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            Update.Occ = obj.Occ; Update.Nb = obj.Nb;
            Ratio = prod(obj.SPO(Diff.pos).'.^Diff.val);
            for d = 1:Diff.num                
                Update.Occ(Diff.pos(d)) = Update.Occ(Diff.pos(d)) + Diff.val(d);
                Update.Nb = Update.Nb + Diff.val(d);
                for m = 1:abs(Diff.val(d))
                    Ratio = Ratio * sqrt(m+min([Update.Occ(Diff.pos(d)),obj.Occ(Diff.pos(d))]))^(-sign(Diff.val(d)));
                end
            end
            for d = 1:abs(sum(Diff.val))
                Ratio = Ratio * sqrt(min([Update.Nb,obj.Nb]))^(-sign(sum(Diff.val(d))));
            end
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'BECR';
            Properties.Graph = obj.Graph.PropertyList; Properties.SPO = obj.SPO;
        end
    end
    
end

