classdef Jast < Modifier
    % Jast - a Modifier subclass that modifies configuration amplitudes
    % using Jastrow density-density or spin-spin correlation factors.
    %   Modifier is the overarching class.
    
    % ---------------------------------
    % Format for Jastrow Modifier object:
    % - Jast.N = number of sites (defined on input).
    % - Jast.Np = number of variational Jastrow parameters.
    % - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
    % - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
    % - Jast.Tj = (N x 1) vector - used to track on-site contributions.
    % - Jast.JsV = (N x N) matrix - contains variational parameter indices for each site.
    % ---------------------------------
    % Format for Update is a vector of new local Jastrow contributions.
    % ---------------------------------
    % Format for dLogp vector is a (Np x 1) vector of relevant two-site terms.
    % ---------------------------------
    % N.B: Jastrow factors here are assumed symmetric i.e.
    % Js(i,j) = JS(j,i).
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties (SetAccess = protected)
        Np = 1; % Number of parameters.
        N = 1; % Number of sites.
        Js = 0; % Placeholder for matrix of Jastrow factors arranged by site.
        JsVar = 0; % Placeholder for vector of Jastrow parameters.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        ParamCap = 10; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(1,1); % Individual parameter flags for variational purposes.
        NormFlag = 1; % Normalisation flag - sets the sum of Js to zero if on.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullBoseCfg; % Function used by Modifer to interface with Cfg structs.
        Tj = 0; % N x 1 vector of local density contributions to overall Jastrow factor.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
        JsV = ones(1); % Matrix of parameter indices for easy reconstruction.
    end
    
    methods
        % Constructor for Jast Modifier:
        function obj = Jast(Hilbert,Graph,Params,VFlag)
            %  Graph can be set to local (LVecs = 0) to give a Jastrow
            %  ansatz without translation invariance.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            obj.N = Hilbert.N;
            if strcmp(Hilbert.Type,'Ferm')
                obj.FFlag = Hilbert.N;
                % For Jastrow, fermionic systems require conversion of Diff
                % (performed before Ratio).
                obj.FullCfg = @FullFermDen; % Only sensitive to density, not spin.
            else
                obj.FFlag = 0; obj.FullCfg = @FullBoseCfg;
            end
            if isempty(Hilbert.Sector)
                obj.NormFlag = 0;
            end
            obj.N = Hilbert.N; obj.Graph = Graph;
            obj = RandomInitPsiJast(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateJast(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Tj = Update;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiJast(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised Jastrow amplitudes for a
        % given set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateJast(obj,Basis);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            if obj.FFlag ~= 0 % Fermionic Diff conversion required.
                [Diff] = Ferm2BoseDiff(Diff);
                % Ferm2BoseDiff converts fermionic Differences to density
                % differences.
            end
            Update = obj.Tj;
            if Diff.num ~= 0
                % For multiple changes in on-site quantum numbers, local dependency is
                % handled by Jast.Tj terms, and 'leftover' terms are contained in the
                % following:
                % Reshape to ensure no issues with different Diff formats.
                DiffVal = reshape(Diff.val,Diff.num,1);
                DiffPos = reshape(Diff.pos,Diff.num,1);
                % Requires the number of elements in DiffVal and DiffPos to be the same -
                % for fermionic treatments, some pre-processing of Diff may be required.
                JsRem = obj.Js(DiffPos,DiffPos) .* (DiffVal * DiffVal');
                Ratio = exp(-sum(JsRem(:))/2);
                for d = 1:Diff.num
                    % Diff.num should reflect number of sites that change local quantum numbers.
                    Ratio = Ratio * exp(-obj.Tj(Diff.pos(d)) * Diff.val(d));
                    for n = 1:obj.N
                        Update(n) = Update(n) + Diff.val(d) * obj.Js(Diff.pos(d),n);
                    end
                end
            else
                Ratio = 1;
            end
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            Cfg_vec = obj.FullCfg(Cfg);          
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            % Parameter indices associated with site pairs contained in JsV.
            for p = 1:obj.Np
                % Collect all two-site contributions related to parameter p.
                DDMatP = - ((obj.JsV == p) .* (Cfg_vec * Cfg_vec'))/2;
                dLogp(p) = sum(DDMatP(:));
            end            
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end 
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = obj.JsVar;
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadJast(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'Jast';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.N = obj.N; Properties.JsV = obj.JsV; 
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
        
    end
    
    
end