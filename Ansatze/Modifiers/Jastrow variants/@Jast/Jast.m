classdef Jast < Modifier
    % Jast - a Modifier subclass that modifies configuration amplitudes
    % using Jastrow density-density or spin-spin correlation factors.
    %   Modifier is the overarching class. Jast itself has subvariants with
    %   symmetries and projections built into them.
    
    % ---------------------------------
    % Format for Jastrow Modifier object:
    % - Jast.N = number of sites (defined on input).
    % - Jast.Np = number of variational Jastrow parameters.
    % - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
    % - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
    % - Jast.Tj = (N x 1) vector - used to track on-site contributions.
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
        Np = 1; % Number of parameters.
        N = 1; % Number of sites.
        Js = 0; % Placeholder for matrix of Jastrow factors arranged by site.
        JsVar = 0; % Placeholder for vector of Jastrow parameters.
    end
    
    properties (Hidden)
        Tj = 0; % N x 1 vector of local density contributions to overall Jastrow factor.
        ParamCap = 10; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(1,1); % Individual parameter flags for variational purposes.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
    end
    
    methods
        % Constructor for unconstrained Jast with no symmetries:
        function obj = Jast(Hilbert,~,Params,VFlag)
            %  Graph necessary for Jast subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            obj.N = Hilbert.N;
            if strcmp(Hilbert.Type,'Ferm')
                obj.FFlag = Hilbert.N;
                % For Jastrow, fermionic systems require conversion of Diff
                % (performed before Ratio) and conversion of Cfg_vec
                % (performed within relevan functions).
            else
                obj.FFlag = 0;
            end
            obj.N = Hilbert.N;
            obj = RandomInitPsiJast(obj,Params);
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,~,dP)
            obj = PsiUpdateJast(obj,dP);
        end
        
        % Update Modifier configuration information inside Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdateJast(obj,Update);
        end
        
        % Initialise Modifier configuration information given a starting Cfg.
        function obj = PrepPsi(obj,Hilbert,Cfg)
            obj = PrepPsiJast(obj,Hilbert,Cfg);
        end
        
    end
    
    methods (Static)   
        % Ratio of amplitudes for two configurations separated by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            if obj.FFlag ~= 0 % Fermionic Diff conversion required.
                [Diff] = Ferm2BoseDiff(Diff);
                % Ferm2BoseDiff converts fermionic Differences to density
                % differences.
            end
            [Ratio,Update] = PsiRatioJast(obj,Diff);
        end
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,~,Cfg)
            [dLogp] = LogDerivJast(obj,Hilbert,Cfg);
        end
    end
    
    
end