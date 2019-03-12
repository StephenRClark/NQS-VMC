classdef JastTI < Jast
    % JastTI - a Jastrow Modifier variant that incorporates translation
    % invariance using a provided Graph.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for Jastrow Modifier object:
    % - Jast.N = number of sites (defined on input).
    % - Jast.Np = number of variational Jastrow parameters.
    % - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
    % - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
    % - Jast.Tj = (N x 1) vector - used to track on-site contributions.
    % Properties added with translation invariance:
    % - Jast.JsV = (N x N) matrix - contains variational parameter indices for each site.
    % ---------------------------------
    % Format for Update is a vector of new local Jastrow contributions.
    % ---------------------------------
    % Format for dLogp vector is a (Np x 1) vector of relevant two-site terms.
    % ---------------------------------
    % N.B: Jastrow factors here are assumed symmetric i.e.
    % Js(i,j) = JS(j,i).
    % ---------------------------------
    
    properties (Hidden) % All other properties inherited from Jast.
        JsV = 0; % N x N matrix marking positions of each of the Np variational terms.
    end
    
    methods
        % Constructor for unconstrained Jast with no symmetries:
        function obj = JastTI(Hilbert,Graph,Params,VFlag)
            %  Graph necessary for Jast subvariants as second argument.
            obj@Jast(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiJastTI(obj,Graph,Params);
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,~,dP)
            obj = PsiUpdateJastTI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from Jast.
        
        % PrepPsi inherited from Jast.        
    end
    
    methods (Static)        
        % PsiRatio inherited from Jast.
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,~,Cfg)
            [dLogp] = LogDerivJastTI(obj,Hilbert,Cfg);
        end
    end
    
    
end