classdef NQSNH < NQS
    % NQSNH - a NQS Modifier variant that uses number-like hidden units
    % with dimension HDim instead of spin-like hidden units.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQS Modifier object with number hidden units:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = Nh + 2*Alpha + 2.
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.A = (Nv x 1) vector - visible site square bias.
    % - NQS.b = (Nh x 1) vector - hidden site bias.
    % - NQS.B = (Nh x 1) vector - hidden site square bias.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv x 1) for d/da.
    % - (Nv x 1) for d/dA.
    % - (Nh x 1) for d/db.
    % - (Nh x 1) for d/dB.
    % - (Nh*Nv x 1) for d/dW.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        A = 0; % Visible site square bias.
        B = 0; % Hidden site square bias.
        HDim = 2; % Hidden unit dimensionality, set to match Hilbert.
    end
    
    properties (Hidden)
        NsqVec = 0; % Squared visible occupancies, Nv x 1 vector.
    end
    
    methods
        
        % Constructor for 1D translation invariant NQS:
        function obj = NQSNH(Hilbert,Graph,Params,VFlag)
            obj@NQS(Hilbert,Graph,Params,VFlag);
            obj.HDim = Hilbert.d; % Set NQS hidden dimension to match visible.
            obj = RandomInitPsiNQSNH(obj,Params);
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,Graph,dP)
            obj = PsiUpdateNQSNH(obj,Graph,dP);
        end
        
        % Update Modifier configuration information inside Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdateNQSNH(obj,Update);
        end
        
        % Initialise Modifier configuration information given a starting Cfg.
        function obj = PrepPsi(obj,Hilbert,Cfg)
            obj = PrepPsiNQSNH(obj,Hilbert,Cfg);
        end
        
    end
    
    methods (Static)
        % Ratio of amplitudes for two configurations separated by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = PsiRatioNQSNH(obj,Diff);
        end
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,Graph,Cfg)
            [dLogp] = LogDerivNQSNH(obj,Hilbert,Graph,Cfg);
        end
    end
    
end