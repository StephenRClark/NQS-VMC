classdef NQSTISS < NQS
    % NQSTISS - a NQS Modifier variant that incorporates translation
    % invariance using a provided Graph as well as spin symmetry. For the
    % same parameter number characterised by Alpha,
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQS Modifier object with translation invariance:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = Nh + Alpha + 1. (computed here).
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.b = (Nh x 1) vector - hidden site bias.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % Properties added with translation invariance:
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (1 x 1) for d/da.
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv x 1) for d/dWv.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        ati = 0; % Visible site bias, 1 x 1 scalar.
        bti = 0; % Hidden site biases, Alpha x 1 vector.
        Wv = 0; % Hidden-visible couplings, Alpha x Nv matrix.
        Alpha = 1; % Hidden unit density / number of unique sets of W couplings.
    end
    
    methods
        
        % Constructor for 1D translation invariant NQS:
        function obj = NQSTISS(Hilbert,Graph,Params,VFlag)
            obj@NQS(Hilbert,Graph,Params,VFlag);
            if obj.SFlag == 0
                error('Provided Hilbert object does not describe a spin-symmetric system. Try NQSTI as an alternative.');
            end
            obj = RandomInitPsiNQSTISS(obj,Graph,Params);
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,Graph,dP)
            obj = PsiUpdateNQSTISS(obj,Graph,dP);
        end
        
        % PsiCfgUpdate inherited from NQS.
        
        % PrepPsi inherited from NQS.
    end
    
    methods (Static)
        % PsiRatio inherited from NQS.
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,Graph,Cfg)
            [dLogp] = LogDerivNQSTISS(obj,Hilbert,Graph,Cfg);
        end
    end
    
end