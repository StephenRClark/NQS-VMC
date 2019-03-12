classdef NQSNHTI < NQSNH
    % NQSNHTI - a NQS Modifier variant that incorporates translation
    % invariance using a provided Graph with number-like hidden units and
    % square biases A and B.
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
    % Properties added with translation invariance:
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (1 x 1) for d/da.
    % - (1 x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha x 1) for d/dB.
    % - (Alpha*Nv x 1) for d/dWv.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        ati = 0; % Visible site bias, 1 x 1 scalar.
        Ati = 0; % Visible site square bias, 1 x 1 scalar.
        bti = 0; % Hidden site biases, Alpha x 1 vector.
        Bti = 0; % Hidden site square bias, 1 x 1 scalar.
        Wv = 0; % Hidden-visible couplings, Alpha x Nv matrix.
        Alpha = 1; % Hidden unit density / number of unique sets of W couplings.
    end
    
    methods
        
        % Constructor for 1D translation invariant NQS:
        function obj = NQSNHTI(Hilbert,Graph,Params,VFlag)
            obj@NQSNH(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiNQSNHTI(obj,Graph,Params);
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,Graph,dP)
            obj = PsiUpdateNQSNHTI(obj,Graph,dP);
        end
        
        % PsiCfgUpdate inherited from NQSNH.
        
        % PrepPsi inherited from NQSNH.
    end
    
    methods (Static)
        % PsiRatio inherited from NQSNH.
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,Graph,Cfg)
            [dLogp] = LogDerivNQSNHTI(obj,Hilbert,Graph,Cfg);
        end
    end
    
end