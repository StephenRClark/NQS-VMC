classdef NQSOHTI < NQSOH
    % NQSTI - a NQSOH Modifier variant that incorporates translation
    % invariance using a provided Graph.
    %   NQSOH is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQSOH Modifier object:
    % - NQSOH.Nv = number of "visible" spins.
    % - NQSOH.Nh = number of "hidden" spins.
    % - NQSOH.Np = number of parameters in the ansatz = VDim + Alpha + (VDim*Nv * Alpha).
    % - NQSOH.VDim = dimensions of the visible units.
    % - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
    % - NQSOH.b = (Nh x 1) vector - hidden site bias.
    % - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
    % - NQSOH.Theta = (Nh x 1) vector - effective angles.
    % - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
    % Properties added with translation invariance:
    % - NQS.ati = (VDim x 1) vector - reduced parameter set for TI.
    % - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Wv = (Alpha x VDim*Nv) matrix - reduced parameter set for TI.
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP and
    % new one-hot vector OHVecP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (VDim x 1) for d/da.
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv*VDim x 1) for d/dW.
    % Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        ati = 0; % Visible site bias, 1 x 1 scalar.
        bti = 0; % Hidden site biases, Alpha x 1 vector.
        Wv = 0; % Hidden-visible couplings, Alpha x Nv matrix.
        Alpha = 1; % Hidden unit density / number of unique sets of W couplings.
    end
    
    methods
        
        % Constructor for 1D translation invariant NQS:
        function obj = NQSOHTI(Hilbert,Graph,Params,VFlag)
            obj@NQSOH(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiNQSOHTI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSOHTI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from NQSOH.
        
        % PrepPsi inherited from NQSOH.
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSOHTI(obj,Params);
        end
   
        % PsiRatio inherited from NQSOH.
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = LogDerivNQSOHTI(obj,Cfg);
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSOHTI(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSOHTI(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSOHTI';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Alpha = obj.Alpha; 
            Properties.VDim = obj.VDim; Properties.VList = obj.VList;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end