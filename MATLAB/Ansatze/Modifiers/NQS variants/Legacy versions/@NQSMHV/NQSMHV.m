classdef NQSMHV < NQSMH
    % NQSMHV - a NQSMH Modifier variant that incorporates translation
    % invariance using the provided graph. Variant that fixes some long
    % range W/X to be equal, allowing close range values to vary.
    %   NQSMH is overarching class, which is itself a subclass of Modifier.
    %   While related to NQS, enough features differ that NQSMH cannot be
    %   instanced as a subclass of NQS.
    
    % ---------------------------------
    % Format for NQS Modifier object with multiplon-holon interactions:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
    % - NQS.A = (Nv x 1) vector - visible site square bias.
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.BH = (Nh x 1) vector - hidden holon bias.
    % - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
    % - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
    % - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
    % - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
    % - NQS.Hv = (Nv x 1) vector - vector of visible holons.
    % - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
    % - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
    % Properties added with translation invariance:
    % - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.BHti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.BMti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % ---------------------------------
    % Format for Update is a struct with five fields:
    % Update.ThetaH - vector of new effective angles ThetaHP.
    % Update.ThetaM - vector of new effective angles ThetaMP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % Update.Hv - vector of new holon operator values HvP.
    % Update.Mv - vector of new multiplon operator values MvP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (1 x 1) for d/da.
    % - (1 x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha x 1) for d/dB.
    % - (Alpha*Nv x 1) for d/dW.
    % - (Alpha*Nv x 1) for d/dX.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        ati = 0; % Visible site bias, 1 x 1 scalar.
        Ati = 0; % Visible site square bias, 1 x 1 scalar.
        BHti = 0; % Hidden site holon bias, Alpha x 1 vector.
        BMti = 0; % Hidden site multiplon bias, Alpha x 1 vector.
        Wv = 0; % Hidden-visible MM/HH couplings, Alpha x Nv matrix.
        Xv = 0; % Hidden-visible MH/HM couplings, Alpha x Nv matrix.
        Alpha = 1; % Hidden unit density / number of unique sets of couplings.
    end
    
    methods
        % Constructor for multiplon-holon NQS with translation invariance:
        function obj = NQSMHV(Hilbert,Graph,Params,VFlag)
            obj@NQSMH(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiNQSMHV(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSMHV(obj,dP);
        end
        
        % PsiCfgUpdate: inherited from NQSMH.
        
        % PrepPsi: inherited from NQSMH.
        
        % PsiGenerate: inherited from NQSMH.
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSMHV(obj,Params);
        end
        
        % PsiRatio: inherited from NQSMH.
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = LogDerivNQSMHV(obj,Cfg);
        end
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            [Params] = ParamListNQSMHV(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSMHV(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSMHV';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds; 
            Properties.Nv = obj.Nv; Properties.Alpha = obj.Alpha; Properties.HDim = obj.HDim;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end