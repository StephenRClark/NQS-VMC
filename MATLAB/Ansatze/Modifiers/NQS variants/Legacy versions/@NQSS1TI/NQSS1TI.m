classdef NQSS1TI < NQSS1
    % NQSS1 - a NQS Modifier variant that uses binary hidden units but
    % visible interaction terms suited for spin-1 systems. TI for
    % translation invariance.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % Format for NQS Modifier object modified for spin-1:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = 2*Nv*Alpha + 2 + Alpha.
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.A = (Nv x 1) vector - visible site square bias.
    % - NQS.b = (Nh x 1) vector - hidden site bias.
    % - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % Properties added with translation invariance:
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.wv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (1 x 1) for d/da.
    % - (1 x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv x 1) for d/dw.
    % - (Alpha*Nv x 1) for d/dW.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        ati = 0; % Visible site linear bias.
        Ati = 0; % Visible site square bias.
        bti = 0; % Hidden site linear bias.
        wv = 0; % Hidden-visible linear coupling.
        Wv = 0; % Hidden-visible square coupling.
        Alpha = 1; % Hidden unit density.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSS1TI(Hilbert,Graph,Params,VFlag)
            obj@NQSS1(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiNQSS1TI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSS1TI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from NQSS1.
        
        % PrepPsi inherited from NQSS1.
        
        % PsiGenerate inherited from NQSS1.
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSS1TI(obj,Params);
        end
        
        % PsiRatio inherited from NQSS1.
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
            % translates by some combination of Graph.Lvecs.
            Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
            Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            Cfg_sqr = Cfg_vec.^2;
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.
            dTheta = tanh(obj.Theta);
            if obj.OptInds(1) == 1
                dLogp(1) = sum(Cfg_vec); % Insert d/da.
            end
            if obj.OptInds(2) == 1
                dLogp(2) = sum(Cfg_sqr); % Insert d/dA.
            end
            for a=1:obj.Alpha % Derivatives need to be computed by Alpha sector
                if obj.OptInds(2+a) == 1
                    dLogp(2+a) = sum(dTheta((1:Ntr)+(a-1)*Ntr)); % Insert d/db.
                end
                for v = 1:obj.Nv
                    PIndw = 2 + obj.Alpha + v + (a-1)*obj.Nv;
                    PIndW = PIndw + obj.Alpha*obj.Nv;
                    % For each layer labelled by a, find the indices of the associated translates.
                    for b = 1:Ntr
                        TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
                        if VInd ~= 0
                            if obj.OptInds(PIndw) == 1
                                dLogp(PIndw) = dLogp(PIndw) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dw.
                            end
                            if obj.OptInds(PIndW) == 1
                                dLogp(PIndW) = dLogp(PIndW) + (Cfg_sqr(VInd)*dTheta(TInd)); % Insert d/dW.
                            end
                        end
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSS1TI(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSS1TI(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSS1TI';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Alpha = obj.Alpha;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end