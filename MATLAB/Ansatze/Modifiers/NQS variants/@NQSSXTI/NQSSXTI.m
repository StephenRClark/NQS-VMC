classdef NQSSXTI < NQSSX
    % NQSSXTI - a NQS Modifier variant that incorporates translation
    % invariance using a provided Graph with a square-square interaction
    % term X, multinomial visible and hidden units.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQS Modifier object with square-square interaction:
    % - NQSSX.Nv = number of "visible" spins.
    % - NQSSX.Nh = number of "hidden" spins.
    % - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
    % - NQSSX.a = (Nv x 1) vector - visible site bias.
    % - NQSSX.A = (Nv x 1) vector - visible site square bias.
    % - NQSSX.b = (Nh x 1) vector - hidden site bias.
    % - NQSSX.B = (Nh x 1) vector - hidden site square bias.
    % - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
    % - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
    % - NQSSX.HDim = dimension of the hidden units.
    % - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
    % - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
    % - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
    % - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
    % - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % Properties added with translation invariance:
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Wv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Xv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.VisVec - vector of new visible occupancies.
    % Update.ThetaSq - vector of new effective angles ThetaSqP.
    % Update.NsqVec - vector of new squared visible occupancies.
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
        bti = 0; % Hidden site biases, Alpha x 1 vector.
        Bti = 0; % Hidden site square bias, Alpha x 1 vector.
        Wv = 0; % Hidden-visible couplings, Alpha x Nv matrix.
        Xv = 0; % Hidden-visible square couplings, Alpha x Nv matrix.
        Alpha = 1; % Hidden unit density / number of unique sets of couplings.
    end
    
    methods
        % Constructor for number-hidden translation invariant NQS:
        function obj = NQSSXTI(Hilbert,Graph,Params,VFlag)
            obj@NQSSX(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiNQSSXTI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSSXTI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from NQSSH.
        
        % PrepPsi inherited from NQSSH.
        
        % PsiGenerate: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSSXTI(obj,Params);
        end
        
        % PsiRatio inherited from NQSNH.
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            % Make local copies to reduce notation in code below.
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
            % translates by some combination of Graph.Lvecs.
            Ng = GraphObj.N; % Number of actual sites in Graph - necessary if NQS uses enlarged lattice.
            Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.            
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            Cfg_sqr = Cfg_vec.^2;
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            if obj.OptInds(1) == 1
                dLogp(1) = sum(Cfg_vec); % Insert d/da.
            end
            if obj.OptInds(2) == 1
                dLogp(2) = sum(Cfg_sqr); % Insert d/dA.
            end
            dTheta = dL_SqTrace(obj.Theta,obj.ThetaSq,obj.HVal);
            dThetaSq = dS_SqTrace(obj.Theta,obj.ThetaSq,obj.HVal);
            % Accounting for shift structure of W matrix requires either construction
            % of shifted Theta matrix or shifted Cfg vector - the latter is done here
            for a=1:obj.Alpha % Derivatives need to be computed by Alpha sector
                if obj.OptInds(2+a) == 1
                    dLogp(2+a) = sum(dTheta((1:Ntr)+(a-1)*Ntr)); % Insert d/db.
                end
                if obj.OptInds(2+obj.Alpha+a) == 1
                    dLogp(2+obj.Alpha+a) = sum(dThetaSq((1:Ntr)+(a-1)*Ntr)); % Insert d/dB.
                end
                for v = 1:obj.Nv
                    PIndW = 2 + 2*obj.Alpha + v + (a-1)*obj.Nv;
                    PIndX = PIndW + obj.Alpha*obj.Nv;
                    % For each layer labelled by a, find the indices of the associated translates.
                    for b = 1:Ntr
                        TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
                        if VInd ~= 0
                            if obj.OptInds(PIndW) == 1
                                dLogp(PIndW) = dLogp(PIndW) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dw.
                            end
                            if obj.OptInds(PIndX) == 1
                                dLogp(PIndX) = dLogp(PIndX) + (Cfg_sqr(VInd)*dThetaSq(TInd)); % Insert d/dW.
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
            Params = ParamListNQSSXTI(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSSXTI(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSSXTI';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Alpha = obj.Alpha; Properties.HDim = obj.HDim;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end