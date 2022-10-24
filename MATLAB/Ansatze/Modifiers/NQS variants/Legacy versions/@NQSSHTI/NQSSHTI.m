classdef NQSSHTI < NQSSH
    % NQSSHTI - a NQS Modifier variant that incorporates translation
    % invariance using a provided Graph with higher dimension spin-like
    % hidden units and square biases A and B.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQS Modifier object with high dimension hidden units:
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
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.bti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.Bti = (Alpha x 1) vector - reduced parameter set for TI.
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
    % - (Alpha x 1) for d/dB.
    % - (Alpha*Nv x 1) for d/dWv.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        ati = 0; % Visible site bias, 1 x 1 scalar.
        Ati = 0; % Visible site square bias, 1 x 1 scalar.
        bti = 0; % Hidden site biases, Alpha x 1 vector.
        Bti = 0; % Hidden site square bias, Alpha x 1 vector.
        Wv = 0; % Hidden-visible couplings, Alpha x Nv matrix.
        Alpha = 1; % Hidden unit density / number of unique sets of W couplings.
    end
    
    methods
        % Constructor for number-hidden translation invariant NQS:
        function obj = NQSSHTI(Hilbert,Graph,Params,VFlag)
            obj@NQSSH(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiNQSNHTI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSNHTI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from NQSSH.
        
        % PrepPsi inherited from NQSSH.
        
        % PsiGenerate: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSNHTI(obj,Params);
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
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            if obj.OptInds(1) == 1
                dLogp(1) = sum(Cfg_vec); % Insert d/da.
            end
            if obj.OptInds(2) == 1
                dLogp(2) = sum(Cfg_vec.^2); % Insert d/dA.
            end
            dTheta = dT_SHTrace(obj.Theta,obj.B,obj.HDim);
            dB = dB_SHTrace(obj.Theta,obj.B,obj.HDim);            
            % Accounting for shift structure of W matrix requires either construction
            % of shifted Theta matrix or shifted Cfg vector - the latter is done here
            for a=1:obj.Alpha % Derivatives need to be computed by Alpha sector
                if obj.OptInds(2+a) == 1
                    dLogp(2+a) = sum(dTheta((1:Ntr)+(a-1)*Ntr)); % Insert d/db.
                end
                if obj.OptInds(2+obj.Alpha+a) == 1
                    dLogp(2+obj.Alpha+a) = sum(dB((1:Ntr)+(a-1)*Ntr)); % Insert d/dB.
                end
                for v = 1:obj.Nv
                    PInd = 2 + 2*obj.Alpha + v + (a-1)*obj.Nv;
                    % For each layer labelled by a, find the indices of the associated translates.
                    if obj.OptInds(PInd) == 1
                        for b = 1:Ntr
                            TInd = b + (a-1)*Ntr; VInd = BondMap{b}(1+mod(v-1,Ng)) + Ng*(ceil(v/Ng)-1);
                            if VInd ~= 0
                                dLogp(PInd) = dLogp(PInd) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dW.
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
            Params = ParamListNQSNHTI(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSNHTI(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSSHTI';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Alpha = obj.Alpha; Properties.HDim = obj.HDim;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end