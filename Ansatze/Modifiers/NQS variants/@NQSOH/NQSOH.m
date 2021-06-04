classdef NQSOH < Modifier
    % NQSOH - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden spins. OH refers to one-hot encoding, which splits the
    % visible units into d_v binary units.
    %   Modifier is the overarching class. NQSOH itself has subvariants with
    %   symmetries and projections built into them.
    
    % ---------------------------------
    % Format for NQSOH Modifier object:
    % - NQSOH.Nv = number of "visible" spins.
    % - NQSOH.Nh = number of "hidden" spins.
    % - NQSOH.Np = number of parameters in the ansatz = VDim*Nv + Nh + (VDim*Nv * Nh).
    % - NQSOH.VDim = dimensions of the visible units.
    % - NQSOH.a = (VDim*Nv x 1) vector - visible site bias.
    % - NQSOH.b = (Nh x 1) vector - hidden site bias.
    % - NQSOH.W = (Nh x VDim*Nv) matrix - hidden-visible coupling terms.
    % - NQSOH.Theta = (Nh x 1) vector - effective angles.
    % - NQSOH.VList = (VDim x 1) vector - visible site value list for one-hot encoding.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP and
    % new one-hot vector OHVecP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (VDim*Nv x 1) for d/da.
    % Arranged [v, vd], [v, vd+1], ... , [v+1, vd], ...
    % - (Nh x 1) for d/db.
    % - (Nh*Nv*VDim x 1) for d/dW.
    % Arranged [h, v, vd], [h, v, vd+1], ... ,[h, v+1, vd], ...
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties (SetAccess = protected)
        Np = 3; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden spins.
        VDim = 1; % Visible site dimension.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        W = 0; % Hidden-visible coupling terms, Nh x Nv matrix.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        Theta = 0; % Local effective angle, Nh x 1 vector.
        OHVec = 0; % Copy of the configuration converted to a one-hot vector.
        VList = 0; % Ordered list of values the visible site can take.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(3,1); % Individual parameter flags for variational purposes.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
    end
    
    methods
        % Constructor for unconstrained NQS with no symmetries:
        function obj = NQSOH(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('This NQS variant is incompatible with fermionic systems.');
            else
                obj.Nv = Hilbert.N;
                if Hilbert.d == 2
                    error('This NQS variant is inefficient for d = 2 systems. Use the basic NQS Modifier.');
                else
                    obj.VDim = Hilbert.d;
                    if strcmp(Hilbert.Type,'Bose')
                        obj.FullCfg = @FullBoseCfg; obj.VList = (0:(Hilbert.d-1))';
                    elseif strcmp(Hilbert.Type,'Spin')
                        obj.FullCfg = @FullSpinCfg;
                        obj.VList = ((0:(obj.VDim-1)).' - Hilbert.S) * (2-mod(obj.VDim,2));
                    end
                end
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj.Graph = Graph;
            obj = RandomInitPsiNQSOH(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSOH(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.OHVec = Update.OHVec;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSOH(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSOH(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSOH(obj,Params);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            dV = obj.VList(2) - obj.VList(1); % Needed for half-integer spin cases.
            % Convert difference in configuration to difference in one-hot vector.
            dOH = zeros(obj.VDim*obj.Nv,1);
            for d = 1:Diff.num
                SegInds = (1:obj.VDim) + obj.VDim*(Diff.pos(d)-1); Ind0 = sum((1:obj.VDim).'.*obj.OHVec(SegInds));
                Ind1 = 1+mod(Ind0-1 + (Diff.val(d)/dV),obj.VDim) +  obj.VDim*(Diff.pos(d)-1);
                dOH(SegInds) = -obj.OHVec(SegInds); dOH(Ind1) = 1;
            end
            ThetaP = obj.Theta + obj.W*dOH; OHVecP = obj.OHVec + dOH;
            Ratio = exp(sum(obj.a .* dOH)) * prod(cosh(ThetaP)./cosh(obj.Theta));
            Update.OHVec = OHVecP; Update.Theta = ThetaP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            OHVecT = zeros(obj.VDim,obj.Nv);
            for v = 1:obj.VDim
                OHVecT(v,:) = (Cfg_vec.' == obj.VList(v));
            end
            obj.OHVec = reshape(OHVecT,obj.VDim*obj.Nv,1);            
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.
            dTheta = tanh(obj.Theta);            
            dLogp(1:(obj.VDim*obj.Nv)) = obj.OHVec; % Insert d/da.
            dLogp((1:obj.Nh)+(obj.VDim*obj.Nv)) = dTheta; % Insert d/db.
            dLogp((1:(obj.Nh*obj.Nv*obj.VDim))+obj.Nh+(obj.VDim*obj.Nv)) = ...
                reshape((obj.OHVec * (dTheta.')),obj.Nh*obj.Nv*obj.VDim,1);
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp = dLogp .* obj.OptInds;
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end
        
        % ParamList: outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = ParamListNQSOH(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSOH(obj,P);
        end
        
       % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSOH';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; 
            Properties.VDim = obj.VDim; Properties.VList = obj.VList;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end