classdef NQSNH < Modifier
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
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        A = 0; % Visible site square bias.
        B = 0; % Hidden site square bias.
        HDim = 2; % Hidden unit dimensionality, set to match Hilbert.
        a = 0;
        b = 0;
        W = 0;
        Nv = 0;
        Nh = 0;
        Np = 0;
        Graph
    end
    
    properties (Hidden)
        NsqVec = 0; % Squared visible occupancies, Nv x 1 vector.
        ParamCap = 5;
        Theta = 0;
        OptInds = 0;
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSNH(Hilbert,Graph,Params,VFlag)
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('Base NQS compatibility with fermionic systems has been retired. Use NQSF modifier instead.');
            else
                obj.Nv = Hilbert.N;
                if strcmp(Hilbert.Type,'Bose')
                    obj.FullCfg = @FullBoseCfg;
                elseif strcmp(Hilbert.Type,'Spin')
                    obj.FullCfg = @FullSpinCfg;
                end
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj.Graph = Graph;
            obj.HDim = Params.HDim; % Set NQS hidden dimension in Params.
            if (obj.HDim < 2) || (floor(obj.HDim) ~= obj.HDim)
                error('Hidden dimension must be an integer no less than 2.');
            end
            obj = RandomInitPsiNQSNH(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSNH(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.NsqVec = Update.NsqVec;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSNH(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSNH(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSNH(obj,Params);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            Ratio = exp(sum(Diff.val.'.*obj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
            Theta_shift = zeros(obj.Nh,1); % Initialise effective angle shift.
            Nsq_shift = zeros(obj.Nv,1);
            % Only loop over the sites where there are differences:
            for i=1:Diff.num
                Theta_shift = Theta_shift + Diff.val(i)*obj.W(:,Diff.pos(i));
                Nsq_shift(Diff.pos(i)) = 2*sqrt(obj.NsqVec(Diff.pos(i)))*Diff.val(i) + (Diff.val(i)^2);
            end
            NsqP = obj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
            Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*obj.A(Diff.pos))); % Compute visible square bias contribution.
            ThetaP = obj.Theta + Theta_shift; % Update the effective angle for the proposed configuration.
            Ratio = Ratio * prod(NHTrace(ThetaP,obj.B,obj.HDim) ./ ...
                NHTrace(obj.Theta,obj.B,obj.HDim)); % Compute full ratio.
            % Collect new configuration information into Update.
            Update.Theta = ThetaP; Update.NsqVec = NsqP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.            
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            dLogp(1:obj.Nv) = Cfg_vec; % Insert d/da.
            dLogp((1:obj.Nv)+obj.Nv) = Cfg_vec.^2; % Insert d/dA.
            dLogp((1:(obj.Nv*obj.Nh))+2*(obj.Nv+obj.Nh)) = ...
                reshape((dTSBTrace(obj.Theta,obj.HDim,obj.B)*Cfg_vec.'),obj.Nh*obj.Nv,1); % Insert d/dW.
            dLogp((1:obj.Nh)+2*obj.Nv) = dTSBTrace(obj.Theta,obj.HDim,obj.B); % Insert d/db.
            dLogp((1:obj.Nh)+2*obj.Nv+obj.Nh) = dBSBTrace(obj.Theta,obj.HDim,obj.B); % Insert d/dB.
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;            
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSNH(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSNH(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSNH';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.HDim = obj.HDim;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end