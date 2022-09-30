classdef NQSB < Modifier
    % NQSB - a NQS Modifier variant that introduces square visible biases
    % and hidden biases. This variant is tuned for systems with non-binary 
    % visible dimension and hidden dimension.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    %   NQSB draws translational symmetries from the provided Graph object.
    
    % ---------------------------------
    % Format for NQSB Modifier:
    % - NQSB.Nv = number of "visible" units.
    % - NQSB.Nh = number of "hidden" units.
    % - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
    % - NQSB.a = (Nv x 1) vector - visible site bias.
    % - NQSB.av = (Nsl x 1) vector - visible bias parameters.
    % - NQSB.A = (Nv x 1) vector - visible site square bias.
    % - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
    % - NQSB.b = (Nh x 1) vector - hidden site bias.
    % - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
    % - NQSB.B = (Nh x 1) vector- hidden site square bias.
    % - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
    % - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
    % - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
    % - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
    % - NQSB.HDim = dimension of the hidden units.
    % - NQSB.Theta = (Nh x 1) vector - effective angles.
    % - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nsl x 1) for d/da.
    % - (Nsl x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha x 1) for d/dB
    % - (Alpha*Nv x 1) for d/dW.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        Np = 2; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden units.
        Alpha = 1; % Hidden unit density.
        HDim = 1; % Hidden unit dimension.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        av = 0; % Visible bias parameters, Nsl x 1 vector.
        A = 0; % Visible square bias terms, Nv x 1 vector.
        Av = 0; % Visible square bias parameters, Nsl x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        bv = 0; % Hidden bias parameters, Alpha x 1 vector,
        B = 0; % Hidden square bias terms, Nh x 1 vector.
        Bv = 0; % Hidden square bias parameters, Alpha x 1 vector.
        W = 0; % Hidden-visible coupling terms, Nh x Nv matrix.
        Wm = 0; % Coupling parameters, Alpha x Nv matrix.
        Graph; % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        NsqVec = 0; % Squared visible occupancies, Nv x 1 vector.
        Theta = 0; % Local effective angle, Nh x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(1,2); % Individual parameter flags for variational purposes.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullBoseCfg; % Default case assumes bosonic Hilbert.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSB(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin == 4 % Assume variational if no VFlag specified.
                obj.VFlag = VFlag;
            elseif nargin < 4
                obj.VFlag = 1;
            end
            if nargin < 3 % Assume instance with zero starting parameters.
                Params.nmag = 0; Params.nphs = 0;
                Params.Alpha = 1; Params.HDim = Hilbert.d;
                Params.a = 0; Params.b = 0; Params.W = 0;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('This NQS variant is incompatible with fermionic systems.');
            else
                obj.Nv = Hilbert.N;
                if isfield(Params,'HDim') == 0
                    Params.HDim = Hilbert.d; % Assume equal dimension if not specified.
                end
                obj.HDim = Params.HDim;
                if strcmp(Hilbert.Type,'Bose')
                    obj.FullCfg = @FullBoseCfg;
                elseif strcmp(Hilbert.Type,'Spin')
                    obj.FullCfg = @FullSpinCfg;
                end
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N; obj.Alpha = Params.Alpha;
            else
                obj.Alpha = ceil(Params.Nh/Hilbert.N); obj.Nh = obj.Alpha*Hilbert.N;
            end
            obj.Graph = Graph;
            obj = RandomInitPsiNQSB(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSB(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.NsqVec = Update.NsqVec;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSB(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSB(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSB(obj,Params);
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
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
            Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
            Nsl = max(SLInds); % Number of sublattices for da.
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.
            for s = 1:Nsl
                if sum(obj.OptInds(s,:))~=0
                    dLogp(s) = sum(Cfg_vec(SLInds==s)); % Insert d/da.
                end
                if sum(obj.OptInds(s+Nsl,:))~=0
                    dLogp(s+Nsl) = sum(Cfg_vec(SLInds==s).^2); % Insert d/dA.
                end
            end
            dTheta = dT_NHTrace(obj.Theta,obj.B,obj.HDim);
            dB = dB_NHTrace(obj.Theta,obj.B,obj.HDim);
            for al = 1:obj.Alpha
                bInd = 2*Nsl + al; BInd = bInd + obj.Alpha;
                if sum(obj.OptInds(bInd,:)) ~= 0
                    dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
                end
                if sum(obj.OptInds(BInd,:)) ~= 0
                    dLogp(BInd) = sum(dB((1:Ntr)+(al-1)*Ntr)); % Insert d/dB.
                end
                for v = 1:obj.Nv
                    PInd = 2*Nsl + 2*obj.Alpha + (al-1)*obj.Nv + v;
                    if sum(obj.OptInds(PInd,:)) ~= 0
                        for bd = 1:numel(BondMap)
                            HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                            if VInd ~= 0
                                dLogp(PInd) = dLogp(PInd) + Cfg_vec(VInd)*dTheta(HInd);
                            end
                        end
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp = real(dLogp).*obj.OptInds(:,1) + 1i*imag(dLogp).*obj.OptInds(:,2);
            dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSB(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSB(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSB';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.Alpha = obj.Alpha;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end