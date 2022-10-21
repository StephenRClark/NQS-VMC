classdef NQSS1 < Modifier
    % NQSS1 - a NQS Modifier variant that uses binary hidden units but
    % visible interaction terms suited for spin-1 systems.
    %   NQS is overarching class, which is itself a subclass of Modifier.

    % Format for NQS Modifier object modified for spin-1:
    % - NQSS1.Nv = number of "visible" spins.
    % - NQSS1.Nh = number of "hidden" spins.
    % - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
    % - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
    % - NQSS1.a = (Nv x 1) vector - visible site bias.
    % - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
    % - NQSS1.A = (Nv x 1) vector - visible site square bias.
    % - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
    % - NQSS1.b = (Nh x 1) vector - hidden site bias.
    % - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
    % - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
    % - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
    % - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
    % - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
    % - NQSS1.Theta = (Nh x 1) vector - effective angles.
    % - NQSS1.VisVec = (Nv x 1) vector - visible occupancies vector.
    % - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.VisVec - vector of new visible occupancies.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nsl x 1) for d/da.
    % - (Nsl x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv x 1) for d/dw.
    % - (Alpha*Nv x 1) for d/dW.
    % ---------------------------------

    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end

    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        Np = 5; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden units.
        Alpha = 1; % Hidden unit density.
        a = 0; % Visible site bias, Nv x 1 vector.
        av = 0; % Visible bias parameters, Nsl x 1 vector.
        A = 0; % Visible site square bias, Nv x 1 vector.
        Av = 0; % Visible square bias parameters, Nsl x 1 vector.
        b = 0; % Hidden site bias, Nh x 1 vector.
        bv = 0; % Hidden bias parameters, Alpha x 1 vector.
        w = 0; % Hidden-visible linear coupling terms, Nh x Nv matrix.
        wm = 0; % Hidden-visible linear coupling parameters, Alpha x Nv matrix.
        W = 0; % Hidden-visible square coupling terms, Nh x Nv matrix.
        Wm = 0; % Hidden-visible square coupling parameters, Alpha x Nv matrix.
        Graph; % Details connectivity of lattice - used to include symmetries.
    end

    properties (Hidden)
        VisVec = 0; % Visible occupancies, Nv x 1 vector.
        NsqVec = 0; % Squared visible occupancies, Nv x 1 vector.
        Theta = 0; % Local effective angle, Nh x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(5,2); % Individual parameter flags for variational purposes.
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullBoseCfg; % Default case assumes bosonic Hilbert.
    end

    methods
        % Constructor for general number hidden NQS:
        function obj = NQSS1(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin == 4 % Assume variational if no VFlag specified.
                obj.VFlag = VFlag;
            elseif nargin < 4
                obj.VFlag = 1;
            end
            if nargin < 3 % Assume instance with zero starting parameters.
                Params.nmag = 0; Params.nphs = 0;
                Params.Alpha = 1;
                Params.a = 0; Params.b = 0; Params.W = 0;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('This NQS variant is incompatible with fermionic systems.');
            elseif strcmp(Hilbert.Type,'Bose')
                obj.FullCfg = @(cfg) FullBoseCfg(cfg) - 1;
                disp('NQSS1 modifier is intended for spins - may have unexpected behaviour for bosons.');
            else
                obj.FullCfg = @FullSpinCfg;
            end
            obj.Nv = Hilbert.N; obj.Graph = Graph;
            obj = RandomInitPsiNQSS1(obj,Params);
        end

        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSS1(obj,dP);
        end

        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.VisVec = Update.VisVec; obj.NsqVec = Update.NsqVec;
        end

        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSS1(obj,Cfg);
        end

        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSS1(obj,Basis);
        end

        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSS1(obj,Params);
        end

        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            Ratio = exp(sum(Diff.val.'.*obj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
            Theta_shift = zeros(obj.Nh,1); % Initialise effective angle shift.
            Vis_shift = zeros(obj.Nv,1); Nsq_shift = zeros(obj.Nv,1);
            % Only loop over the sites where there are differences:
            for i=1:Diff.num
                Vis_shift(Diff.pos(i)) = Diff.val(i);
                Theta_shift = Theta_shift + Diff.val(i)*obj.w(:,Diff.pos(i));
                Nsq_shift(Diff.pos(i)) = 2*Diff.val(i)*obj.VisVec(Diff.pos(i)) + Diff.val(i)^2;
                Theta_shift = Theta_shift + Nsq_shift(Diff.pos(i))*obj.W(:,Diff.pos(i));
            end
            VisP = obj.VisVec + Vis_shift; NsqP = obj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
            Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*obj.A(Diff.pos))); % Compute visible square bias contribution.
            ThetaP = obj.Theta + Theta_shift; % Update the effective angle for the proposed configuration.
            Ratio = Ratio * prod(cosh(ThetaP) ./ cosh(obj.Theta)); % Compute full ratio.
            % Collect new configuration information into Update.
            Update.Theta = ThetaP; Update.VisVec = VisP; Update.NsqVec = NsqP;
        end

        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
            Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
            Nsl = max(SLInds); % Number of sublattices for da.
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            Cfg_sqr = Cfg_vec.^2;
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.
            for s = 1:Nsl
                if sum(obj.OptInds(s,:))~=0
                    dLogp(s) = sum(Cfg_vec(SLInds==s)); % Insert d/da.
                end
                if sum(obj.OptInds(s+Nsl,:))~=0
                    dLogp(s+Nsl) = sum(Cfg_vec(SLInds==s).^2); % Insert d/dA.
                end
            end
            dTheta = tanh(obj.Theta);
            for al=1:obj.Alpha % Derivatives need to be computed by Alpha sector
                bInd = 2*Nsl + al;
                if sum(obj.OptInds(bInd,:)) ~= 0
                    dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
                end
                for n = 1:obj.Nv
                    PIndw = 2*Nsl + obj.Alpha + n + (al-1)*obj.Nv;
                    PIndW = PIndw + obj.Alpha*obj.Nv;
                    % For each layer labelled by a, find the indices of the associated translates.
                    for bd = 1:Ntr
                        TInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                        if VInd ~= 0
                            if sum(obj.OptInds(PIndw,:)) ~= 0
                                dLogp(PIndw) = dLogp(PIndw) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dw.
                            end
                            if sum(obj.OptInds(PIndW,:)) ~= 0
                                dLogp(PIndW) = dLogp(PIndW) + (Cfg_sqr(VInd)*dTheta(TInd)); % Insert d/dW.
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
            Params = ParamListNQSS1(obj);
        end

        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSS1(obj,P);
        end

        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSS1';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end

end