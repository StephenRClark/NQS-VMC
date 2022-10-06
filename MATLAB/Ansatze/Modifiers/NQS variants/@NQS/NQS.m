classdef NQS < Modifier
    % NQS - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden units.
    %   Modifier is the overarching class. NQS itself has subvariants with
    %   symmetries and projections built into them.

    % ---------------------------------
    % Format for NQS Modifier object:
    % - NQS.Nv = number of "visible" units.
    % - NQS.Nh = number of "hidden" units.
    % - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
    % - NQS.Alpha = number of unique coupling sets or "hidden unit density".
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.av = (Nsl x 1) vector - visible bias parameters.
    % - NQS.b = (Nh x 1) vector - hidden site bias.
    % - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
    % - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nsl x 1) for d/da.
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv x 1) for d/dW.
    % Arranged [a, v], [a, v+1] ... [a+1, v] ...
    % ---------------------------------

    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end

    properties (SetAccess = protected)
        Np = 3; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden units.
        Alpha = 1; % Density of hidden units.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        av = 0; % Visible bias parameters, Nsl x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        bv = 0; % Hidden bias parameters, Alpha x 1 vector.
        W = 0; % Hidden-visible coupling terms, Nh x Nv matrix.
        Wm = 0; % Hidden-visible coupling parameters, Alpha x Nv matrix.
        Graph % Details connectivity of lattice - used to include symmetries.
    end

    properties (Hidden)
        Theta = 0; % Local effective angle, Nh x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(3,2); % Individual parameter flags for variational purposes.
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
    end

    methods
        % Constructor for unconstrained NQS with no symmetries:
        function obj = NQS(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
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
            obj = RandomInitPsiNQS(obj,Params);
        end

        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQS(obj,dP);
        end

        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update;
        end

        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQS(obj,Cfg);
        end

        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQS(obj,Basis);
        end

        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQS(obj,Params);
        end

        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            Ratio = exp(sum(Diff.val.'.*obj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
            Theta_shift = zeros(obj.Nh,1); % Initialise effective angle shift.
            % Only loop over the sites where there are differences:
            for i=1:Diff.num
                Theta_shift = Theta_shift + Diff.val(i)*obj.W(:,Diff.pos(i));
            end
            Update = obj.Theta + Theta_shift; % Update the effective angle for the proposed configuration.
            Ratio = Ratio * prod(cosh(Update)./cosh(obj.Theta));
        end

        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            % Extract information on translational symmetries from Graph.
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
            Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
            Nsl = max(SLInds); % Number of sublattices for da.
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.
            for s = 1:Nsl
                if sum(obj.OptInds(s,:)) ~= 0
                    dLogp(s) = sum(Cfg_vec(SLInds==s)); % Insert d/da.
                end
            end
            dTheta = tanh(obj.Theta);
            % Accounting for shift structure of W matrix requires either construction
            % of shifted Theta matrix or shifted Cfg vector - the latter is done here
            for al=1:obj.Alpha % Derivatives need to be computed by Alpha sector
                bInd = Nsl+al;
                if sum(obj.OptInds(Nsl+al,:)) ~= 0
                    dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr)); % Insert d/db.
                end
                for v = 1:obj.Nv
                    PInd = Nsl + obj.Alpha + v + (al-1)*obj.Nv;
                    % For each layer labelled by a, find the indices of the associated translates.
                    if sum(obj.OptInds(PInd,:)) ~= 0
                        for bd = 1:Ntr
                            TInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                            if VInd ~= 0
                                dLogp(PInd) = dLogp(PInd) + (Cfg_vec(VInd)*dTheta(TInd)); % Insert d/dW.
                            end
                        end
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp = real(dLogp).*obj.OptInds(:,1) + 1i*imag(dLogp).*obj.OptInds(:,2);
            dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
        end

        % ParamList: outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = ParamListNQS(obj);
        end

        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQS(obj,P);
        end

        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQS';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end

end