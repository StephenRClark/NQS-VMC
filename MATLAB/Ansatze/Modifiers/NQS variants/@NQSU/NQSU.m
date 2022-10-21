classdef NQSU < Modifier
    % NQSU - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden units. U is a parameter-reduced version of OH.
    %   Modifier is the overarching class. NQSU pulls translational
    %   symmetries from the provided Graph object.

    % ---------------------------------
    % Format for NQSU Modifier object:
    % - NQSU.Nv = number of "visible" units.
    % - NQSU.Nh = number of "hidden" units.
    % - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
    % - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
    % - NQSU.VDim = dimensions of the visible units.
    % - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
    % - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
    % - NQSU.b = (Nh x 1) vector - hidden site bias.
    % - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
    % - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
    % - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
    % - NQSU.Theta = (Nh x 1) vector - effective angles.
    % - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP and
    % new unary vector UVecP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nmax*Nsl x 1) for d/da.
    % Arranged [sl, vd], [sl, vd+1], ... , [sl+1, vd], ...
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv*Nmax x 1) for d/dW.
    % Arranged [a, v, vd], [a, v, vd+1], ... ,[a, v+1, vd], ...
    % ---------------------------------

    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end

    properties (SetAccess = protected)
        Np = 1; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden units.
        Alpha = 1; % Hidden unit density.
        VDim = 1; % Visible site dimension.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        av = 0; % Visible bias parameters, Nsl x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        bv = 0; % Hidden bias parameters, Alpha x 1 vector,
        W = 0; % Hidden-visible coupling terms, Nh x Nv matrix.
        Wm = 0; % Coupling parameters, Alpha x Nv matrix.
        Graph % Details connectivity of lattice - used to include symmetries.
    end

    properties (Hidden)
        Theta = 0; % Local effective angle, Nh x 1 vector.
        UVec = 0; % Copy of the configuration converted to a unary vector.
        VList = 0; % Ordered list of values the visible site can take.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(1,1); % Individual parameter flags for variational purposes.
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
    end

    methods
        % Constructor for NQS with randomly initialised parameters:
        function obj = NQSU(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin == 4 % Assume variational if no VFlag specified.
                obj.VFlag = VFlag;
            elseif nargin < 4
                obj.VFlag = 1;
            end
            if nargin < 3 % Assume instance with zero starting parameters.
                Params.nmag = 0; Params.nphs = 0;
                Params.Alpha = 1; Params.VDim = Hilbert.d;
                Params.a = 0; Params.b = 0; Params.W = 0;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('This NQS variant is incompatible with fermionic systems.');
            else
                obj.Nv = Hilbert.N;
                obj.VDim = Hilbert.d;
                if strcmp(Hilbert.Type,'Bose')
                    obj.FullCfg = @FullBoseCfg; obj.VList = (0:(Hilbert.d-1))';
                elseif strcmp(Hilbert.Type,'Spin')
                    obj.FullCfg = @FullSpinCfg;
                    obj.VList = ((0:(obj.VDim-1)).' - Hilbert.S) * (2-mod(obj.VDim,2));
                end
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N; obj.Alpha = Params.Alpha;
            else
                obj.Alpha = ceil(Params.Nh/Hilbert.N); obj.Nh = obj.Alpha*Hilbert.N;
            end
            obj.Graph = Graph;
            obj = RandomInitPsiNQSU(obj,Params);
        end

        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSU(obj,dP);
        end

        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.UVec = Update.UVec;
        end

        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSU(obj,Cfg);
        end

        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSU(obj,Basis);
        end

        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSU(obj,Params);
        end

        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            dV = obj.VList(2) - obj.VList(1); % Needed for half-integer spin cases.
            Nmax = obj.VDim-1;
            % Convert difference in configuration to difference in unary vector.
            dU = zeros(Nmax*obj.Nv,1);
            for d = 1:Diff.num
                SegInds = (1:Nmax) + Nmax*(Diff.pos(d)-1);
                Ind0 = sum((1:Nmax).'.*obj.UVec(SegInds)); IndP = Ind0+(Diff.val(d)/dV);
                dU(SegInds) = -obj.UVec(SegInds);
                if IndP > 0
                    dU(IndP + Nmax*(Diff.pos(d)-1)) = 1;
                end
            end
            ThetaP = obj.Theta + obj.W*dU; UVecP = obj.UVec + dU;
            Ratio = exp(sum(obj.a .* dU)) * prod(cosh(ThetaP)./cosh(obj.Theta));
            Update.UVec = UVecP; Update.Theta = ThetaP;
        end

        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(NQSObj,Cfg)
            Nmax = NQSObj.VDim - 1;
            % Extract information on translational symmetries from Graph.
            GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
            Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
            Nsl = max(SLInds); % Number of sublattices for da.
            Cfg_vec = NQSObj.FullCfg(Cfg); % Build the spin configuration vector.
            UMat = zeros(Nmax,NQSObj.Nv);
            dLogp = zeros(NQSObj.Np,1); % Initialise full vector of derivatives.
            for v = 1:Nmax
                UMat(v,:) = (Cfg_vec.' == NQSObj.VList(v+1));
                for s = 1:Nsl
                    PInd = v + (s-1)*Nmax;
                    if sum(NQSObj.OptInds(PInd,:))~=0
                        % Use sublattice indices to determine which sites contribute
                        dLogp(PInd) = sum(UMat(PInd,GraphObj.SLInds==s));
                    end
                end
            end
            dTheta = tanh(NQSObj.Theta);
            for al = 1:NQSObj.Alpha
                bInd = Nmax*Nsl + al;
                if sum(NQSObj.OptInds(bInd,:)) ~= 0
                    dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr));
                end
                for n = 1:NQSObj.Nv
                    for v = 1:Nmax
                        PInd = Nmax*Nsl + NQSObj.Alpha + Nmax*((n-1)+(al-1)*NQSObj.Nv) + v;
                        if sum(NQSObj.OptInds(PInd,:)) ~= 0
                            for bd = 1:numel(BondMap)
                                HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                                dLogp(PInd) = dLogp(PInd) + UMat(v,VInd)*dTheta(HInd);
                            end
                        end
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp = real(dLogp).*NQSObj.OptInds(:,1) + 1i*imag(dLogp).*NQSObj.OptInds(:,2);
            dLogp(isnan(dLogp)) = 0; dLogp(isinf(dLogp)) = 0;
        end

        % ParamList: outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = ParamListNQSU(obj);
        end

        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSU(obj,P);
        end

        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSU';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.Alpha = obj.Alpha;
            Properties.VDim = obj.VDim; Properties.VList = obj.VList;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end

end