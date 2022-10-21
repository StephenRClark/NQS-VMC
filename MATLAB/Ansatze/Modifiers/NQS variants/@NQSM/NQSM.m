classdef NQSM < Modifier
    % NQSM - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden units. M is a parameter-reduced version of U.
    %   Modifier is the overarching class. NQSM pulls translational
    %   symmetries from the provided Graph object.

    % ---------------------------------
    % Format for NQSM Modifier object:
    % - NQSM.Nv = number of "visible" units.
    % - NQSM.Nh = number of "hidden" units.
    % - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
    % - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
    % - NQSM.a = (3*Nv x 1). vector - visible site bias.
    % - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
    % - NQSM.b = (Nh x 1) vector - hidden site bias.
    % - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
    % - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
    % - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
    % - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
    % - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
    % - NQSM.Theta = (Nh x 1) vector - effective angles.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP, new
    % visible occupancies VisVec and new species matrix OMatP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (3*Nsl x 1) for d/da.
    % Arranged [sl, H], [sl+1, H] ... [sl, D] ... [sl, M] ...
    % - (Alpha x 1) for d/db.
    % - (Alpha*Nv x 1) for d/dW.
    % Arranged [a, v], [a, v+1] ... [a+1, v] ...
    % - (Alpha*Nv x 1) for d/dX.
    % Arranged [a, v], [a, v+1] ... [a+1, v] ...
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
        a = 0; % Visible site bias terms, 3*Nv x 1 vector.
        av = 0; % Visible bias parameters, 3*Nsl x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        bv = 0; % Hidden bias parameters, Alpha x 1 vector,
        W = 0; % Holon coupling terms, Nh x Nv matrix.
        Wm = 0; % Holon coupling parameters, Alpha x Nv matrix.
        X = 0; % Doublon coupling terms, Nh x Nv matrix.
        Xm = 0; % Doublon coupling parameters, Alpha x Nv matrix.
        Graph % Details connectivity of lattice - used to include symmetries.
    end

    properties (Hidden)
        Theta = 0; % Local effective angle, Nh x 1 vector.
        VisVec = 0; % Visible occupancies, Nv x 1 vector.
        OMat = 0; % Matrix of holon, doublon and multiplon values, Nv x 3 matrix.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(1,1); % Individual parameter flags for variational purposes.
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullBoseCfg; % Default case assumes bosonic Hilbert.
    end

    methods
        % Constructor for NQS with randomly initialised parameters:
        function obj = NQSM(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin == 4 % Assume variational if no VFlag specified.
                obj.VFlag = VFlag;
            elseif nargin < 4
                obj.VFlag = 1;
            end
            if nargin < 3 % Assume instance with zero starting parameters.
                Params.nmag = 0; Params.nphs = 0;
                Params.Alpha = 1;
                Params.a = 0; Params.b = 0; Params.W = 0; Params.X = 0;
            end
            if strcmp(Hilbert.Type,'Bose')==0
                error('This NQS variant is only intended for bosonic systems.');
            else
                obj.Nv = Hilbert.N;
                obj.VDim = Hilbert.d;
                obj.FullCfg = @FullBoseCfg;
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N; obj.Alpha = Params.Alpha;
            else
                obj.Alpha = ceil(Params.Nh/Hilbert.N); obj.Nh = obj.Alpha*Hilbert.N;
            end
            obj.Graph = Graph;
            obj = RandomInitPsiNQSM(obj,Params);
        end

        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSM(obj,dP);
        end

        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.OMat = Update.OMat; obj.VisVec = Update.VisVec;
        end

        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSM(obj,Cfg);
        end

        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSM(obj,Basis);
        end

        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSM(obj,Params);
        end

        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            VisVecP = obj.VisVec; OMatP = obj.OMat; ThetaP = obj.Theta;
            for d = 1:Diff.num
                VisVecP(Diff.pos(d)) = VisVecP(Diff.pos(d)) + Diff.val(d);
                switch VisVecP(Diff.pos(d))
                    case 0
                        OMatP(Diff.pos(d),:) = [1 0 0];
                    case 1
                        OMatP(Diff.pos(d),:) = [0 0 0];
                    case 2
                        OMatP(Diff.pos(d),:) = [0 1 0];
                    otherwise
                        OMatP(Diff.pos(d),:) = [0 0 1];
                end
            end
            dHDM = reshape((OMatP-obj.OMat),3*obj.Nv,1); dHD = dHDM(1:(2*obj.Nv));
            ThetaP = ThetaP + [obj.W, obj.X] * dHD;
            Ratio = exp(sum(obj.a .* dHDM)) * prod(cosh(ThetaP)./cosh(obj.Theta));
            Update.OMat = OMatP; Update.VisVec = VisVecP; Update.Theta = ThetaP;
        end

        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            % Extract information on translational symmetries from Graph.
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
            Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
            Nsl = max(SLInds); % Number of sublattices for da.
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            OMatP = [(Cfg_vec==0), (Cfg_vec==2), (Cfg_vec>2)];
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.
            for s = 1:Nsl
                for v = 1:3
                    PInd = v + (s-1)*3;
                    if sum(obj.OptInds(PInd,:))~=0
                        dLogp((1:3) + (s-1)*3) = sum(OMatP(SLInds==s,:),1).'; % Insert d/da
                    end
                end
            end
            dTheta = tanh(obj.Theta);
            for al = 1:obj.Alpha
                bInd = 3*Nsl + al;
                if sum(obj.OptInds(bInd,:)) ~= 0
                    dLogp(bInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr));
                end
                for n = 1:obj.Nv
                    PIndW = 3*Nsl + obj.Alpha + (al-1)*obj.Nv + n;
                    PIndX = PIndW + obj.Alpha*obj.Nv;
                    % Big if / elseif to avoid extraneous if checking or loops
                    if (sum(obj.OptInds(PIndW,:))~=0) && (sum(obj.OptInds(PIndX))~=0)
                        for bd = 1:numel(BondMap)
                            HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                            dLogp(PIndW) = dLogp(PIndW) + OMatP(VInd,1)*dTheta(HInd);
                            dLogp(PIndX) = dLogp(PIndX) + OMatP(VInd,2)*dTheta(HInd);
                        end
                    elseif (sum(obj.OptInds(PIndW,:))~=0) && (sum(obj.OptInds(PIndX))==0)
                        for bd = 1:numel(BondMap)
                            HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                            dLogp(PIndW) = dLogp(PIndW) + OMatP(VInd,1)*dTheta(HInd);
                        end
                    elseif (sum(obj.OptInds(PIndW,:))==0) && (sum(obj.OptInds(PIndX))~=0)
                        for bd = 1:numel(BondMap)
                            HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(n);
                            dLogp(PIndX) = dLogp(PIndX) + OMatP(VInd,2)*dTheta(HInd);
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
            Params = ParamListNQSM(obj);
        end

        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSM(obj,P);
        end

        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSM';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.Alpha = obj.Alpha;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end

end