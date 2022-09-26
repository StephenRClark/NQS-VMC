classdef NQSP < NQS
    % NQSP - a NQS Modifier variant that uses number-like hidden units
    % with dimension HDim instead of spin-like hidden units. The hidden and
    % visible units are scaled to the interval [0 1] and polynomial terms
    % are included up to orders VOrder, HOrder.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    %   NQSP draws translational symmetries from the provided Graph object.
    
    % ---------------------------------
    % Format for NQSP Modifier:
    % - NQSP.Nv = number of "visible" spins.
    % - NQSP.Nh = number of "hidden" spins.
    % - NQSP.Np = number of parameters in the ansatz = (Nsl x VOrder) + (Alpha x
    % HOrder) + (Nv x VOrder)(Alpha x HOrder)
    % - NQSP.VDim = dimension of the visible units.
    % - NQSP.HDim = dimension of the hidden units.
    % - NQSP.VOrder = highest power of visible unit interactions. Max value VDim-1.
    % - NQSP.HOrder = highest power of hidden unit interactions. Max value HDim-1.
    % - NQSP.a = (Nv x VOrder) matrix - visible site biases.
    % - NQSP.av = (Nsl x VOrder) matrix - visible bias parameters
    % - NQSP.b = (Nh x HOrder) matrix - hidden site bias.
    % - NQSP.bv = (Alpha x HOrder) matrix - hidden bias parameters.
    % - NQSP.W = (Nh x Nv x HOrder x VOrder) array - hidden-visible coupling terms.
    % - NQSP.Wm = (Alpha x Nv x HOrder x VOrder) array - hidden-visible coupling parameters
    % - NQSP.Alpha = number of unique coupling sets or "hidden unit density".
    % - NQSP.Theta = (Nh x HOrder) matrix - effective angles by hidden order.
    % - NQSP.VisVec = (Nv x 1) vector - visible occupancies.
    % - NQSP.Rescale = flag for visible unit rescaling to [0 1] interval.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - matrix of new effective angles ThetaP.
    % Update.VisVec - vector of new visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nsl x VOrder) x 1 for d/da. Group by Sublattice > Visible order
    % > [sl,vo], [sl, vo+1] ... [sl+1, vo]
    % - (Alpha x HOrder) x 1 for d/db. Group by Alpha > Hidden order
    % > [al, ho], [al, ho+1] ... [al+1, ho]
    % - (Alpha x Nv) x (HOrder x VOrder) for d/dW. Group by Alpha > Position > Hidden order > Visible order
    % > [al,v,ho,vo], [al,v,ho,vo+1] ... [al,v,ho+1,vo] ... [al,v+1,ho,vo] ... [al+1,v,ho,vo] 
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        Np = 2; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden spins.
        Alpha = 1; % Hidden unit density.
        VDim = 1; % Visible unit dimension.
        HDim = 1; % Hidden unit dimension.
        VOrder = 0; % Visible interaction orders.
        HOrder = 0; % Hidden interaction orders.
        a = 0; % Visible site bias terms, Nv x VOrder matrix.
        av = 0; % Visible bias parameters, Nsl x VOrder matrix.
        b = 0; % Hidden spin bias terms, Nh x HOrder matrix.
        bv = 0; % Hidden bias parameters, Alpha x HOrder matrix,
        W = 0; % Hidden-visible coupling terms, Nh x Nv x HOrder x VOrder array.
        Wm = 0; % Coupling parameters, Alpha x Nv x HOrder x VOrder array.
        Graph; % Details connectivity of lattice - used to include symmetries.
        Rescale = 1; % Rescale by default.
    end
    
    properties (Hidden)
        VisVec = 0; % Rescaled visible occupancies n_i/n_max, Nv x 1 vector.
        Theta = 0; % Local effective angle, Nh x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(1,2); % Individual parameter flags for variational purposes.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSP(Hilbert,Graph,Params,VFlag)
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
                Params.VOrder = 1; Params.HOrder = 1;
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
            obj.VDim = round(Hilbert.d); obj.VOrder = round(Params.VOrder);
            obj.HDim = round(Params.HDim); obj.HOrder = round(Params.HOrder);
            if (obj.HDim < 2)
                error('Hidden dimension must be an integer no less than 2.');
            end
            if (obj.VOrder >= obj.VDim)
                disp('Visible interaction order must be less than visible dimension - setting to VDim - 1.');
                obj.VOrder = obj.VDim - 1;
            end
            if (obj.HOrder >= obj.HDim)
                disp('Hidden interaction order must be less than hidden dimension - setting to HDim - 1.');
                obj.HOrder = obj.HDim - 1;
            end
            obj = RandomInitPsiNQSP(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSP(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.VisVec = Update.VisVec;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSP(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSP(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSP(obj,Params);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            VisVecP = obj.VisVec(:);
            VisPow_shift = zeros(Diff.num,obj.VOrder); % Shift in visible occupancy powers for each altered site.
            VisPow_shift(:,1) = Diff.val(:)*((obj.VDim-1)^(-obj.Rescale));
            % Calculate power shifts and visible bias contributions.
            Ratio = exp(sum(obj.a(Diff.pos(:),1).*VisPow_shift(:,1)));
            for vo = 2:obj.VOrder
                VisPow_shift(:,vo) = (VisVecP(Diff.pos(:)+VisPow_shift(:,1))).^vo - VisVecP(Diff.pos(:)).^vo;
                Ratio = Ratio * exp(sum(obj.a(Diff.pos(:),vo).*VisPow_shift(:,vo)));
            end
            Theta_shift = zeros(obj.Nh,obj.HOrder); % Initialise effective angle shift.
            % Theta(h,ho) = sum(W(h,:,ho,vo)*VisVec.^vo)+b(h,ho);
            % Only loop over the sites where there are differences:
            for d=1:Diff.num
                for ho = 1:obj.HOrder
                    Theta_shift(:,ho) = Theta_shift(:,ho) + sum(reshape(...
                        obj.W(:,Diff.pos(d),ho,:),obj.Nh,obj.VOrder).*VisPow_shift(d,:),2);
                end
            end
            ThetaP = obj.Theta + Theta_shift; % Update the effective angles for the proposed configuration.
            % Calculate hidden unit contributions to ratio.
            HidPow = ((0:(obj.HDim-1)).'*((obj.HDim-1)^(-obj.Rescale))) .^ (1:obj.HOrder);
            HidAmp = zeros(obj.Nh,1); HidAmpP = zeros(obj.Nh,1);
            for h = 1:obj.HDim
                HidAmp = HidAmp + exp(obj.Theta*(HidPow(h,:).'));
                HidAmpP = HidAmpP + exp(ThetaP*(HidPow(h,:).'));
            end
            Ratio = Ratio * prod(HidAmpP./HidAmp);
            VisVecP(Diff.pos(:)) = VisVecP(Diff.pos(:))+Diff.val(:)*((obj.VDim-1)^(-obj.Rescale));
            % Add new VisVec and Theta to Update.
            Update.Theta = ThetaP; Update.VisVec = VisVecP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
            Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
            Nsl = max(SLInds); % Number of sublattices for da.
            Cfg_vec = obj.FullCfg(Cfg)*((obj.VDim-1)^-obj.Rescale); % Build the spin configuration vector.
            obj.VisVec = Cfg_vec; % Ensure VisVec is properly assigned.
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            VisPow = Cfg_vec(:) .^ (1:obj.VOrder); % Nv x VOrder
            HidPow = ((0:(obj.HDim-1)).'*((obj.HDim-1)^(-obj.Rescale))) .^ (1:obj.HOrder); % HDim x HOrder
            HidArray = reshape(HidPow.',1,obj.HOrder,obj.HDim);
            dTheta = (sum(exp(obj.Theta.*HidArray).*HidArray,3))./sum(exp(obj.Theta.*HidArray),3); % Nh x HOrder
            for s = 1:Nsl
                for vo = 1:obj.VOrder
                    PInd = vo + (s-1)*obj.VOrder;
                    if sum(obj.OptInds(PInd,:))~=0
                        dLogp(PInd) = sum(VisPow(SLInds==s,vo)); % Insert d/da.
                    end
                end
            end
            for al = 1:obj.Alpha
                for ho = 1:obj.HOrder
                    PInd = Nsl*obj.VOrder + ho + (al-1)*obj.HOrder;
                    if sum(obj.OptInds(PInd,:))~=0
                        dLogp(PInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr,ho));
                    end
                    for v = 1:obj.Nv
                        for vo = 1:obj.VOrder
                            PInd = Nsl*obj.VOrder + obj.Alpha*obj.HOrder + ...
                                vo + obj.VOrder(ho-1 + obj.HOrder*(v-1 + (al-1)*obj.Nv));
                            if sum(obj.OptInds(PInd,:))~=0
                                for bd = 1:numel(BondMap)
                                    HInd = bd + (al-1)*Ntr; VInd = BondMap{bd}(v);
                                    if VInd ~= 0
                                        dLogp(PInd) = dLogp(PInd) + VisPow(VInd,vo)*dTheta(HInd,ho); % Insert d/dW
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp = real(dLogp).*obj.OptInds(:,1) + 1i*imag(dLogp).*obj.OptInds(:,2);
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSP(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSP(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSP';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.Rescale = obj.Rescale;
            Properties.VDim = obj.VDim; Properties.VOrder = obj.VOrder;
            Properties.HDim = obj.HDim; Properties.HOrder = obj.HOrder;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end