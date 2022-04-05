classdef NQSP < NQS
    % NQSP - a NQS Modifier variant that uses number-like hidden units
    % with dimension HDim instead of spin-like hidden units. The hidden and
    % visible units are scaled to the interval [0 1] and polynomial terms
    % are included up to orders VOrder, HOrder.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQS Modifier object with number hidden units:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
    % - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
    % - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
    % interactions. Max value VDim-1.
    % - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
    % interactions. Max value HDim-1.
    % - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
    % HOrder) + (Nv x VOrder)(Nh x HOrder)
    % - NQS.a = (Nv x VOrder) matrix - visible site biases.
    % - NQS.b = (Nh x HOrder) matrix - hidden site bias.
    % - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % - NQS.VisVec = (Nv x 1) vector - visible occupancies.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - matrix of new effective angles ThetaP.
    % Update.VisVec - vector of new visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv x VOrder) x 1 for d/da.
    % - (Nh x HOrder) x 1 for d/db.
    % - (Nh x Nv) x (HOrder x VOrder) for d/dW.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        VDim = 2; % Visible unit dimension.
        VOrder = 1; % Highest visible interaction order.
        HDim = 2; % Hidden unit dimension.
        HOrder = 1; % Highest hidden interaction order.
    end
    
    properties (Hidden)
        VisVec = 0; % Rescaled visible occupancies n_i/n_max, Nv x 1 vector.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSP(Hilbert,Graph,Params,VFlag)
            obj@NQS(Hilbert,Graph,Params,VFlag);
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
            VisPow_shift(:,1) = Diff.val(:)/(obj.VDim-1);
            % Calculate power shifts and visible bias contributions.
            Ratio = exp(sum(a(Diff.pos(:),1).*VisPow_shift(:,1)));
            for vo = 2:obj.VOrder
                VisPow_shift(:,vo) = (VisVecP(Diff.pos(:)+VisPow_shift(:,1))).^vo - VisVecP(Diff.pos(:)).^vo;
                Ratio = Ratio * exp(sum(a(Diff.pos(:),vo).*VisPow_shift(:,vo)));
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
            HidPow = ((0:(obj.HDim-1))/(obj.HDim-1).') .^ (1:obj.HOrder);
            HidAmp = zeros(obj.Nh,1); HidAmpP = zeros(obj.Nh,1);
            for h = 1:obj.HDim
                HidAmp = HidAmp + exp(obj.Theta*(HidPow(h,:).'));
                HidAmpP = HidAmpP + exp(ThetaP*(HidPow(h,:).'));
            end
            Ratio = Ratio * prod(HidAmpP./HidAmp);
            VisVecP(Diff.pos(:)) = VisVecP(Diff.pos(:))+Diff.val(:)/(obj.VDim-1);
            % Add new VisVec and Theta to Update.
            Update.Theta = ThetaP; Update.VisVec = VisVecP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            Cfg_vec = obj.FullCfg(Cfg)/(obj.VDim-1); % Build the spin configuration vector.
            obj.VisVec = Cfg_vec; % Ensure VisVec is properly assigned.
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            VisPow = Cfg_vec(:) .^ (1:obj.VOrder);
            HidPow = ((0:(obj.HDim-1))/(obj.HDim-1).') .^ (1:obj.HOrder);
            HidArray = reshape(HidPow.',1,obj.HOrder,obj.HDim);
            dTheta = (sum(exp(obj.Theta.*HidArray).*HidArray,3))./sum(exp(obj.Theta.*HidArray),3);
            for vo = 1:obj.VOrder
                dLogp(1:Nv+(vo-1)*Nv) = VisPow(:,vo); % Insert d/da.
            end
            for ho = 1:obj.HOrder
                dLogp(1:Nh+(ho-1)*Nh+Nv*obj.VOrder) = dTheta(:,ho);
                for v = 1:Nv
                    for vo = 1:obj.VOrder
                        IndP = Nv*obj.VOrder + Nh*obj.HOrder + Nh*(v-1) + Nh*Nv*(ho-1) + Nh*Nv*obj.HOrder*(vo-1);
                        dLogp(1:Nh+IndP) = dTheta(:,ho)*VisPow(v,vo); % Insert d/dW.
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp = dLogp * (obj.OptInds~=0);
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
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh;
            Properties.VDim = obj.VDim; Properties.VOrder = obj.VOrder;
            Properties.HDim = obj.HDim; Properties.HOrder = obj.HOrder;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end