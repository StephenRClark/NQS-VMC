classdef NQSS1 < NQS
    % NQSS1 - a NQS Modifier variant that uses binary hidden units but
    % visible interaction terms suited for spin-1 systems
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % Format for NQS Modifier object modified for spin-1:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.A = (Nv x 1) vector - visible site square bias.
    % - NQS.b = (Nh x 1) vector - hidden site bias.
    % - NQS.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv x 1) for d/da.
    % - (Nv x 1) for d/dA.
    % - (Nh x 1) for d/db.
    % - (Nh*Nv x 1) for d/dw.
    % - (Nh*Nv x 1) for d/dW.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        A = 0; % Visible site square bias.
        B = 0; % Hidden site square bias.
        w = 0; % Hidden-visible square coupling.
    end
    
    properties (Hidden)
        VisVec = 0; % Visible occupancies, Nv x 1 vector.
        NsqVec = 0; % Squared visible occupancies, Nv x 1 vector.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSS1(Hilbert,Graph,Params,VFlag)
            obj@NQS(Hilbert,Graph,Params,VFlag);
            if strcmp(Hilbert.Type,'Bose')
                obj.FullCfg = @(cfg) FullBoseCfg(cfg) - 1;
                disp('NQSS1 modifier is intended for spins - may have unexpected behaviour for bosons.');
            end
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
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            Cfg_sqr = Cfg_vec.^2;            
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            dTheta = tanh(obj.Theta);            
            dLogp(1:obj.Nv) = Cfg_vec; % Insert d/da.
            dLogp((1:obj.Nv)+obj.Nv) = Cfg_sqr; % Insert d/dA.
            dLogp((1:obj.Nh)+2*obj.Nv) = dTheta; % Insert d/db.
            dLogp((1:(obj.Nv*obj.Nh))+2*obj.Nv+obj.Nh) = ...
                reshape((dTheta*Cfg_vec.'),obj.Nh*obj.Nv,1); % Insert d/dw.
            dLogp((1:(obj.Nv*obj.Nh))+2*obj.Nv+obj.Nh+obj.Nv*obj.Nh) = ...
                reshape((dTheta*Cfg_sqr.'),obj.Nh*obj.Nv,1); % Insert d/dW.            
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;            
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSS1(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSS1(obj,P);
        end
        
        % ChangeFullCfg: change out the configuration reading function.
        function [obj] = ChangeFullCfg(obj,F_new)
            disp(['Configuration reading function changed from ' func2str(obj.FullCfg) ...
                ' to ' func2str(F_new)]);
            obj.FullCfg = F_new;            
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