classdef NQSSX < NQS
    % NQSSX - a NQS Modifier variant that uses a square-square interaction
    % term X, with multinomial visible and hidden units.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for NQS Modifier object with square-square interaction:
    % - NQSSX.Nv = number of "visible" spins.
    % - NQSSX.Nh = number of "hidden" spins.
    % - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
    % - NQSSX.a = (Nv x 1) vector - visible site bias.
    % - NQSSX.A = (Nv x 1) vector - visible site square bias.
    % - NQSSX.b = (Nh x 1) vector - hidden site bias.
    % - NQSSX.B = (Nh x 1) vector - hidden site square bias.
    % - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
    % - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
    % - NQSSX.HDim = dimension of the hidden units.
    % - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
    % - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
    % - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
    % - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
    % - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % ---------------------------------
    % Format for Update is a struct with two fields:
    % Update.Theta - vector of new effective angles ThetaP.
    % Update.VisVec - vector of new visible occupancies.
    % Update.ThetaSq - vector of new effective angles ThetaSqP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv x 1) for d/da.
    % - (Nv x 1) for d/dA.
    % - (Nh x 1) for d/db.
    % - (Nh x 1) for d/dB.
    % - (Nh*Nv x 1) for d/dW.
    % - (Nh*Nv x 1) for d/dX.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        A = 0; % Visible site square bias.
        B = 0; % Hidden site square bias.
        X = 0; % Hidden-visible square coupling.
        HVal = 0; % Hidden unit values, 1 x HDim vector.
        HDim = 1; % Hidden unit dimension.
    end
    
    properties (Hidden)
        VisVec = 0; % Visible occupancies, Nv x 1 vector.
        ThetaSq = 0; % Effective square-hidden angles, Nh x 1 vector.
        NsqVec = 0; % Squared visible occupancies, Nv x 1 vector.
    end
    
    methods
        % Constructor for general number hidden NQS:
        function obj = NQSSX(Hilbert,Graph,Params,VFlag)
            obj@NQS(Hilbert,Graph,Params,VFlag);
            if strcmp(Hilbert.Type,'Bose')
                obj.FullCfg = @(cfg) FullBoseCfg(cfg) - 1;
            end
            obj = RandomInitPsiNQSSX(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSSX(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.VisVec = Update.VisVec; 
            obj.ThetaSq = Update.ThetaSq; obj.NsqVec = Update.NsqVec;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSSX(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSSX(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSSX(obj,Params);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            Ratio = exp(sum(Diff.val.'.*obj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
            Theta_shift = zeros(obj.Nh,1); ThetaSq_shift = zeros(obj.Nh,1); % Initialise effective angle shift.
            Vis_shift = zeros(obj.Nv,1); Nsq_shift = zeros(obj.Nv,1);
            % Only loop over the sites where there are differences:
            for i=1:Diff.num
                Vis_shift(Diff.pos(i)) = Diff.val(i);
                Theta_shift = Theta_shift + Diff.val(i)*obj.W(:,Diff.pos(i));
                Nsq_shift(Diff.pos(i)) = 2*Diff.val(i)*obj.VisVec(Diff.pos(i)) + Diff.val(i)^2;
                ThetaSq_shift = ThetaSq_shift + Nsq_shift(Diff.pos(i))*obj.X(:,Diff.pos(i));
            end
            VisP = obj.VisVec + Vis_shift; NsqP = obj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
            Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*obj.A(Diff.pos))); % Compute visible square bias contribution.
            ThetaP = obj.Theta + Theta_shift; ThetaSqP = obj.ThetaSq + ThetaSq_shift;% Update the effective angle for the proposed configuration.
            Ratio = Ratio * prod(SqTrace(ThetaP,ThetaSqP,obj.HVal) ./ ...
                SqTrace(obj.Theta,obj.ThetaSq,obj.HVal)); % Compute full ratio.
            % Collect new configuration information into Update.
            Update.Theta = ThetaP; Update.ThetaSq = ThetaSqP; Update.VisVec = VisP; Update.NsqVec = NsqP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            Cfg_sqr = Cfg_vec.^2;            
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            dTheta = dL_SqTrace(obj.Theta,obj.ThetaSq,obj.HVal);
            dThetaSq = dS_SqTrace(obj.Theta,obj.ThetaSq,obj.HVal);
            dLogp(1:obj.Nv) = Cfg_vec; % Insert d/da.
            dLogp((1:obj.Nv)+obj.Nv) = Cfg_sqr; % Insert d/dA.
            dLogp((1:obj.Nh)+2*obj.Nv) = dTheta; % Insert d/db.
            dLogp((1:obj.Nh)+2*obj.Nv+obj.Nh) = dThetaSq; % Insert d/db.
            dLogp((1:(obj.Nv*obj.Nh))+2*obj.Nv+2*obj.Nh) = ...
                reshape((dTheta*Cfg_vec.'),obj.Nh*obj.Nv,1); % Insert d/dW.
            dLogp((1:(obj.Nv*obj.Nh))+2*obj.Nv+2*obj.Nh+obj.Nv*obj.Nh) = ...
                reshape((dThetaSq*Cfg_sqr.'),obj.Nh*obj.Nv,1); % Insert d/dX.            
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;            
        end
        
        % ParamList: outputs a Np x 1 vector of parameters.
        function [Params] = ParamList(obj)
            Params = ParamListNQSSX(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSSX(obj,P);
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
            Properties.Type = 'NQSSX';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end