classdef NQSMH < Modifier
    % NQSMH - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible and
    % hidden neurons, intended for systems where the visible on-site
    % dimension is greater than 2. The interaction terms are phrased in
    % terms of holons and multiplons.
    %   Modifier is the overarching class. NQSMH itself has subvariants with
    %   symmetries and projections built into them. While related to NQS,
    %   enough features differ that NQSMH cannot be instanced as a
    %   subclass.
    
    % ---------------------------------
    % Format for NQS Modifier object with multiplon-holon interactions:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = 2*Nv*Nv + 2*Nh + 2*Nv.
    % - NQS.A = (Nv x 1) vector - visible site square bias.
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.BH = (Nh x 1) vector - hidden holon bias.
    % - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible MM/HH coupling terms.
    % - NQS.X = (Nh x Nv) matrix - hidden-visible MH/HM coupling terms.
    % - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
    % - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
    % - NQS.Hv = (Nv x 1) vector - vector of visible holons.
    % - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
    % - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
    % ---------------------------------
    % Format for Update is a struct with five fields:
    % Update.ThetaH - vector of new effective angles ThetaHP.
    % Update.ThetaM - vector of new effective angles ThetaMP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % Update.Hv - vector of new holon operator values HvP.
    % Update.Mv - vector of new multiplon operator values MvP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv x 1) for d/da.
    % - (Nv x 1) for d/dA.
    % - (Nh x 1) for d/dBH.
    % - (Nh x 1) for d/dBM.
    % - (Nh*Nv x 1) for d/dW.
    % - (Nh*Nv x 1) for d/dX.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties (SetAccess = protected)
        Np = 6; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden spins.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        A = 0; % Visible site square bias terms, Nv x 1 vector.
        BH = 0; % Hidden holon bias terms, Nh x 1 vector.
        BM = 0; % Hidden multiplon bias terms, Nh x 1 vector.
        W = 0; % Hidden-visible MM/HH coupling terms, Nh x Nv matrix.
        X = 0; % Hidden-visible MH/HM coupling terms, Nh x Nv matrix.
        HDim = 2; % Hidden unit dimension.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        ThetaH = 0; % Local effective hidden holon angle, Nh x 1 vector.
        ThetaM = 0; % Local effective hidden multiplon angle, Nh x 1 vector.
        Hv = 0; % Visible holon operator values, Nv x 1 vector.
        Mv = 0; % Visible multiplon operator values, Nv x 1 vector.
        NsqVec = 0; % Visible squared values, Nv x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(6,1); % Individual parameter flags for variational purposes.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullBoseCfg; % Default case assumes bosonic Hilbert.
    end
    
    methods
        % Constructor for multiplon-holon NQS with no symmetries:
        function obj = NQSMH(Hilbert,Graph,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('Fermionic systems are not currently supported for this Modifier.');
            else
                obj.Nv = Hilbert.N;
                obj.FullCfg = Hilbert.FullCfgFunc;
                obj.HDim = Hilbert.d;
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj.Graph = Graph;
            obj = RandomInitPsiNQSMH(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSMH(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.ThetaH = Update.ThetaH; obj.ThetaM = Update.ThetaM;
            obj.NsqVec = Update.NsqVec; obj.Hv = Update.Hv; obj.Mv = Update.Mv;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSMH(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSMH(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSMH(obj,Params);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            Ratio = exp(sum(Diff.val.'.*obj.a(Diff.pos))); % Initialise the ratio with the a-vector contribution.
            ThetaH_shift = zeros(obj.Nh,1); ThetaM_shift = zeros(obj.Nh,1); % Initialise effective angle shifts.
            Nsq_shift = zeros(obj.Nv,1); dHv = zeros(obj.Nv,1); dMv = zeros(obj.Nv,1);
            % Only loop over the sites where there are differences:
            for i=1:Diff.num
                dHv(Diff.pos(i)) = ((obj.Mv(Diff.pos(i)) + Diff.val(i)) < 0) - (Diff.val(i) > 0)*obj.Hv(Diff.pos(i));
                dMv(Diff.pos(i)) = Diff.val(i) + dHv(Diff.pos(i)); % Holon / multiplon operator differences.
                Nsq_shift(Diff.pos(i)) = 2*sqrt(obj.NsqVec(Diff.pos(i)))*Diff.val(i) + (Diff.val(i)^2);
                ThetaH_shift = ThetaH_shift + dHv(Diff.pos(i))*obj.W(:,Diff.pos(i)) + ...
                    dMv(Diff.pos(i))*obj.X(:,Diff.pos(i)); % Change in hidden holon effective angles.
                ThetaM_shift = ThetaM_shift + dMv(Diff.pos(i))*obj.W(:,Diff.pos(i)) + ...
                    dHv(Diff.pos(i))*obj.X(:,Diff.pos(i)); % Change in hidden multiplon effective angles.
            end
            NsqP = obj.NsqVec + Nsq_shift;% Update the squared occupancy vector for the proposed configuration.
            HvP = obj.Hv + dHv; MvP = obj.Mv + dMv; % Update the holon / multiplon operator vectors for the proposed configuration.
            Ratio = Ratio * exp(sum(Nsq_shift(Diff.pos).*obj.A(Diff.pos))); % Compute visible square bias contribution.
            ThetaHP = obj.ThetaH + ThetaH_shift; ThetaMP = obj.ThetaM + ThetaM_shift;
            % Update the effective angles for the proposed configuration.
            Ratio = Ratio * prod(MHTrace(ThetaHP,ThetaMP,obj.HDim) ./ ...
                MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim)); % Compute full ratio.
            % Collect new configuration information into Update.
            Update.ThetaH = ThetaHP; Update.ThetaM = ThetaMP;
            Update.NsqVec = NsqP; Update.Hv = HvP; Update.Mv = MvP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)          
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.
            obj.Hv = (Cfg_vec == 0); obj.Mv = (Cfg_vec - 1) .* (Cfg_vec > 0); % Build the holon / multiplon vectors.            
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            dLogp(1:obj.Nv) = Cfg_vec; % Insert d/da.
            dLogp((1:obj.Nv)+obj.Nv) = Cfg_vec.^2; % Insert d/dA.            
            dLogp((1:obj.Nh)+2*obj.Nv) = dTH_MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim); % Insert d/dBH.
            dLogp((1:obj.Nh)+2*obj.Nv+obj.Nh) = dTM_MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim); % Insert d/dBM.            
            dLogp((1:(obj.Nv*obj.Nh))+2*(obj.Nv+obj.Nh)) = ...
                reshape( ((dTH_MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim)*obj.Hv.') + ...
                (dTM_MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim)*obj.Mv.')), obj.Nh*obj.Nv, 1); % Insert d/dW.            
            dLogp((1:(obj.Nv*obj.Nh))+2*(obj.Nv+obj.Nh)+(obj.Nv*obj.Nh)) = ...
                reshape( ((dTH_MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim)*obj.Mv.') + ...
                (dTM_MHTrace(obj.ThetaH,obj.ThetaM,obj.HDim)*obj.Xv.')), obj.Nh*obj.Nv, 1);% Insert d/dX.            
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end
        
        % ParamList: outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = ParamListNQSMH(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadNQSMH(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NQSMH';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds; 
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.HDim = obj.HDim;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end