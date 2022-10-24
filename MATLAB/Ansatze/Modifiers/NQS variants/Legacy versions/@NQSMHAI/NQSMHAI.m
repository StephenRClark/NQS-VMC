classdef NQSMHAI < Modifier
    % NQSMHAI - a NQSMH Modifier variant that incorporates translation
    % invariance using the provided graph. Variant that breaks
    % multiplon-holon symmetry, adding more interactions.
    %   NQS is overarching class, which is itself a subclass of Modifier.
    %   While related to NQSMH enough features differ that NQSMHAI cannot
    %   be instanced as a subclass of NQS.
    
    % ---------------------------------
    % Format for NQS Modifier object with multiplon-holon interactions:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = 2*Alpha*Nv + 2*Alpha + 2.
    % - NQS.A = (Nv x 1) vector - visible site square bias.
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.BH = (Nh x 1) vector - hidden holon bias.
    % - NQS.BM = (Nh x 1) vector - hidden multiplon bias.
    % - NQS.WH = (Nh x Nv) matrix - hidden-visible HH coupling terms.
    % - NQS.WM = (Nh x Nv) matrix - hidden-visible MM coupling terms.
    % - NQS.XH = (Nh x Nv) matrix - hidden-visible HM coupling terms.
    % - NQS.XM = (Nh x Nv) matrix - hidden-visible MH coupling terms.
    % - NQS.ThetaH = (Nh x 1) vector - effective angles for hidden holons.
    % - NQS.ThetaM = (Nh x 1) vector - effective angles for hidden multiplons.
    % - NQS.Hv = (Nv x 1) vector - vector of visible holons.
    % - NQS.Mv = (Nv x 1) vector - vector of visible multiplons.
    % - NQS.NsqVec = (Nv x 1) vector - squared visible occupancies.
    % - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
    % Properties added with translation invariance:
    % - NQS.Ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.ati = (1 x 1) scalar - reduced parameter set for TI.
    % - NQS.BHti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.BMti = (Alpha x 1) vector - reduced parameter set for TI.
    % - NQS.WHv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.WMv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.XHv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.XMv = (Alpha x Nv) matrix - reduced parameter set for TI.
    % - NQS.Alpha = hidden unit density / number of unique sets of W couplings.
    % ---------------------------------
    % Format for Update is a struct with five fields:
    % Update.ThetaH - vector of new effective angles ThetaHP.
    % Update.ThetaM - vector of new effective angles ThetaMP.
    % Update.NsqVec - vector of new squared visible occupancies.
    % Update.Hv - vector of new holon operator values HvP.
    % Update.Mv - vector of new multiplon operator values MvP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (1 x 1) for d/da.
    % - (1 x 1) for d/dA.
    % - (Alpha x 1) for d/db.
    % - (Alpha x 1) for d/dB.
    % - (Alpha*Nv x 1) for d/dWH.
    % - (Alpha*Nv x 1) for d/dWM.
    % - (Alpha*Nv x 1) for d/dXH.
    % - (Alpha*Nv x 1) for d/dXM.
    % ---------------------------------
    
    properties (SetAccess = protected) % Default to one visible, one hidden plus state with no input.
        Np = 8; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden spins.
        HDim = 2; % Hidden unit dimension.
        Graph % Details connectivity of lattice - used to include symmetries.
        ati = 0; % Visible site bias, 1 x 1 scalar.
        a = 0; % Visible site bias vector, Nv x 1 vector.
        Ati = 0; % Visible site square bias, 1 x 1 scalar.
        A = 0; % Visible site square bias vector, Nv x 1 vector.
        BHti = 0; % Hidden site holon bias, Alpha x 1 vector.
        BH = 0; % Hidden site holon biases, Nh x 1 vector.
        BMti = 0; % Hidden site multiplon bias, Alpha x 1 vector.
        BM = 0; % Hidden site multiplon biases, Nh x 1 vector.
        WHv = 0; % Hidden-visible HH couplings, Alpha x Nv matrix.
        WH = 0; % Hidden-visible HH coupling matrix, Nh x Nv.
        WMv = 0; % Hidden-visible MM couplings, Alpha x Nv matrix.
        WM = 0; % Hidden-visible MM coupling matrix, Nh x Nv.
        XHv = 0; % Hidden-visible HM couplings, Alpha x Nv matrix.
        XH = 0; % Hidden-visible HM coupling matrix, Nh x Nv.
        XMv = 0; % Hidden-visible MH couplings, Alpha x Nv matrix.
        XM = 0; % Hidden-visible MH coupling matrix, Nh x Nv.
        Alpha = 1; % Hidden unit density / number of unique sets of couplings.
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
        % Constructor for multiplon-holon NQS with translation invariance:
        function obj = NQSMHAI(Hilbert,Graph,Params,VFlag)
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Bose') == 0
                error('This NQS variant is only implemented for bosonic systems.');
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj.Nv = Hilbert.N; obj.HDim = Hilbert.d; obj.VFlag = VFlag;
            obj.Graph = Graph; obj.Hilbert = Hilbert;
            obj = RandomInitPsiNQSMHAI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateNQSMHAI(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdateNQSMH(obj,Update);
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiNQSMHAI(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateNQSMHAI(obj,Basis);
        end
        
        % AddHidden: Generate additional hidden units and associated
        % parameters.
        function [obj] = AddHidden(obj,Params)
            obj = AddHiddenNQSMHAI(obj,Params);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = PsiRatioNQSMHAI(obj,Diff);
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = LogDerivNQSMHAI(obj,Cfg);
        end
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            [Params] = ParamListNQSMHAI(obj);
        end
    end
    
end