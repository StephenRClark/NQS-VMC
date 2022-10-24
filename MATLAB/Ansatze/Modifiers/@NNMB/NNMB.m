classdef  NNMB < Modifier
    % NNMB - a Modifier subclass that modifies configuration amplitudes
    % with a nearest neighbour many body (NNMB) term that counts isolated
    % holons and doublons.
    %   Modifier is overarching class. NNMB is not implemented for spin
    %   Hilbert spaces.
    
    % ---------------------------------
    % Format for NNMB Modifier object:
    % - NNMB.GMB = density fluctuation cost.
    % - NNMB.DenFluc = (N x 1) vector of on-site particle density fluctuations.
    % - NNMB.Nmean = mean on-site particle density.
    % ---------------------------------
    % Format for Update is a vector of new density fluctations.
    % ---------------------------------
    % Format for dLogp is a single scalar, dGMB.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
        GMB = 0.05; % Default starting value for G.
    end
    
    properties (SetAccess = protected)
        Np = 1; % Number of parameters.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        OptInds = 1; % One parameter to vary.
        ParamCap = 25; % Parameter cap to avoid runaway errors.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullFermCfg; % Function used by Modifer to interface with Cfg structs.
        Nmean = 0; % Placeholder for mean on-site particle density.
        DenFluc = 1; % Placeholder for vector of density fluctuations.
        Xi = 1; % Placeholder for vector of many body operator values. 
        NNList = []; % (N x z) nearest neighbour list generated from assigned Graph.
    end
    
    methods
        % Constructor for GutzF Modifier object:
        function obj = NNMB(Hilbert,Graph,Params,VFlag)
            % Graph is second argument - required to find nearest
            % neighbours.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Spin') == 1
                error('NNMB Modifier subtype is not compatible with Spin Hilbert objects.')
            end
            if strcmp(Hilbert.Type,'Bose') == 1
                obj.FullCfg = @FullBoseCfg;
            end
            obj.GMB = Params.GMB; obj.Np = 1; obj.OptInds = [VFlag, imag(Params.GMB)~=0];
            % Prepare placeholders for later initialisation in PrepPsi.
            obj.DenFluc = zeros(Hilbert.N,1); obj.Xi = zeros(Hilbert.N,1); obj.Nmean = 0;
            % Generate list of nearest neighbours from Graph.Bonds.
            obj.NNList = ReverseBond(Graph.Bonds); obj.Graph = Graph;
            % Graph defines only nearest neighbours in one direction along
            % each dimension, so require ReverseBond to add others.
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            dP = real(dP)*obj.OptInds(1) + 1i*imag(dP)*obj.OptInds(2);
            obj.GMB = obj.GMB + dP;
            if abs(obj.GMB) > obj.ParamCap
                obj.GMB = sign(obj.GMB) * obj.ParamCap;
            end
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.DenFluc = Update.DenFluc; obj.Xi = Update.MBO;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            Cfg_vec = sum(obj.FullCfg(Cfg),2); N = obj.Graph.N;
            PDen = round(sum(Cfg_vec(:))/obj.Graph.N); obj.Nmean = PDen;
            obj.DenFluc = Cfg_vec - PDen; MBO = zeros(Cfg.N,1);
            for n = 1:N
                MBO(n) = (obj.DenFluc(n) <= -1)*prod(obj.DenFluc(obj.NNList(n,:)) <= 0) ...
                    + (obj.DenFluc(n) >= 1)*prod(obj.DenFluc(obj.NNList(n,:)) >= 0);
            end
            obj.Xi = MBO;
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            DenFlucP = obj.DenFluc; N = numel(DenFlucP); MBO0 = obj.Xi; List = obj.NNList;
            DenFlucP(Diff.pos) = DenFlucP(Diff.pos) + reshape(Diff.val,numel(Diff.val),1);
            % Only need to check sites neighbouring changed sites.
            EList = zeros(N,1); % Flags to check which MBO values need recalculating.
            for d = 1:Diff.num
                EList([Diff.pos(d) List(Diff.pos(d),:)]) = 1;
            end
            MBOP = MBO0;
            for n = 1:N
                if EList(n) == 1
                    MBOP(n) = (DenFlucP(n) <= -1)*prod(DenFlucP(List(n,:)) <= 0) ...
                        + (DenFlucP(n) >= 1)*prod(DenFlucP(List(n,:)) >= 0);
                end
            end
            Ratio = exp(obj.GMB*sum(MBOP - MBO0));
            Update.DenFluc = DenFlucP; Update.MBO = MBOP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,~)
            dLogp = sum(obj.Xi);
        end
        
        % PsiGenerate: Generate full normalised NNMB amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            [Psi] = PsiGenerateNNMB(obj,Basis);
        end
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = obj.GMB;
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj.GMB = P;
            if real(obj.GMB) > obj.ParamCap
                obj.GMB = obj.ParamCap + imag(obj.GMB);
            end
            if real(obj.GMB) < 0
                obj.GMB = 0 + imag(obj.GMB);
            end
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'NNMB';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds; 
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end