classdef NQS < Modifier
    % NQS - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden spins. Fermionic NQS will split configurations into a
    % bipartite spin system architecture.
    %   Modifier is the overarching class. NQS itself has subvariants with
    %   symmetries and projections built into them.
    
    % ---------------------------------
    % Format for NQS Modifier object:
    % - NQS.Nv = number of "visible" spins.
    % - NQS.Nh = number of "hidden" spins.
    % - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
    % - NQS.a = (Nv x 1) vector - visible site bias.
    % - NQS.b = (Nh x 1) vector - hidden site bias.
    % - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
    % - NQS.Theta = (Nh x 1) vector - effective angles.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv x 1) for d/da.
    % - (Nh x 1) for d/db.
    % - (Nh*Nv x 1) for d/dW.
    % ---------------------------------
    
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties (SetAccess = protected)
        Np = 3; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden spins.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        W = 0; % Hidden-visible coupling terms, Nh x Nv matrix.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        Theta = 0; % Local effective angle, Nh x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(3,1); % Individual parameter flags for variational purposes.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
        SFlag = 0; % Spin symmetry flag for fermionic case.
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
                obj.Nv = 2*Hilbert.N; obj.FFlag = Hilbert.N; obj.FullCfg = @BiPtFermCfg;
                if (isempty(Hilbert.Sector) == 0) && (Hilbert.Sector(1) == Hilbert.Sector(2))
                    obj.SFlag = 1;
                end
                % For NQS, fermionic extension simply doubles the size of
                % the effective lattice, though if N_up = N_dn, spin
                % symmetry is imposed as well.
            else
                obj.Nv = Hilbert.N; obj.FFlag = 0;
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
            if obj.FFlag ~= 0 % Fermionic Diff conversion required.
                [Diff] = Ferm2BiPtDiff(Diff,obj.FFlag);
                % Ferm2BiPtDiff converts differences for pairs of
                % occupation numbers to differences in occupation on a
                % double sized lattice.
            end
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
            Cfg_vec = obj.FullCfg(Cfg); % Build the spin configuration vector.            
            dLogp = zeros(obj.Np,1); % Initialise full vector of derivatives.            
            dLogp(1:obj.Nv) = Cfg_vec; % Insert d/da.
            dLogp((obj.Nv+obj.Nh+1):(obj.Nv+obj.Nh+obj.Nv*obj.Nh)) = ...
                reshape((tanh(obj.Theta)*Cfg_vec.'),obj.Nh*obj.Nv,1); % Insert d/dW.
            dLogp((obj.Nv+1):(obj.Nv+obj.Nh)) = tanh(obj.Theta); % Insert d/db.
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;            
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