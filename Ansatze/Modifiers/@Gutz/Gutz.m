classdef  Gutz < Modifier
    % Gutz - a Modifier subclass that modifies configuration amplitudes
    % according to the local particle density at each site. Incompatible
    % with Spin and Boson Hilbert spaces.
    %   Modifier is overarching class. Gutz is only implemented for
    %   fermionic Hilbert spaces.
    
    % ---------------------------------
    % Format for Gutz Modifier object:
    % - Gutz.G = on-site density suppression factor. 0 < G < 1 for fermions.
    % - Gutz.Den = (N x 1) vector of particle number per site.
    % - Gutz.Nmean = mean on-site particle density.
    % ---------------------------------
    % Format for Update is a vector of new density fluctations.
    % ---------------------------------
    % Format for dLogp is a single scalar, dG.
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
        Np = 1; % Number of parameters.
        G = 0.05; % Default starting value for G.
        Den = 1; % Placeholder for vector of density fluctuations.
        Nmean = 0; % Placeholder for mean on-site particle density.
    end
    
    properties (Hidden)
        OptInds = 1; % One parameter to vary.
    end
    
    methods
        % Constructor for GutzF Modifier object:
        function obj = Gutz(Hilbert,Params,VFlag)
            % Graph is second argument - used for other Modifiers but not
            % necessary here.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm') == 0
                error('Gutz Modifier subtype is only compatible with fermionic Hilbert objects.')
            end
            obj.G = Params.G; obj.Np = 1; obj.OptInds = VFlag;
            % Prepare placeholders for later initialisation in PrepPsi.
            obj.Den = zeros(Hilbert.N,1); obj.Nmean = 0;
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,~,dP)
            obj.G = obj.G + dP;
            if real(obj.G) < 0
                obj.G = 0 + imag(obj.G);
            elseif real(obj.G) > 1
                obj.G = 1 + imag(obj.G);
            end
        end
        % Update Modifier configuration information inside Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Den = Update;
        end
        
        % Initialise Modifier configuration information given a starting Cfg.
        function obj = PrepPsi(obj,Hilbert,Cfg)
            Cfg_vec = Hilbert.FullCfgRef(Cfg);
            PDen = sum(Cfg_vec(:))/Hilbert.N; obj.Nmean = PDen;
            DbOc = sum(Cfg_vec,2); obj.Den = DbOc; % dN determines whether Gutzwiller factor kicks in.
            % In fermionic case, Gutzwiller factor only relevant when there
            % are double occupancies.
        end
    end
    
    methods (Static)
        % Ratio of amplitudes for two configurations separated by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            % Convention used here is that PG = prod(1-gD(i)) i.e. a factor
            % of 1-g for every doubly occupied site.
            if Diff.type == 0
                % Fermion swap, no double occupancy changes occur from this.
                Ratio = 1; Update = obj.Den;
            elseif Diff.type == 2
                % Fermion pair move, no change in double occupancy number but change in position.
                Ratio = 1; Update = obj.Den; Update(Diff.pos(1)) = 0; Update(Diff.pos(2)) = 2;
                % First position is vacated site, second is destination.
            elseif Diff.type == 1 || Diff.type == -1
                % Single fermion move.
                Update = obj.Den;
                for d = 1:numel(Diff.pos)
                    Update(Diff.pos(d)) = Update(Diff.pos(d)) + Diff.val(d);
                end
                Ratio = (1-obj.G)^(sum(Update==2) - sum(obj.Den==2));
                if abs(Ratio)>1e10 
                    % Attempted fix for numerical instabilities in cases where G = 1.
                    Ratio = 1e10; 
                end
            end
        end
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,~,Cfg)
            % Convention used here is that PG = prod(1-gD(i)) i.e. a factor
            % of 1-g for every doubly occupied site.
            % Construct full configuration vector.
            Cfg_vec = Hilbert.FullCfgRef(Cfg);
            DbOc = sum(prod(Cfg_vec,2)); % Number of doubly occupied sites.
            % Only one parameter to optimise here, so dLogp is a scalar.           
            dLogp = -DbOc;
        end
    end
    
end