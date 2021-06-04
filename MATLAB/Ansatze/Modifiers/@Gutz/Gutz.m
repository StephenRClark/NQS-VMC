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
        G = 0.05; % Default starting value for G.
    end
    
    properties (SetAccess = protected)
        Np = 1; % Number of parameters.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        OptInds = 1; % One parameter to vary.
        ParamCap = 100; % Parameter cap.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullFermCfg; % Function used by Modifer to interface with Cfg structs.
        RatioFunc = @PsiRatioGFerm; % PsiRatio function handle - default to fermionic case.
        PrepFunc = @PrepPsiGFerm; % PrepPsi function handle - default to fermionic case.
        DerivFunc = @LogDerivGFerm; % LogDeriv function handle - default to fermionic case.
        Den = 1; % Placeholder for vector of density fluctuations.
        Nmean = 0; % Placeholder for mean on-site particle density.
    end
    
    methods
        % Constructor for GutzF Modifier object:
        function obj = Gutz(Hilbert,Graph,Params,VFlag)
            % Graph is second argument - used for other Modifiers but not
            % necessary here.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Spin') == 1
                error('Gutz Modifier subtype is not compatible with spin Hilbert objects.')
            end
            obj.G = Params.G; obj.Np = 1; obj.OptInds = VFlag;
            % Prepare placeholders for later initialisation in PrepPsi.
            obj.Den = zeros(Hilbert.N,1); obj.Nmean = 0; obj.Graph = Graph;
            if strcmp(Hilbert.Type,'Ferm') == 1
                obj.DerivFunc = @LogDerivGFerm;
                obj.PrepFunc = @PrepPsiGFerm;
                obj.RatioFunc = @PsiRatioGFerm;
                obj.FullCfg = @FullFermCfg;
            elseif strcmp(Hilbert.Type,'Bose') == 1
                obj.DerivFunc = @LogDerivGBose;
                obj.PrepFunc = @PrepPsiGBose;
                obj.RatioFunc = @PsiRatioGBose;
                obj.FullCfg = @FullBoseCfg;
            end
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj.G = obj.G + dP;
            if real(obj.G) > obj.ParamCap
                obj.G = obj.ParamCap + imag(obj.G);
            end
            if real(obj.G) < 0
                obj.G = 0 + imag(obj.G);
            end
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Den = Update;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            [obj] = obj.PrepFunc(obj,Cfg);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio, Update] = obj.RatioFunc(obj,Diff);
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = obj.DerivFunc(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised wavefunction amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            if strcmp(func2str(obj.FullCfg),'FullFermCfg')
                [Psi] = PsiGenerateGFerm(obj,Basis);
            elseif strcmp(func2str(obj.FullCfg),'FullBoseCfg')
                [Psi] = PsiGenerateGBose(obj,Basis);
            end
        end
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = obj.G;
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj.G = P;
            if real(obj.G) > obj.ParamCap
                obj.G = obj.ParamCap + imag(obj.G);
            end
            if real(obj.G) < 0
                obj.G = 0 + imag(obj.G);
            end
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            if obj.FullCfg == @FullBoseCfg
                Properties.Type = 'GutzB';
            else
                Properties.Type = 'GutzF';
            end
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds; 
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
    
end

% --- Fermionic Gutzwiller wave function amplitude ratio function ---

function [Ratio, Update] = PsiRatioGFerm(obj,Diff)
% Convention used here is that PG = exp(-gD(i)) i.e. a factor
% of exp(-g) for every doubly occupied site.
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
    Ratio = exp(-obj.G * (sum(Update==2) - sum(obj.Den==2)));
    if abs(Ratio)>1e10
        % Attempted fix for numerical instabilities in cases where G = 1.
        Ratio = 1e10*sign(Ratio);
    end
end
end

% --- Fermionic Gutzwiller wave function preparation function ---

function [obj] = PrepPsiGFerm(obj,Cfg)
Cfg_vec = obj.FullCfg(Cfg);
PDen = sum(Cfg_vec(:))/obj.Graph.N; obj.Nmean = PDen;
DbOc = sum(Cfg_vec,2); obj.Den = DbOc; % dN determines whether Gutzwiller factor kicks in.
% In fermionic case, Gutzwiller factor only relevant when there are double occupancies.
end

% --- Fermionic Gutzwiller logarithmic derivative function ---

function [dLogp] = LogDerivGFerm(obj,Cfg)
% Convention used here is that PG = exp(-gD(i)) i.e. a factor
% of exp(-g) for every doubly occupied site.
% Construct full configuration vector.
Cfg_vec = obj.FullCfg(Cfg);
DbOc = sum(prod(Cfg_vec,2)); % Number of doubly occupied sites.
% Only one parameter to optimise here, so dLogp is a scalar.
dLogp = -DbOc;
end

% --- Bosonic Gutzwiller wave function amplitude ratio function ---

function [Ratio, Update] = PsiRatioGBose(obj,Diff)
% Convention used here is that PG = exp(-g(n-Nmean)^2) i.e. a factor
% of exp(-g) for every fluctuation from mean occupation.
Update = obj.Den; 
for d = 1:numel(Diff.pos)
    Update(Diff.pos(d)) = Update(Diff.pos(d)) + Diff.val(d);
end
Ratio = exp(-obj.G * ( sum((Update.^2) - (obj.Den.^2)) ) );
if abs(Ratio)>1e10
    % Attempted fix for numerical instabilities.
    Ratio = 1e10*sign(Ratio);
end

end

% --- Bosonic Gutzwiller wave function preparation function ---

function [obj] = PrepPsiGBose(obj,Cfg)
Cfg_vec = obj.FullCfg(Cfg);
PDen = sum(Cfg_vec(:))/obj.Graph.N; obj.Nmean = PDen;
obj.Den = Cfg_vec - obj.Nmean; % dN determines whether Gutzwiller factor kicks in.
% In bosonic case, Gutzwiller factor affects deviations from mean density.
end

% --- Bosonic Gutzwiller logarithmic derivative function ---

function [dLogp] = LogDerivGBose(obj,Cfg)
% Convention used here is that PG = exp(-g(n-Nmean)^2) i.e. a factor
% of exp(-g) for every fluctuation from mean occupation.
% Construct full configuration vector.
Cfg_vec = obj.FullCfg(Cfg);
DenFluc = sum((Cfg_vec - obj.Nmean).^2);
dLogp = -DenFluc;
end

% --- Exact normalised fermionic Gutzwiller amplitude generating function ---

function [Psi] = PsiGenerateGFerm(obj,Basis)
N = size(Basis,2)/2;
Basis = Basis(:,1:N) + Basis(:,N+(1:N));
Psi = exp(-obj.G * sum(Basis==2,2));
Psi = Psi / sum(abs(Psi).^2);
end

% --- Exact normalised bosonic Gutzwiller amplitude generating function ---

function [Psi] = PsiGenerateGBose(obj,Basis)
N = size(Basis,2); Nmean = sum(Basis(1,:))/N;
Basis = Basis - Nmean;
Psi = exp(-obj.G * sum(Basis.^2,2));
Psi = Psi / sum(abs(Psi).^2);
end