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
        Np = 3; % Number of parameters.
        Nv = 1; % Number of visible neurons.
        Nh = 1; % Number of hidden spins.
        a = 0; % Visible site bias terms, Nv x 1 vector.
        b = 0; % Hidden spin bias terms, Nh x 1 vector.
        W = 0; % Hidden-visible coupling terms, Nh x Nv matrix.
    end
    
    properties (Hidden)
        Theta = 0; % Local effective angle, Nh x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(3,1); % Individual parameter flags for variational purposes.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
        SFlag = 0; % Spin symmetry flag for fermionic case.
    end
    
    methods
        % Constructor for unconstrained NQS with no symmetries:
        function obj = NQS(Hilbert,~,Params,VFlag)
            % Graph necessary for NQS subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm')
                obj.Nv = 2*Hilbert.N; obj.FFlag = Hilbert.N;
                if (isempty(Hilbert.Sector) == 0) && (Hilbert.Sector(1) == Hilbert.Sector(2))
                    obj.SFlag = 1;
                end
                % For NQS, fermionic extension simply doubles the size of
                % the effective lattice, though if N_up = N_dn, spin
                % symmetry is imposed as well.
            else
                obj.Nv = Hilbert.N; obj.FFlag = 0;
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj = RandomInitPsiNQS(obj,Params);
        end
        
        % Update Modifier variational parameters according to changes dP.
        function obj = PsiUpdate(obj,Graph,dP)
            obj = PsiUpdateNQS(obj,Graph,dP);
        end
        
        % Update Modifier configuration information inside Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdateNQS(obj,Update);
        end
        
        % Initialise Modifier configuration information given a starting Cfg.
        function obj = PrepPsi(obj,Hilbert,Cfg)
            obj = PrepPsiNQS(obj,Hilbert,Cfg);
        end
        
    end
    
    methods (Static)
        
        % Ratio of amplitudes for two configurations separated by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            if obj.FFlag ~= 0 % Fermionic Diff conversion required.
                [Diff] = Ferm2BiPtDiff(Diff,obj.FFlag);
                % Ferm2BiPtDiff converts differences for pairs of
                % occupation numbers to differences in occupation on a
                % double sized lattice.
            end
            [Ratio,Update] = PsiRatioNQS(obj,Diff);
        end
        
        % Logarithmic derivative for the variational parameters in Modifier.
        function [dLogp] = LogDeriv(obj,Hilbert,~,Cfg)
            [dLogp] = LogDerivNQS(obj,Hilbert,Graph,Cfg);
        end
    end
    
end