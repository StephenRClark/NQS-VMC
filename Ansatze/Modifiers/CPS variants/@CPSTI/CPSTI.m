classdef CPSTI < Modifier
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier. Can optimise.
    end

    properties (SetAccess = protected)
        Nv = 1; % Number of visible spins.
        Nh = 1; % Number of hidden spins.
        Np = 8; % Number of parameters. Np = Nv * 2 + Nh * 2 + Nv * Nh * 4
        a = 0; % visible bias array. 2 * Nv elements
        b = 0; % hidden bias array. 2 * Nh elements
        W = 0; % couplings array. 4 * Nv * Nh elements
        Graph % Symmetry purpose.
        ati = 0; % visible bias translational invariant. 2 * 1 elements
        bti = 0; % hidden bias translational invariant. 2 * alpha elements
        Wti = 0; % couplings translational invariant. 4 * Nv * alpha. Nv = Nh.
        Alpha = 1;
    end

    properties (Hidden)
        tempCfgVec = 0;
        Theta = 0; % 2 * Nh elements
        ParamCap = 5;
        OptInds = zeros(8,1);
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
        SFlag = 0; % Spin symmetry flag for fermionic case.
    end

    methods
        function obj = CPSTI(Hilbert,Graph,Params,VFlag)
            % Graph necessary for CPS subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified. What is nargin
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm')
                obj.Nv = 2*Hilbert.N; obj.FFlag = Hilbert.N; obj.FullCfg = @BiPtFermCfg;
                if (isempty(Hilbert.Sector) == 0) && (Hilbert.Sector(1) == Hilbert.Sector(2))
                    obj.SFlag = 1;
                end
                % For CPS, fermionic extension simply doubles the size of
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
            if isfield(Params,'Alpha') %What
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj.Graph = Graph;
            obj = initCPSTI(obj,Params);
        end

        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = logDerivCPSTI(obj,Cfg);
        end

        function obj = PrepPsi(obj,Cfg)
            obj = prepCPSTI(obj,Cfg);
        end

        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta; obj.tempCfgVec = Update.tempCfgVec;
        end

        function [Params] = ParamList(obj)
            Params = paramListCPSTI(obj);
        end

        function [Psi] = PsiGenerate(obj,Basis)
            Psi = zeros(size(Basis,1),1);
            for k = 1:size(Basis,1)
                Cfg_vec = Basis(k,:);
                Psi(k) = generateCPSTI(obj,Cfg_vec);
            end
            ModPsi = sqrt(sum(abs(Psi).^2));
            Psi = Psi/ModPsi;
        end

        function obj = PsiUpdate(obj,dP)
            obj = updateCPSTI(obj,dP);
        end

        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = ratioCPSTI(obj,Diff);
        end
    end
end
