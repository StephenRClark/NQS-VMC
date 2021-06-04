classdef CPS
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier. Can optimise.
    end

    properties (SetAccess = protected)
        Nv = 1; % Number of visible bias.
        Nh = 1; % Number of hidden bias.
        Np = 8; % Number of parameters. Np = Nv * 2 + Nh * 2 + Nv * Nh * 4
        a = 0;
        b = 0;
        W = 0;
        Graph % Symmetry purpose.
        tempCfgVec = 0;
        Theta = 0;
    end

    properties (Hidden)
        ParamCap = 5;
        OptInds = zeros(8,1);
    end

    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
        SFlag = 0; % Spin symmetry flag for fermionic case.
    end

    methods
        function obj = CPS(Hilbert,Graph,Params,VFlag)
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
            obj = RandomInitPsiCPS(obj,Params);
        end

        function [Params] = ParamList(obj)
            Params = ParamListCPS(obj);
        end
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta;
            obj.tempCfgVec = Update.tempCfgVec;%Check with Michael tomorrow
        end
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiCPS(obj,Cfg);
        end
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateCPS(obj,Basis);
        end
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = ratioCPS(Diff);
        end
    end
end
