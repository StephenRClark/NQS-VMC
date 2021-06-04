classdef Pfaf < Reference
    % Pfaf - a fermionic reference state composed of pairing amplitudes
    % that can be individually varied. The configuration amplitudes are
    % Pfaffians of the pairing amplitude matrix.
    %   Reference is overarching class. Variational versions invoke
    %   LogDeriv and PsiUpdate methods. Subvariants include implicit
    %   symmetries.
    
    % ---------------------------------
    % Format for fermionic Pfaffian Reference:
    % - Pfaf.Nf = (1 x 2) vector - number of up/down fermions respectively.
    % - Pfaf.PairMat = (2N x 2N) matrix - contains all pairing terms.
    % - Pfaf.PfI = (Nf x Nf) matrix - inverse of reduced PfFull matrix.
    % - Pfaf.PfG = (2N x Nf) matrix - matrix used for ratio calculations.
    % - Pfaf.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
    % - Pfaf.Np = number of variational parameters associated with Pfaf Reference.
    % Pfaf properties used in variational version:
    % - Pfaf.PfV = (2N x 2N) array - logs which variational parameters make up the elements of Pfaf.PairMat.
    % - Pfaf.PfVR = (Nf x Nf) array - reduced matrix constructed from PfV.
    % - Pfaf.PfVar = (Np x 1) vector - variational parameters in PfFull.
    % ---------------------------------
    % Format for Update is a struct containing updates for PfI, PfG and FermLoc.
    % ---------------------------------
    % Format for dLogp vector is a (Np x 1) vector of parameter derivatives.
    % ---------------------------------
    
    properties
        VFlag = 1; % Variational flag, activated by default.
    end
    
    properties (SetAccess = protected)
        Type = 'Ferm'; % Identifier for the reference state.
        N = 1; % Number of sites, inherited from input Hilbert.
        Np = 1; % Number of parameters associated with Pfaf if made variational.
        Nf = [1 1]; % Number of up and down fermions, inherited from input Hilbert.
        PairMat = [0 1; -1 0]; % Antisymmetric matrix of pairing amplitudes.
        Graph % Details connectivity of sites - used to include symmetries.
    end
    properties (Hidden)
        ParamCap = 10; % Parameter magnitude cap.
        OptInds = 1; % Np x 1 vector - if entry n is 1, parameter derivative for p will be calculated.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullFermCfg; % Function used by Reference to interface with Cfg structs.
        PfI = [0 1; -1 0]; % Placeholder for reduced pairing amplitude matrix inverse, used in ratios.
        PfG = eye(2); % Placeholder for modified matrix used for ratios.
        FermLoc = [1 2]; % Details locations of fermions by index for sign tracking purposes.
        % Fermions and their creation operators are listed from left to
        % right in order of ascending site index, with down fermions
        % occupying 'sites' N+1 -> 2N.
        
        % Following hidden properties only get called in variational versions.
        PfV = zeros(2,2); % 2N x 2N array used to reconstruct PairMat if variational.
        PfVR = zeros(2,2); % Nf x Nf array constructed from PfV.
        PfVar = 1; % Np x 1 vector of variational parameters in the pairing matrix.
    end
    
    methods
        % Constructor for Pfaf Reference object:
        function obj = Pfaf(Hilbert,Graph,Params,VFlag)
            % Required fields in Params:
            % PfVar - Np x 1 vector of parameters desired in pairing amplitude matrix.
            % PfV - 2N x 2N array detailing parameters associated with each pairing term.
            if nargin == 3
                disp('No variational flag specified - assuming fixed Slater determinant reference.')
                obj.Np = 0;
            elseif nargin == 4
                obj.Np = numel(Params.PfVar); obj.VFlag = VFlag;
            end
            N = Hilbert.N; obj.N = N;
            if isempty(Hilbert.Sector)
                % Placeholder N_up and N_dn if no Sector specified.
                Nf = [round(N/2) round(N/2)];
            else
                Nf = Hilbert.Sector;
            end
            obj.Nf = Nf;
            
            % Pfaf can be initialised in two ways - from the eigenstates of
            % a quadratic Hamiltonian, and with random parameters with
            % symmetries imposed by Graph.
            % First option will initialise with FixedInitPsiPfaf, though if
            % variational, will then take the individual terms and optimise
            % without regard for initial Hamiltonian parameters.
            % Second option will require Graph with desired symmetry, and
            % will populate a matrix with small random parameters.
            % Both have same format and utilise same LogDerivPfaf and
            % PsiUpdate functions - only difference is initialisation.
            
            if isfield(Params,'CArr') && isfield(Params,'HVar') % Connectivity array used in reference Hamilonian ansatze.
                [obj] = FixedInitPsiPfaf(obj,Params);
            else % Will assume random parameters subject to Graph symmetries if no CArr present.
                [obj] = RandomInitPsiPfaf(obj,Graph,Params);
            end
        end
        
        % PrepPsi: Initialise Reference configuration values given a starting Cfg.
        function [obj] = PrepPsi(obj,Hilbert,Cfg)
            obj = PrepPsiPfaf(obj,Hilbert,Cfg);
        end
        
        % PsiCfgUpdate: Update Reference configuration information according to Update.
        function [obj] = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdatePfaf(obj,Update);
        end
        
        % PsiUpdate: Update Reference variational parameters according to changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdatePfaf(obj,dP);
        end
        
        % RndBatchSelect: Randomly select some proportion of parameters to
        % optimise and disable optimisation of the rest. Used for random
        % batch optimisation schemes.
        function [obj] = RndParamSelect(obj,OFrac)
            NpR = round(min(max(1,OFrac*obj.Np),obj.Np)); % Ensure at least 1 or at most Np are selected.
            OptIndsP = zeros(obj.Np,1); OptIndsP(randperm(obj.Np,NpR,1)) = 1;
            obj.OptInds = OptIndsP;
        end
        
        % ParamList; outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = obj.PfVar;
        end
        
        % PsiRatio: Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = PsiRatioPfaf(obj,Diff);
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters in Reference.
        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = LogDerivPfaf(obj,Cfg);
        end
    end
    
end