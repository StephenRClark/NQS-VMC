classdef SDet < Reference
    % SDet - a Slater determinant fermionic state that populates the lowest
    % energy eigenstates of the provided quadratic fermionic Hamiltonian.
    %   Reference is overarching class. Variational versions invoke
    %   LogDeriv and PsiUpdate methods.
    
    % ---------------------------------
    % Format for SDet Reference:
    % - SDet.Nf = (1 x 2) vector - number of up/down fermions respectively.
    % - SDet.Orbitals = (2N x 2N) matrix - contains all available single particle orbitals.
    % - SDet.OrbMat = (2N x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
    % - SDet.DetMat = (2N x Nf) matrix - elements are used for determinants in PsiRatio.
    % - SDet.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
    % - SDet.Np = number of variational parameters associated with SDet Reference.
    % SDet properties used in variational version:
    % - SDet.CArr = (2N x 2N x Np) array - connectivity array for the reference Hamiltonian.
    % - SDet.WArr = (2N x 2N x Np) array - transformed connectivity array for the reference Hamiltonian.
    % - SDet.EnFac = (2N x 2N) matrix - elements are used in LogDeriv function.
    % - SDet.HVar = (Np x 1) vector - variational parameters in the reference Hamiltonian used.
    % ---------------------------------
    % Format for Update is a struct containing updates for DetMat and FermLoc.
    % ---------------------------------
    % Format for dLogp is a Np x 1 vector of parameter derivatives.
    % ---------------------------------
    % N.B. First parameter in HVar will automatically be excluded from the
    % optimisation and act as the energy scale of the reference, as a
    % preventative measure against runaway HVar scaling. It is recommended
    % to set the first term in HVar as the hopping term.
    % ---------------------------------
    
    properties
        Type = 'SDet'; % Identifier for the reference state.
        VFlag = 0; % Variational flag, deactivated by default.
        N = 1; % Number of sites, inherited from input Hilbert.
        Nf = [1 1]; % Number of up and down fermions, inherited from input Hilbert.
        Orbitals = eye(2); % All eigenstates of the provided Hamiltonian.
        % Eigenstates are arranged by 'spin' (i.e. if orbital has greater
        % probability of being occupied by an up or down fermion) then
        % energy in ascending order.
        Np = 0; % Number of variational parameters associated with SDet if made variational.
    end
     
    properties (Hidden)
        ParamCap = 10; % Parameter magnitude cap.
        OrbMat = eye(2); % Reduced orbital matrix used for computing ratios.
        DetMat = eye(2); % Matrix of determinant ratios calculated from OrbMat and Orbitals.
        FermLoc = [1 2]; % Details locations of fermions by index for sign tracking purposes.
        % Fermions and their creation operators are listed from left to
        % right in order of ascending site index, with down fermions
        % occupying 'sites' N+1 -> 2N.
        
        % Following hidden properties only get called in variational versions.
        CArr = eye(2); % 2N x 2N x Np array used to reconstruct Hamiltonian if variational.
        WArr = eye(2); % 2N x 2N x Np array transformed from CArr for calculating derivatives.
        HVar = []; % Np x 1 vector of variational terms in the quadratic Hamiltonian.
        EnFac = eye(2); % 2N x 2N matrix of energy factors used for calculating derivatives.
        OptInds = 0; % Np x 1 vector - if entry n is 1, parameter derivative for p will be calculated.
    end
    
    methods
        % Constructor for SDet Reference object:
        function obj = SDet(Hilbert,Params,VFlag)
            % Required fields in Params:
            % HVar - Np x 1 vector of parameters desired in quadratic Hamiltonian.
            % CArr - 2N x 2N x Np array detailing connectivity of each parameter.
            if nargin == 2
                disp('No variational flag specified - assuming fixed Slater determinant reference.')
                obj.Np = 0;
            elseif nargin == 3
                obj.Np = numel(Params.HVar); obj.VFlag = VFlag;
            end
            N = Hilbert.N; obj.N = Hilbert.N;
            if isempty(Hilbert.Sector)
                % Placeholder N_up and N_dn if no Sector specified.
                Nf = [round(N/2) round(N/2)];
            else
                Nf = Hilbert.Sector;
            end
            obj.Nf = Nf;
            [obj] = FixedInitPsiSDet(obj,Params);            
        end
        
        % Initialise Reference configuration values given a starting Cfg.
        function [obj] = PrepPsi(obj,Hilbert,Cfg)
            obj = PrepPsiSDet(obj,Hilbert,Cfg);
        end
        
        % Update Reference configuration information according to Update.
        function [obj] = PsiCfgUpdate(obj,Update)
            obj = PsiCfgUpdateSDet(obj,Update);
        end
        
        % Update Reference variational parameters according to changes dP.
        function obj = PsiUpdate(obj,~,dP)
            obj = PsiUpdateSDet(obj,0,dP);
        end
    end
    
    methods (Static)
        % Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            [Ratio,Update] = PsiRatioSDet(obj,Diff);
        end
        
        % Logarithmic derivative for the variational parameters in Reference.
        function [dLogp] = LogDeriv(obj,Hilbert,~,Cfg)
            [dLogp] = LogDerivSDet(obj,Hilbert,Cfg);
        end
    end
    
end