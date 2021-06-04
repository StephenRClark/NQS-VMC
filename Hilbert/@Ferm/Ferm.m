classdef Ferm < Hilbert
    % Ferm - Object containing details of a spinful fermionic configuration
    % Hilbert space.
    %   Hilbert is overarching class.
    
    properties (SetAccess = protected)
        Type = 'Ferm';
        N = 1; % Default to single site if no input.
        Nf = 2; % Populate with one fermion of each spin if no input.
        Sz = 0; % Total z-spin projection.
        d = 4; % Default to single qubit if no input.
        Sector = [1 1]; % First entry signifies fixed particle number, second is for fixed projection.
    end
    
    properties (Hidden, SetAccess = protected)
        PropMoveFunc = @MoveFermCfg; % Default proposed configuration move.
        FullCfgFunc = @FullFermCfg; % Default configuration vector for reference state.
        RandomCfgFunc = @RandomFermZeroMag; % Default configuration generating function.
        Diff2CfgFunc = @Diff2CfgFerm; % Default difference to configuration conversion function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Ferm(N,Nferm,SzTotal)
            if nargin < 2
                Nferm = [];
            end
            if nargin < 3
                SzTotal = [];
            end
            if mod(N,1)~=0
                error('Number of sites N should be integer.');
            end
            if nargin == 3
                if mod(Nferm+SzTotal,2) ~= 0
                    error('Number of fermions and total spin projection are incompatible.');
                end
            end            
            % Assign RandomCfgFun and starting CParams according to Sector.
            obj.N = N; CParams.N = N; 
            if isempty(Nferm) % Fermion number not fixed.
                obj.RandomCfgFunc = @RandomFermZeroMag; 
                obj.Sector = [0 0]; obj.Sz = 0;
                CParams.Nf = N; obj.PropMoveFunc = @AddFermCfg;
            elseif isempty(SzTotal) % Spin projection not fixed.
                obj.RandomCfgFunc = @RandomFermZeroMag;
                obj.Sector = [1 0]; obj.Sz = 0; CParams.Nf = Nferm;
            else % Spin projection and number fixed.
                obj.RandomCfgFunc = @RandomFermFixedPop; 
                obj.Sector = [1 1]; obj.Sz = SzTotal;
                CParams.N_up = 0.5*(Nferm + SzTotal); CParams.N_dn = 0.5*(Nferm - SzTotal);
            end
            obj.CParams = CParams;
        end
        
        % RandomCfg: Generate Cfg struct subject to Hilbert parameters.
        function [Cfg] = RandomCfg(obj)
            [Cfg] = obj.RandomCfgFunc(obj.CParams);
        end
        
        % PropMove: Given a Cfg struct, generate a proposed CfgP and
        % corresponding difference struct Diff.
        function [Diff,CfgP] = PropMove(obj,Cfg)
            [Diff,CfgP] = obj.PropMoveFunc(Cfg);
        end
        
        % FullCfg: Given a Cfg struct, create a vector representation of
        % the configuration.
        function [Cfg_vec] = FullCfg(obj,Cfg)
            [Cfg_vec] = obj.FullCfgFunc(Cfg);
        end
        
        % Diff2Cfg: Given a Cfg and Diff struct, combine to create a new
        % CfgP struct.
        function [CfgP] = Diff2Cfg(obj,Diff,Cfg)
            [CfgP] = obj.Diff2CfgFunc(Diff,Cfg);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.N = obj.N; Properties.Nf = obj.Nferm; 
            Properties.Sz = obj.Sz; Properties.Sector = obj.Sector;
        end        
    end
    
end