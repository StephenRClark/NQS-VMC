classdef Ferm < Hilbert
    % Ferm - Object containing details of a spinful fermionic configuration
    % Hilbert space.
    %   Hilbert is overarching class.
    
    properties
        Type = 'Ferm';
        N = 1; % Default to single site if no input.
        Nf = 2; % Populate with one fermion of each spin if no input.
        d = 4; % Default to single qubit if no input.
        Sector = [1 1]; % Fixed to one up and one down fermion if no input.
    end
    
    properties (Hidden)
        PropMove = @MoveFermCfg; % Default proposed configuration move.
        FullCfgRef = @FullFermCfg; % Default configuration vector for reference state.
        FullCfgMod = @BiPtFermCfg; % Default configuration vector for amplitude modifier.
        RandomCfgFun = @RandomFermZeroMag; % Default configuration generating function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Ferm(N,N_up,N_dn)
            if nargin == 1
                % Assume half-filling, zero magnetisation
                if mod(N,2) == 1
                    error('Require even number of sites for half-filled zero magnetisation configurations');
                else
                    obj.N = N; obj.Nf = N; obj.Sector = [N/2 N/2];
                end
            elseif nargin == 3
                obj.N = N; obj.Nf = N_up + N_dn; obj.Sector = [N_up N_dn];
            else
                % Revert to defaults.
                obj.N = 1; obj.Nf = 2; obj.Sector = [1 1];
            end
            if mod(obj.N,1)~=0
                error('Number of sites N should be integer.')
            end
            % Assign PropMove according to whether Sector is defined or not.
            if isempty(obj.Sector)
                obj.PropMove = @AddFermCfg;
            else
                obj.PropMove = @MoveFermCfg;
            end
            % Assign RandomCfgFun and starting CParams according to Sector.
            CParams.N = N;
            if isempty(obj.Sector)
                % If not specified, start in half-filled zero magnetisation case.
                obj.RandomCfgFun = @RandomFermZeroMag;
                CParams.Nf = N;
            elseif obj.Sector(1) == obj.Sector(2)
                obj.RandomCfgFun = @RandomFermZeroMag;
                CParams.Nf = sum(obj.Sector);
            else
                obj.RandomCfgFun = @RandomFermFixedPop;
                CParams.N_up = obj.Sector(1); CParams.N_dn = obj.Sector(2);
            end
            obj.CParams = CParams;
        end
    end
    
    methods (Static)
        % Call RandomCfg function according to sector.
        function [Cfg] = RandomCfg(obj)
            [Cfg] = obj.RandomCfgFun(obj.CParams);
        end
    end
end

