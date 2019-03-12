classdef Bose < Hilbert
    % Bose - Object containing details of a spinless bosonic configuration 
    % Hilbert space.
    %   Hilbert is overarching class.
    
    properties
        Type = 'Bose'; 
        N = 1; % Default to single site if no input.
        Nb = 1; % Populate with one boson if no input.
        d = 2; % Default to single qubit if no input.
        Sector = []; % Empty sector - assume no fixed boson number.        
    end
    
    properties (Hidden)        
        PropMove = @MoveBoseCfg; % Default proposed configuration move.
        FullCfgRef = @FullBoseCfg; % Default configuration vector for reference state.
        FullCfgMod = @FullBoseCfg; % Default configuration vector for amplitude modifier.
        RandomCfgFun = @RandomBoseZeroSpin; % Default configuration generating function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Bose(N,Nb,Nmax)
            if nargin == 3 
                % Third input is maximum permitted number of bosons per site.
                obj.N = N; obj.Nb = Nb; obj.Sector = Nb; obj.d = Nmax+1; 
                if Nb/N >= Nmax
                    error('Requested boson density is higher than requested tolerance.')
                end
            elseif nargin == 2
                % Assume input Nb is the desired fixed number sector.
                obj.N = N; obj.Nb = Nb; obj.Sector = Nb; obj.d = Nb+1; Nmax = Nb;
            elseif nargin == 1 
                % If no Nb specified, assume density of 1 and no fixed number.
                obj.N = N; obj.Nb = N; obj.Sector = []; obj.d = obj.Nb+1; Nmax = N;
            else
                % Revert to defaults with no inputs;
                obj.N = 1; obj.Nb = 1; obj.Sector = []; obj.d = 2; Nmax = 1;
            end
            if mod(obj.N,1)~=0
                error('Number of sites N should be integer.')
            end
            % Assign PropMove and CParams according to whether Sector is defined or not.
            CParams.N = N;
            if isempty(obj.Sector)
                obj.PropMove = @AddBoseCfg;
                CParams.Nb = N; % Start at uniform density if no fixed density.
                Nmax = N^2; % Some hign peak occupation to allow for plenty of bosons to be added.
            else
                obj.PropMove = @MoveBoseCfg;
                CParams.Nb = obj.Sector; 
            end
            CParams.Nmax = Nmax; obj.CParams = CParams;
            % Reference and modifiers always use same configuration vector for bosons.
        end        
    end
    
    methods (Static)
        % Call RandomCfg function according to sector.
        function [Cfg] = RandomCfg(obj)
            [Cfg] = obj.RandomCfgFun(obj.CParams);
        end
    end
end

