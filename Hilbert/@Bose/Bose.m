classdef Bose < Hilbert
    % Bose - Object containing details of a spinless bosonic configuration
    % Hilbert space.
    %   Hilbert is overarching class.
    
    properties (SetAccess = protected)
        Type = 'Bose';
        N = 1; % Default to single site if no input.
        Nb = 1; % Populate with one boson if no input.
        d = 2; % Default to single qubit if no input.
        Sector = []; % Empty sector - assume no fixed boson number.
    end
    
    properties (Hidden, SetAccess = protected)
        PropMoveFunc = @MoveBoseCfg; % Default proposed configuration move.
        FullCfgFunc = @FullBoseCfg; % Default configuration vector for reference state.
        RandomCfgFunc = @RandomBoseZeroSpin; % Default configuration generating function.
        Diff2CfgFunc = @Diff2CfgBose; % Default difference to configuration conversion function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Bose(N,Nb,Nmax)
            if nargin == 3
                % Third input is maximum permitted number of bosons per site.
                obj.N = N; obj.Nb = Nb; obj.Sector = Nb; obj.d = Nmax+1;
                if Nb/N > Nmax
                    error('Requested boson density is higher than requested tolerance.')
                end
            elseif nargin == 3 && isempty(Nb)
                % Empty Nb suggests no fixed number sector i.e. a grand
                % canonical calculation. Third argument is Nmax.
                obj.N = N; obj.Nb = Nb; obj.Sector = []; obj.d = Nmax + 1;
            elseif nargin == 2
                % Assume input Nb is the desired fixed number sector. Set
                % Nmax to be average density plus four.
                obj.N = N; obj.Nb = Nb; obj.Sector = Nb; obj.d = ceil(Nb/N) + 5; Nmax = ceil(Nb/N) + 4;
            elseif nargin == 1
                % If no Nb specified, assume starting density of 1 and no
                % fixed number. Set Nmax to be five.
                obj.N = N; obj.Nb = N; obj.Sector = []; obj.d = 6; Nmax = 5;
            else
                % Revert to defaults with no inputs;
                obj.N = 1; obj.Nb = 1; obj.Sector = []; obj.d = 2; Nmax = 1;
            end
            if mod(obj.N,1)~=0
                error('Number of sites N should be integer.')
            end
%             if obj.N == obj.Nb
%                 obj.RandomCfgFunc = @RandomBoseUnifFill;
%             end
            % Assign PropMove and CParams according to whether Sector is defined or not.
            CParams.N = N;
            if isempty(obj.Sector)
                obj.PropMoveFunc = @AddBoseCfg;
                CParams.Nb = N; % Start at uniform density if no fixed density.
            else
                obj.PropMoveFunc = @MoveBoseCfg;
                CParams.Nb = obj.Sector;
            end
            CParams.Nmax = Nmax; obj.CParams = CParams;
            % Reference and modifiers always use same configuration vector for bosons.
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
            Properties.N = obj.N; Properties.Nb = obj.Nb; 
            Properties.Nmax = obj.d-1; Properties.Sector = obj.Sector;
        end
        
    end
    
end