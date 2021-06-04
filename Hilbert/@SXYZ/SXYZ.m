classdef SXYZ < Hilbert
    % SXYZ - object containing details of a spin-1 configuration Hilbert
    % space in the xyz basis with the numerical mapping:
    % x -> +1, y -> 0, z -> -1.
    %   Hilbert is overarching class.
    
    properties (SetAccess = protected)
        Type = 'Spin';
        N = 1; % Default to single site if no input.
        S = 1; % Default to spin-1
        d = 3; % Single site Hilbert space dimension.
        Sector = 0; % Sector determines 'parity' of the system (0 = even, 1 = odd).
        % If Sector is empty, then no parity is preserved.
    end
    
    properties (Hidden, SetAccess = protected)
        PropMoveFunc = @SwapSxyzCfg; % Default proposed configuration move.
        FullCfgFunc = @FullSpinCfg; % Default configuration vector for reference state.
        RandomCfgFunc = @RandomSxyzFixedParity; % Default random configuration generating function.
        Diff2CfgFunc = @Diff2CfgSpin1; % Default difference to configuration conversion function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = SXYZ(N,Sector)
            % Currently only spin-1/2 and spin-1 are supported.
            if nargin == 1 % Default to no Sector if not specified.
                Sector = [];
            end
            if Sector > 1
                disp('Sector should be 0 or 1 (even or odd) - taking floor and mod of provided value.');
                Sector = floor(mod(Sector,2));
            end
            obj.Sector = Sector;
            if isempty(obj.Sector)
                obj.PropMoveFunc = @FlipSpin1Cfg;
                obj.RandomCfgFunc = @RandomSpin1NoSector;
            end
            CParams.N = N; CParams.Sector = Sector;
            obj.CParams = CParams; obj.N = N;
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
        
        % ChangePropMove: Change an existing Hilbert object's peoposed move
        % function to a new compatible one. Requires preservation of
        % existing symmetry constraints.
        function obj = ChangePropMove(obj,FuncHandle)
            FuncStr = func2str(FuncHandle);
            if isempty(obj.Sector)
                if strcmp(FuncStr,'FlipSpin1Cfg')
                    obj.PropMoveFunc = FuncHandle;
                    disp(['Proposed move function changed to ' FuncStr '.']);
                elseif strcmp(FuncStr,'FlipMultiSpin1')
                    obj.PropMoveFunc = FuncHandle;
                    disp(['Proposed move function changed to ' FuncStr '.']);
                else
                    disp('Invalid proposed move function - no change applied.');
                end
            else
                if strcmp(FuncStr,'MoveSpin1Cfg')
                    obj.PropMoveFunc = FuncHandle;
                    disp(['Proposed move function changed to ' FuncStr '.']);
                elseif strcmp(FuncStr,'MoveMultiSpin1')
                    obj.PropMoveFunc = FuncHandle;
                    disp(['Proposed move function changed to ' FuncStr '.']);
                else
                    disp('Invalid proposed move function - no change applied.');
                end
            end
        end
        
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'SXYZ';
            Properties.N = obj.N; Properties.Sector = obj.Sector;
        end
    end
    
end