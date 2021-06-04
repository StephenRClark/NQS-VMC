classdef Spin < Hilbert
    % Spin - object containing details of a spin configuration Hilbert
    % space.
    %   Hilbert is overarching class.
    
    properties (SetAccess = protected)
        Type = 'Spin';
        N = 1; % Default to single site if no input.
        S = 1/2; % Default to spin-1/2 if no input.
        d = 2; % Single site Hilbert space dimension.
        Sector = []; % Empty sector - assume all Sz sectors permitted if no input.
    end
    
    properties (Hidden, SetAccess = protected)
        PropMoveFunc = @FlipSpinCfg; % Default proposed configuration move.
        FullCfgFunc = @FullSpinCfg; % Default configuration vector for reference state.
        RandomCfgFunc = @RandomSpinZeroMag; % Default random configuration generating function.
        Diff2CfgFunc = @Diff2CfgSpinH; % Default difference to configuration conversion function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Spin(N,S,Sector)
            % Currently only spin-1/2 and spin-1 are supported.
            if nargin == 1 % Default to spin-1/2 if no S defined.
                obj.S = 1/2;
                if mod(N,1)
                    error('Number of sites N should be integer.')
                else
                    obj.N = N; obj.Sector = [];
                end
            elseif nargin == 2 % Default to no Sz projection if no Sector input.
                obj.S = S; obj.N = N; obj.Sector = [];
            elseif nargin == 3
                obj.S = S; obj.N = N; obj.Sector = Sector;
            elseif nargin == 0
                obj.S = 1/2; obj.N = 1; obj.Sector = [];
            end
            if obj.S ~= 1/2 && obj.S ~= 1
                error('Only spin-1/2 and spin-1 are currently supported. Ensure the second argument is 1/2 or 1.')
            end
            obj.d = 2*obj.S + 1;
            % Assign PropMove according to whether sector is defined or not.
            if isempty(obj.Sector)
                if obj.S == 1/2
                    obj.PropMoveFunc = @FlipSpinCfg;
                elseif obj.S == 1
                    obj.PropMoveFunc = @FlipSpin1Cfg;
                end
            else
                if obj.S == 1/2
                    obj.PropMoveFunc = @SwapSpinCfg;
                elseif obj.S == 1
                    obj.PropMoveFunc = @MoveSpin1Cfg;
                end
            end
            % Assign RandomCfgFun and starting CParams according to S and
            % Sector respectively.
            if obj.S == 1/2
                obj.RandomCfgFunc = @RandomSpinHFixedMag;
                obj.Diff2CfgFunc = @Diff2CfgSpinH;
            elseif obj.S == 1
                obj.RandomCfgFunc = @RandomSpin1FixedMag;
                obj.Diff2CfgFunc = @Diff2CfgSpin1;
            end
            CParams.N = N;
            if isempty(obj.Sector)
                if obj.S == 1
                    obj.RandomCfgFunc = @RandomSpin1NoSector;
                elseif obj.S == 1/2
                    obj.RandomCfgFunc = @RandomSpinHalf;
                    CParams.SzT = 0;
                end
            else
                CParams.SzT = obj.Sector;
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
        
        % ChangePropMove: Change an existing Hilbert object's peoposed move
        % function to a new compatible one. Requires preservation of
        % existing symmetry constraints.
        function obj = ChangePropMove(obj,FuncHandle)
            FuncStr = func2str(FuncHandle);
            if obj.S == 1/2
                if isempty(obj.Sector)
                    if strcmp(FuncStr,'FlipSpinCfg')
                        obj.PropMoveFunc = FuncHandle;
                        disp(['Proposed move function changed to ' FuncStr '.']);
                    elseif strcmp(FuncStr,'FlipMultiSpin')
                        obj.PropMoveFunc = FuncHandle;
                        disp(['Proposed move function changed to ' FuncStr '.']);
                    elseif strcmp(FuncStr,'FlipMultiSpinGS')
                        obj.PropMoveFunc = FuncHandle;
                        disp(['Proposed move function changed to ' FuncStr '.']);
                    else
                        disp('Invalid proposed move function - no change applied.');
                    end
                end
            elseif obj.S == 1
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
        end
        
        % ChangeStartCfg: Change an existing Hilbert object's proposed move
        % function to a new compatible one. Need both function handle and
        % configuration parameters.
        function [obj] = ChangeStartCfg(obj,FuncHandle,CParams)
            FuncStr = func2str(FuncHandle); obj.CParams = CParams;
            if obj.S == 1/2
                if isempty(obj.Sector)
                    if strcmp(FuncStr,'RandomSpinHFixedMag')
                        obj.RandomCfgFunc = FuncHandle;
                        disp(['Starting configuration function changed to ' FuncStr '.']);
                    elseif strcmp(FuncStr,'RandomSpinHalf')
                        obj.RandomCfgFunc = FuncHandle;
                        disp(['Starting configuration function changed to ' FuncStr '.']);
                    else
                        disp('Invalid starting configuration function - no change applied.');
                    end
                end
            elseif obj.S == 1
                if isempty(obj.Sector)
                    if strcmp(FuncStr,'RandomSpin1FixedMag')
                        obj.RandomCfgFunc = FuncHandle;
                        disp(['Starting configuration function changed to ' FuncStr '.']);
                    elseif strcmp(FuncStr,'RandomSpin1NoSector')
                        obj.RandomCfgFunc = FuncHandle;
                        disp(['Starting configuration function changed to ' FuncStr '.']);
                    else
                        disp('Invalid starting configuration function - no change applied.');
                    end
                end
            end
        end        
             
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.N = obj.N; Properties.S = obj.S;
            Properties.Sector = obj.Sector;
        end
    end
    
end