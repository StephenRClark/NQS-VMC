classdef Hamiltonian
    % Hamiltonian - an object containing the details of the Hamiltonian
    % being applied to the physical system.
    %   Hamiltonian objects are fed into Sampler objects for sampling with
    %   an Ansatz.
    
    properties
        Operator % List of Operator objects that compose the Hamiltonian.
        HParams % List of parameters that characterise the Hamiltonian.
        TFuncs % Cell list of time-dependent functions that modify the HParams. Can be left empty.
        % Note that each Operator needs its own Graph that specifies their
        % connectivity on the lattice.
    end
    
    properties (Hidden)
        HParams0 % Original input Hamiltonian parameters - used as a reference when time evolving.
    end
    
    methods
        % Constructor for the Hamiltonian object:
        function [obj] = Hamiltonian(Operators,HParams,TFuncs)
            if nargin == 2 % if no TFuncs specified, assume no time evolution functions.            
                TFuncs = cell(numel(Operators),1);                
            end
            obj.TFuncs = TFuncs;
            if numel(obj.TFuncs) ~= numel(Operators)
                error(['Number of time evolution functions does not match number of parameters. '...
                    '3rd input should be a cell list with an equal number of entries (empty cells are fine)']);
            end
            if numel(HParams) ~= numel(Operators)
                error('Number of parameters does not match number of operators.');
            end
            obj.Operator = Operators; obj.HParams = HParams;      
        end
        
        % HamMatEls: generate the matrix elements and the differences
        % created by the action of the Hamiltonian on a configuration Cfg.
        function [Diff,HMatEls] = HamMatEls(obj,Cfg)
            for h = 1:numel(obj.HParams)
                if h == 1
                    [Diff,OpMatEls] = obj.Operator{h}.CorrMatEls(Cfg);
                    HMatEls = OpMatEls * obj.HParams(h); Diff = reshape(Diff,numel(Diff),1);
                else
                    [DiffT,OpMatElsT] = obj.Operator{h}.CorrMatEls(Cfg);
                    HMatEls = [HMatEls; OpMatElsT*obj.HParams(h)];
                    DiffT = reshape(DiffT,numel(DiffT),1); Diff = [Diff; DiffT]; 
                end
            end
        end
        
        % TimeEvolveH: Modify HParams according to any TFuncs specified at
        % initialisation. Only called during time-dependent optimisations.
        function obj = TimeEvolveH(obj,time)
            for p = 1:numel(obj.HParams)
                if isempty(obj.TFuncs{p}) == 0
                    obj.HParams(p) = obj.TFuncs{p}(obj.HParams0(p),time);
                end
            end
        end
        
        % EnLoc: given a Cfg struct, sample the individual operators
        % composing the Hamiltonian and output the local energy EnLoc.
        function [EnLoc] = EnergySample(obj,Cfg,Ansatz)
            EnLoc = 0;
            for h = 1:numel(obj.HParams)
                % Call in turn each Operator's LocalSample function over
                % their associated Graph's Bonds, then weight by J
                % contribution.
                EnLoc = EnLoc + obj.HParams(h)*obj.Operator{h}.LocalSample(Cfg,0,0,Ansatz);
            end
        end
        
        % RWMatEls: Given a Cfg struct, output the matrix elements weighted
        % by PsiRatios and the configurations linked by this Operator.
        % Found in folder.
        [CfgP, RWOpCells, CfgIndsP] = RWMatEls(obj,Ansatz,Cfg,InitCells,CfgInds);
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.HParams = obj.HParams;
            Properties.Operator = cell(numel(obj.Operator),1);
            for o = 1:numel(obj.Operator)
                Properties.Operator{o} = obj.Operator{o}.PropertyList;
            end            
        end
    end
    
end