classdef Operator2S < Operator
    % Operator2S - Operator subclass that deals with two-site correlations.
    %   Operator is overarching class.
    
    properties (SetAccess = protected)
        CorrMatElsFunc % Matrix element function of the operator when acting on a configuration.
    end
    
    methods
        % General constructor for two-site Operator subclass:
        function obj = Operator2S(Hilbert,Graph,CorrMatEls)
            obj@Operator(Hilbert,Graph);
            obj.CorrMatElsFunc = CorrMatEls;
        end
        
        % LocalSample: found in folder. 
        [CorrSamp] = LocalSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % GraphSample: found in folder.
        [CorrSamp] = GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % CorrMatEls: simply call the associated CorrMatElsFunc.
        function [Diff, OpMatEls] = CorrMatEls(obj,Cfg)
            [Diff, OpMatEls] = obj.CorrMatElsFunc(obj.Hilbert,Cfg,obj.Graph);
        end
        
        % RWMatEls: found in folder.
        [CfgP, RWOpCells, CfgInds] = RWMatEls(obj,~,Cfg,InitCells,CfgInds,Ind);
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Hilbert = obj.Hilbert.PropertyList; 
            Properties.Graph = obj.Graph.PropertyList; 
            Properties.FuncHandle = func2str(obj.CorrMatElsFunc);
        end
    end
    
end