classdef OperatorDg < Operator
    % OperatorDg - an Operator subclass where the operator is diagonal in
    % the configuration basis.
    %   Operator is overarching class.
    
    properties (SetAccess = protected)
        CfgVal % Some function that evaluates the operator value for a particular Cfg.
    end
    
    methods
        function obj = OperatorDg(Hilbert,Graph,CfgVal)
            obj@Operator(Hilbert,Graph);
            obj.CfgVal = CfgVal; % CfgVal should be a function handle.
        end
        
        % LocalSample: found in folder. 
        [CorrSamp] = LocalSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % GraphSample: found in folder.
        [CorrSamp] = GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % CorrMatEls: found in folder.
        [Diff, OpMatEls] = CorrMatEls(obj,Cfg)
                
        % RWMatEls: found in folder.
        [CfgP, RWOpCells, CfgInds] = RWMatEls(obj,~,Cfg,InitCells,CfgInds,Ind);
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Hilbert = obj.Hilbert.PropertyList; 
            Properties.Graph = obj.Graph.PropertyList; 
            Properties.FuncHandle = func2str(obj.CfgVal);
        end
    end
    
end