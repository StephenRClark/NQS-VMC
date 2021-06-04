classdef OperatorMS < Operator
    % OperatorMS - Operator subclass that deals with multi-site 
    % correlations, requiring the site groups to be listed in
    % Graph.ExtraLabels.
    %   Operator is overarching class.
    
    properties (SetAccess = protected)
        CorrMatElsFunc % Matrix element function of the operator when acting on a configuration.
        ListIndex % Index detailing which list in Graph.ExtraLabels to use.
    end
    
    methods
        % General constructor for two-site Operator subclass:
        function obj = OperatorMS(Hilbert,Graph,CorrMatEls,ListIndex)
            obj@Operator(Hilbert,Graph);
            obj.CorrMatElsFunc = CorrMatEls;
            if (ListIndex <= 0 || ListIndex > numel(Graph.ExtraLabels))
                error(['Invalid list index - should be an integer between 1 and ' numel(Graph.ExtraLabels) '.']);
            end
            obj.ListIndex = ListIndex;
        end
        
        % LocalSample: found in folder. 
        [CorrSamp] = LocalSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % GraphSample: found in folder.
        [CorrSamp] = GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % CorrMatEls: simply call the associated CorrMatElsFunc.
        function [Diff, OpMatEls] = CorrMatEls(obj,Cfg)
            [Diff, OpMatEls] = obj.CorrMatElsFunc(obj.Hilbert,Cfg,obj.Graph,obj.ListIndex);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Hilbert = obj.Hilbert.PropertyList; 
            Properties.Graph = obj.Graph.PropertyList; 
            Properties.FuncHandle = func2str(obj.CorrMatElsFunc);
        end
    end
    
end