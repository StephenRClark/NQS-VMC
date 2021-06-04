classdef Operator
    % Operator - acts within a Hilbert space of configuration vectors and
    % links a given configuration to those related by an operator.
    %   Operator is overarching class - subclasses will define specific
    %   functions for the methods listed below.
    
    properties 
        Hilbert % Details Hilbert space in which Ansatz is represented.
        Graph % Details nearest neighbour connectivity of sites.
    end
    
    methods 
        % General constructor for class - just assigns input Hilbert and Graph.
        function obj = Operator(Hilbert,Graph)
           obj.Hilbert = Hilbert; obj.Graph = Graph;
        end        
    end
    
    methods (Abstract)
        % LocalSample: Sample the operator over the primary Bonds in the
        % Graph associated with the Operator object.
        [CorrSamp] = LocalSample(obj,Cfg,EnLoc,dLogp,AnsatzObj);
        
        % GraphSample: Sample the operator over the entire Graph associated
        % with the Operator object.
        [CorrSamp] = GraphSample(obj,Cfg,EnLoc,dLogp,AnsatzObj);
        
        % CorrMatEls: Given a Cfg struct, output the matrix elements and
        % configuration differences linked by this Operator.
        [Diff, OpMatEls] = CorrMatEls(obj,Cfg);
    end
    
end

