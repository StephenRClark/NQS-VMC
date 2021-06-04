classdef OperatorC 
    % OperatorC - Operator variant that evaluates combinations of
    % subordinate Operator objects in order(s) specified by OpOrder.
    %   OperatorC takes Operator objects and shares some methods but is not
    %   an Operator subclass.
    
    properties
        Hilbert % Details Hilbert space in which Ansatz is represented.
        Graph % Primary graph associated with the composite Operator.
        SubOperator % Subordinate Operator that is to be evaluated with the Hamiltonian.
        OpOrder % Cell array with three entries per operator detailing composition order.
    end
    
    methods
        % General constructor for composite Operator subclass:
        function obj = OperatorC(Hilbert,SubOperator,OpOrder)
            obj.Hilbert = Hilbert; obj.SubOperator = SubOperator; obj.OpOrder = OpOrder;
        end
        
        % GraphSample: sample the operator over the entire Graph associated
        % with the Operator object. Note that no LocalSample option is
        % available, as composite Operators are only designed for
        % outputting end observables, not for use in Hamiltonian objects.
        [CorrSamp] = obj.GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
    end
end

