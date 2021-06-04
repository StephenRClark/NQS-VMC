classdef OperatorH
    % OperatorH - Operator variant that evaluates combinations of
    % subordinate Operator objects in order(s) specified by OpOrder,
    % including the Hamiltonian H.
    %   OperatorH takes Operator objects and shares some methods but is not
    %   an Operator subclass.
    
    properties
        Hilbert % Details Hilbert space in which Ansatz is represented.
        SubOperator % Subordinate Operator that is to be evaluated.
        OpOrder % Cell array with three entries per operator detailing composition order.
        Hamiltonian % Relevant Hamiltonian object.
        HOrder % {2 x 1} cell array with entries detailing Hamiltonian position in
               % left and right operator sampling.
    end
    
    methods
        % General constructor for composite Operator subclass:
        function obj = OperatorH(Hilbert,Hamiltonian,SubOperator,OpOrder,HOrder)
            obj.Hilbert = Hilbert; obj.SubOperator = SubOperator; obj.OpOrder = OpOrder;
            obj.Hamiltonian = Hamiltonian; obj.HOrder = HOrder;
        end
        
        % GraphSample: sample the operator over the entire Graph associated
        % with the Operator object. Note that no LocalSample option is
        % available, as composite Operators are only designed for
        % outputting end observables, not for use in Hamiltonian objects.
        [CorrSamp] = obj.GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
    end
end

