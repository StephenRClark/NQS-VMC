classdef OperatorDg < Operator
    % OperatorDg - an Operator subclass where the operator is diagonal in
    % the configuration basis.
    %   Operator is overarching class.
    
    properties
        CfgVal % Some function that evaluates the operator value for a particular Cfg.
    end
    
    methods
        function obj = OperatorDg(Hilbert,Graph,CfgVal)
            obj@Operator(Hilbert,Graph);
            obj.CfgVal = CfgVal; % CfgVal should be a function handle.
        end
        
        % Sample the operator over the primary Bonds in the Graph
        % associated with the Operator object.        
        [CorrSamp] = LocalSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % Sample the operator over the entire Graph associated with the
        % Operator object.
        [CorrSamp] = GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % Generate a zero difference and the configuration value when
        % CorrMatEls is requested.
        [Diff,CorrMatEls] = CorrMatEls(obj,Hilbert,Cfg,Graph);
    end
end

