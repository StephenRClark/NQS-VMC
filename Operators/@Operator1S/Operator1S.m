classdef Operator1S < Operator
    % Operator1S - Operator subclass that deals with one-site correlations
    % that are not diagonal in the configuration basis.
    %   Operator is overarching class.
    
    properties 
        CorrMatEls % Matrix element function of the operator when acting on a configuration.
    end
    
    methods 
        % General constructor for two-site Operator subclass:
        function obj = Operator1S(Hilbert,Graph,CorrMatEls)
            obj@Operator(Hilbert,Graph);
            obj.CorrMatEls = CorrMatEls;            
        end
        
        % Sample the operator over the primary Bonds in the Graph
        % associated with the Operator object.        
        [CorrSamp] = obj.LocalSample(obj,Cfg,EnLoc,dLogp,Ansatz);
        
        % Sample the operator over the entire Graph associated with the
        % Operator object.
        [CorrSamp] = GraphSample(obj,Cfg,EnLoc,dLogp,Ansatz);
    end
end

