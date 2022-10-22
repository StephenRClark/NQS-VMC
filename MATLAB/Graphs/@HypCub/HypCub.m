classdef HypCub < Graph
    % Hypercube - an N-dimensional hypercubic graph.
    %   HypCub is a subclass of Graph, with specific bond construction.
    
    methods
        function obj = HypCub(Dim,Bound,LVecs,SFlag)
            % Dim - Input vector of hypercube lengths in each desired dimension.
            % LVecs - Principal lattice vectors to use as Bonds for Graph.
            % Determine number of sites from Dim.
            N = 1; DimIndVec = ones(1,size(LVecs,1));
            for d = 1:numel(Dim)
                DimIndVec(d) = N;
                N = N * Dim(d);
            end
            Bonds = ones(N,size(LVecs,1)); % Bonds are listed by increasing site index.
            CoOrds = zeros(N,numel(Dim));
            for n = 1:N
                CoOrds(n,:) = mod(floor((n-1)./DimIndVec),Dim);
            end
            for b = 1:size(LVecs,1)
                CoOrdsT = CoOrds + (ones(N,1) * LVecs(b,:)); 
                Bonds(:,b) = (1 + sum(mod(CoOrdsT,Dim).*DimIndVec,2)) .* ...
                    abs(prod(Bound.^((CoOrdsT>=Dim)+(CoOrdsT<0)),2));
            end
            % Initialise using Graph superclass.
            obj@Graph(Dim,Bonds,Bound,LVecs,SFlag);
            % Set Graph sublattice indices.
            obj = FindSublattice(obj);
        end
        
        function [obj] = AddRotations(obj)
            RotatedBonds = RotateBonds(obj);
            obj.BondMap = RotatedBonds;
        end
    end
    
end