classdef Graph
    % Graph - details the connectivity of the sites used.
    %   Graph is overarching class - subclasses will define specific lists
    %   for the properties below.
    
    properties
        Bonds % List of principal links between sites.
        SLInds = 1; % Sublattice indices - used for certain Operators.
        ExtraLabels = {}; % Cell list for extra labelling e.g. plaquettes.
    end
    
    properties (SetAccess = protected)
        N % Number of sites.
        Dim % Vector of number of sites along a nearest neighbour bond direction.
        Bound = 1; % Vector of factors accumulated when crossing boundaries in each direction.
        % 1 - periodic boundary, 0 - closed boundary, -1 - antiperiodic boundary (WIP).
        LVecs % Principal lattice vectors used in Bonds and BondMap.
        BondMap % Lists of all links between all sites.
        Ntr = 1; % Number of translates in BondMap.
        % Can be the same as Bonds if Graph is specified as not spanning.
        VecInds = 1;% Lists of number of lattice vector translations per BondMap entry.
    end
    
    methods
        % General constructor for Graph class:
        function obj = Graph(Dim,Bonds,Bound,LVecs,SFlag)
            if numel(Bound) ~= numel(Dim)
                error('Number of boundary conditions does not match number of dimensions.');
            end
            if sum(mod(Dim,1))~=0
                error('Non-integer site number given.')
            else
                obj.N = prod(Dim); obj.Dim = Dim; obj.LVecs = LVecs;
            end
            obj.Bonds = Bonds;
            for b = 1:numel(Bound)
                if abs(Bound) > 1
                    error('Bound entries must be set to 1 (PBC), 0 (closed) or -1 (antiperiodic).');
                end
            end
            obj.Bound = Bound;
            if SFlag == 1
                obj = GraphMap(obj);
            else
                obj.BondMap = cell(size(Bonds,2),1);
                obj.VecInds = eye(size(Bonds,2));
                for b = 1:size(Bonds,2)
                    obj.BondMap{b} = Bonds(:,b);
                end
                obj.Ntr = size(Bonds,2);
            end
            obj.SLInds = ones(obj.N,1);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Dim = obj.Dim; Properties.Bound = obj.Bound;
            Properties.LVecs = obj.LVecs; Properties.Bonds = obj.Bonds;
            Properties.SLInds = obj.SLInds; Properties.Ntr = obj.Ntr;
            Properties.ExtraLabels = obj.ExtraLabels;
        end
    end
    
end