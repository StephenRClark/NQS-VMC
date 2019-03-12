classdef Graph
    % Graph - details the connectivity of the sites used.
    %   Graph is overarching class - subclasses will define specific lists
    %   for the properties below.
    
    properties
        N % Number of sites.
        Dim % Vector of number of sites along a nearest neighbour bond direction.
        Bonds % List of principal links between sites.
        LVecs % Principal lattice vectors used in Bonds and BondMap.
        BondMap % Lists of all links between all sites.
    end
    
    methods
        % General constructor for Graph class:
        function obj = Graph(Dim,Bonds,LVecs)
            if nargin < 2 % Default to single site if not given arguments in specific order
                obj.N = 1; obj.Bonds = 0; obj.Dim = 1; obj.LVecs = 0;
            else
                if sum(mod(Dim,1))~=0
                    error('Non-integer site number given.')
                else
                    obj.N = prod(Dim); obj.Dim = Dim; obj.LVecs = LVecs;
                end
                obj.Bonds = Bonds;
            end
            obj.BondMap = GraphMap(obj);
        end
    end
    
end