classdef HypCub < Graph
    % Hypercube - an N-dimensional hypercubic graph.
    %   HypCub is a subclass of Graph, with specific bond construction.
    
    properties
        Bound % Specifies the type of boundary.
    end
    
    methods
        function obj = HypCub(Dim,Bound,LVecs)
            % Dim - Input vector of hypercube lengths in each desired dimension.
            % LVecs - Principal lattice vectors to use as Bonds for Graph.
            if nargin > 3
                error('Too many arguments - specify all hypercube lengths in a single vector for input.')
            elseif nargin == 2
                % Assume principal vectors of lattice are being used.
                LVecs = eye(numel(Dim));
            elseif nargin == 0
                Dim = 1; Bound = 'PBC'; LVecs = 0;
            end
            if strcmp(Bound,'PBC') == 0 && strcmp(Bound,'OBC') == 0
                error('Ensure "Bound" entry is either "OBC" or "PBC" string.')
            end
            % Determine number of sites from Dim.
            N = prod(Dim);
            Bonds = ones(N,size(LVecs,1)); % Bonds are listed by increasing site index.
            for n = 1:N
                Bonds(n,:) = Bonds(n,:)*n;
                for l = 1:size(LVecs,1)
                    for d = 1:numel(Dim)
                        if strcmp(Bound,'PBC') % Periodic boundary conditions.
                            Bonds(n,l) = 1 + mod( Bonds(n,l) - 1 + ...
                                (LVecs(l,d) * prod(Dim(1:(d-1)))) , prod(Dim(1:d)) ) + ...
                                (ceil( Bonds(n,l) / prod(Dim(1:d)) ) - 1) * prod(Dim(1:d));
                        elseif strcmp(Bound,'OBC') % Open boundary conditions.
                            Bonds(n,l) = ( 1 + mod( Bonds(n,l) - 1 + ...
                                (LVecs(l,d) * prod(Dim(1:(d-1)))) , prod(Dim(1:d)) ) + ...
                                (ceil( Bonds(n,l) / prod(Dim(1:d)) ) - 1) * prod(Dim(1:d)) ) * ...
                                (mod(ceil(Bonds(n,l)/prod(Dim(1:(d-1)))),Dim(d))>0);
                            % A zero entry in Bonds is a non-existent link.
                        end
                    end
                end
            end
            
            % Initialise using Graph superclass.
            obj@Graph(Dim,Bonds,LVecs);
            
            obj.Bound = Bound;
        end
    end
    
end

