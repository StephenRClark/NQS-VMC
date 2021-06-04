% --- Bipartite hypercube Graph object setup function --- %

function [AGraph, BGraph] = HypCubBP(Dim,Bound,LVecs,SFlag)
if mod(prod(Dim),2) == 1
    error('Cannot create bipartite Graph from lattice with odd number of sites.')
end
% Dim - Input vector of hypercube lengths in each desired dimension.
% LVecs - Principal lattice vectors to use as Bonds for Graph.
if nargin < 4
    SFlag = 1; % Default to lattice spanning map if no input.
elseif nargin == 2
    % Assume principal vectors of lattice are being used.
    LVecs = eye(numel(Dim)); SFlag = 1;
elseif nargin == 0
    Dim = 1; Bound = 'PBC'; LVecs = 0; SFlag = 0;
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
SLInds = zeros(N,size(LVecs,1));
for n = 1:N
    nvec = zeros(numel(Dim),1);
    for d = 1:numel(Dim)
        nvec(d) = 1 + mod(ceil(n/prod(Dim(1:(d-1)))) -1,Dim(d));
    end
    SLInds(n) = (1+mod(sum(nvec)+numel(Dim),2));
end
ABonds = Bonds; BBonds = Bonds;
ABonds(SLInds==2,:) = 0; BBonds(SLInds==1,:) = 0;
AGraph = Graph(Dim,ABonds,LVecs,SFlag);
BGraph = Graph(Dim,BBonds,LVecs,SFlag);
end