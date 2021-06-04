% --- 2D square lattice plaquette labelling function ---

function [PlaqLabels] = SquarePlaqLabel(GraphObj)
% Takes a provided Graph object and outputs a N_plaq x N_vertex list of
% plaquette indices, such that PlaqLabels(i,j) is the jth vertex of
% plaquette i.

% Should use as GraphObj.ExtraLabels{a} = SquarePlaqLabel(GraphObj) to
% store the plaquette label list.

% Will probably have unexpected interactions with non-2d non
% nearest-neighbour HypCub objects.

Dim = GraphObj.Dim; N = GraphObj.N; Bound = GraphObj.Bound; Bonds = GraphObj.Bonds;

% Determine number of square plaquettes from BCs.
PDim = Dim + Bound - 1; Nplaq = prod(PDim); PlaqLabels = zeros(Nplaq,4);
% Loop through all sites and if eligible (no unbroken bonds), list.
p = 0;
for n = 1:N
    List = [n Bonds(n,1) Bonds(n,2)];
    if Bonds(n,1)~=0
        List = [List Bonds(Bonds(n,1),2)];
        if sum(List==0) == 0
            p = p + 1;
            PlaqLabels(p,:) = List;
        end
    end
end
end