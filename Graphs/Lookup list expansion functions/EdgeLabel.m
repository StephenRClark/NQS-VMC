% --- Graph edge labelling function ---

function [EdgeLabels] = EdgeLabel(GraphObj)
% Takes a provided Graph object and outputs a N_edge x 2 list of
% edge indices, such that EdgeLabels(i,j) is the jth vertex of
% edge i.

% Should use as GraphObj.ExtraLabels{a} = EdgeLabel(GraphObj) to store the
% edge label list.

N = GraphObj.N; Bonds = GraphObj.Bonds;

Nedge = sum(Bonds(:)~=0); EdgeLabels = zeros(Nedge,2);
% Loop through all sites and if eligible (no unbroken bonds), list.
e = 0;
for n = 1:N
    for b = 1:size(Bonds,2)
        List = [n Bonds(n,b)];
        if sum(List==0) == 0
            e = e + 1;
            EdgeLabels(e,:) = List;
        end
    end
end