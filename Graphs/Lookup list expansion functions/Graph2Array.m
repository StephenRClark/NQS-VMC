% --- Graph connectivity array construction function ---

function [ConArr] = Graph2Array(Graph,Options)
% Takes a stack of provided Graphs and outputs arrays of matrix elements (1
% if sites are linked, 0 otherwise).

% Options are as follows:
% 0 - spin/boson, one species/orbital per site
% 1 - fermion, uu/dd only, no cross-spin pairing terms
% 2 - fermion, ud/du, cross-spin pairing terms
% One option per graph required. If a zero is included in a list of
% ones/twos, it will act as a one.

Ng = Graph(1).N;
if (sum(Options==1) > 0) || (sum(Options==2) > 0)
    N = 2*Ng;
else
    N = Ng;
end
ConArr = zeros(N,N,numel(Graph));
for g = 1:numel(Graph)
    Bonds = Graph(g).Bonds;
    for s = 1:(N/Ng)
    for n = 1:Ng
        for b = 1:size(Bonds,2)
            if Bonds(n,b) ~= 0
                if Options(g) < 2
                    ConArr(n+Ng*(s-1),Bonds(n,b)+Ng*(s-1),g) = 1; 
                    ConArr(Bonds(n,b)+Ng*(s-1),n+Ng*(s-1),g) = 1;
                else
                    ConArr(n+Ng*(s-1),Bonds(n,b)+Ng*(mod(s,2)),g) = 1; 
                    ConArr(Bonds(n,b)+Ng*(s-1),n+Ng*(mod(s,2)),g) = 1;
                end
            end
        end
    end
    end
end