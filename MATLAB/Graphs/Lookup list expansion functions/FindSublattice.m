function [GraphObj] = FindSublattice(GraphObj)
% Find any applicable sublattices in the provided Graph object, and label
% the sites with their sublattice indices in GraphObj.SLInds;
BondMap = GraphObj.BondMap; Nmap = numel(BondMap);
N = GraphObj.N; KOList = 1:N; SLInds = zeros(N,1);

start = 0; ind = 0;
while (sum(KOList)>0)
    start = start+1;
    if KOList(start)>0
        ind = ind+1;
        for m = 1:Nmap
            site = BondMap{m}(start);
            if KOList(site) > 0
                KOList(site) = 0;
                SLInds(site) = ind;
            else
                break
            end
        end
    end
end
GraphObj.SLInds = SLInds;
end