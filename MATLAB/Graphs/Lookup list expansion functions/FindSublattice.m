function [GraphObj] = FindSublattice(GraphObj)
% Find any applicable sublattices in the provided Graph object, and label
% the sites with their sublattice indices in GraphObj.SLInds;
BondMap = GraphObj.BondMap; Nmap = numel(BondMap);
N = GraphObj.N; KOList = 1:N; SLInds = zeros(N,1); ind = 0;

for n = 1:N
    start = n;
    if SLInds(start) == 0
        ind = ind+1;
        sub_ind = ind;
    else
        sub_ind = SLInds(start);
    end
    sitelist = [];
    for m = 1:Nmap
        site = BondMap{m}(start);
        if site ~= 0
            sitelist = [sitelist; site];
        end
    end
    for s = 1:numel(sitelist)
        if SLInds(sitelist(s))~=0
            sub_ind = SLInds(sitelist(s));
            break
        end
    end
    SLInds(sitelist) = sub_ind;
end
GraphObj.SLInds = SLInds;
end