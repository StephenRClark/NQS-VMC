% --- Graph connectivity mapping function ---

function GraphObj = GraphMap(GraphObj)
% Given a particular Graph object with each site's principal lattice vector
% lookup lists in Bonds, will map out all possible two-site pairs linked by
% some number of these vectors. Automatically accounts for boundary
% conditions.

Bonds = GraphObj.Bonds; Dim = GraphObj.Dim; N = prod(Dim); LVecs = GraphObj.LVecs;
% Need to determine number of translates by LVecs entries until one
% performs a full lattice translation (i.e. n --> m * LVecs --> n).
Nvec = size(LVecs,1); NtrT = 1; IndVec = zeros(Nvec,1); MaxTr = zeros(Nvec,1);
for v = 1:size(LVecs,1)
    w = 1; IndVec(v) = NtrT;
    while sum(abs(rem(w*LVecs(v,:),Dim))) > 0
        w = w + 1;
    end
    MaxTr(v) = w;
    NtrT = NtrT * w;
end
GraphMapT = cell(NtrT,1); VecIndsT = zeros(NtrT,Nvec); N0List = zeros(NtrT);

for n = 1:NtrT
    GraphMapT{n} = zeros(N,1);
    % Entries are indexed by starting site, so the value of GraphMap{n}(m)
    % is the site linked with site m by the nth dR vector (the form of
    % which will depend on the graph).
end

% Dim here is taken to be the number of nearest neighbours along a bond
% direction within the lattice being defined - assumes all bonds of the
% same direction are in the same column of Bonds.

KOList = ones(N,1) * (1:N); % Knockout list to ensure that only distinct mappings remain.
NtrA = 0;
for t = 1:NtrT
    for n = 1:N
        Dest = n;
        for v = 1:Nvec
            VecIndsT(t,v) = mod(floor((t-1)/IndVec(v)),MaxTr(v));
            for x = 1:VecIndsT(t,v)
                if Dest ~= 0
                    Dest = Bonds(Dest,v);
                    if Dest == 0
                        break;
                    end
                end
            end
        end
        GraphMapT{t}(n) = Dest;
        if Dest ~= 0
            if KOList(n,Dest) ~= 0
                KOList(n,Dest) = 0;
                N0List(t) = 1;
            end
        end
    end
    NtrA = NtrA + N0List(t);
end

GraphMap = cell(NtrA,1); VecInds = zeros(NtrA,Nvec); p = 1;
for n = 1:NtrT
    if N0List(n) > 0
        GraphMap{p} = GraphMapT{n}; VecInds(p,:) = VecIndsT(n,:); p = p + 1;        
    end
end

GraphObj.Ntr = NtrA;
GraphObj.BondMap = GraphMap;
GraphObj.VecInds = VecInds;

end