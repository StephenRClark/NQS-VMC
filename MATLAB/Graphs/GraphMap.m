% --- Graph connectivity mapping function ---

function GraphMap = GraphMap(Graph)
% Given a particular Graph object with each site's principal lattice vector
% lookup lists in Bonds, will map out all possible two-site pairs linked by
% some number of these vectors. Automatically accounts for boundary
% conditions.

Bonds = Graph.Bonds; Dim = Graph.Dim; N = prod(Dim); LVecs = Graph.LVecs;

% Need to determine number of translates by LVecs entries until one
% performs a full lattice translation (i.e. n --> m * LVecs --> n).
Nvec = zeros(1,size(LVecs,1));
for l = 1:size(LVecs,1)
    w = 1;
    while sum(abs(rem(w*LVecs(l,:),Dim))) > 0
        w = w + 1;
    end
    Nvec(l) = w;
end
if numel(Nvec) > 1
    GraphMapT = cell(Nvec); Nsep = prod(Nvec); N0List = zeros(Nvec);
else
    GraphMapT = cell(Nvec,1); Nsep = Nvec; N0List = zeros(Nvec,1);
end

for n = 1:Nsep
    GraphMapT{n} = zeros(N,1);
    % Entries are indexed by starting site, so the value of GraphMap{n}(m)
    % is the site linked with site m by the nth dR vector (the form of
    % which will depend on the graph).
end

% Dim here is taken to be the number of nearest neighbours along a bond
% direction within the lattice being defined - assumes all bonds of the
% same direction are in the same column of Bonds.

KOList = ones(N,1) * (1:N); % Knockout list to ensure that only distinct mappings remain.

for v = 1:numel(Nvec)
    for dx = 0:1:(Nvec(v)-1)
        % Bond map is generated iteratively - first, site pairs for all
        % dR(1) are calculated, then from those, all dR(2) and so on.
        for x = 1:prod(Nvec(1:(v-1)))
            if (x*v + dx) == 1
                % Represents dR = zero vector case.
                GraphMapT{1} = (1:N)';
            end
            for n = 1:N
                Dest = GraphMapT{x}(n);
                if Dest ~= 0
                    for m = 1:dx
                        Dest = Bonds(Dest,v);
                        if Dest == 0
                            break
                        end
                    end
                end
                if Dest ~= 0
                    if KOList(n,Dest) ~= 0 % Check knockout list.
                        GraphMapT{x+dx*prod(Nvec(1:(v-1)))}(n) = Dest;
                        KOList(n,Dest) = 0;
                    end
                end
            end
            N0List(x+dx*prod(Nvec(1:(v-1)))) = sum(GraphMapT{x+dx*prod(Nvec(1:(v-1)))});
        end
    end
end

GraphMap = cell(1,sum(N0List(:)>0)); p = 1;
for n = 1:numel(N0List)
    if N0List(n) > 0
        GraphMap{p} = GraphMapT{n}; p = p + 1;
    end
end

end