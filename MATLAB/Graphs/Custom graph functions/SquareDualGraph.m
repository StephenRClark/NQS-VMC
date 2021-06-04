% --- Dual Graph generating function ---

function [DualObj] = SquareDualGraph(GraphObj)
% This function takes an existing Graph object and attempts to create a
% Graph representing the dual of the Graph.

% Current version assumes a 2D square lattice.

Dim = GraphObj.Dim; Bound = GraphObj.Bound; Bonds = GraphObj.Bonds;

N = GraphObj.N; PDim = Dim + Bound - 1;

% Dual lattice site labelling takes all sites linked in one direction as
% 1,...,Ns1, then all sites linked in the other direction as
% Ns1+1,...,Ns1+Ns2.

% Generate Bonds by temporarily creating two HypCub Graphs with the same
% boundaries and altered lengths.

LatA = HypCub([PDim(1),Dim(2)],Bound,eye(2),0); LatB = HypCub([Dim(1),PDim(2)],Bound,eye(2),0);
Ns1 = LatA.N; Ns2 = LatB.N; NsT = Ns1 + Ns2;
DualBonds = zeros(NsT,2); DualBonds(1:Ns1,:) = LatA.Bonds; 
DualBonds((1:Ns2)+Ns1,:) = LatB.Bonds + Ns1*(LatB.Bonds~=0);

DualObj = Graph(NsT,DualBonds,0,0,0); % Generate custom connectionless Graph.

% Store plaquette and vertex site lists in ExtraLabels.

% Convention - first set is vertices, second set is plaquettes.

Vertices = zeros(N,4); % Each original site is a vertex.
Plaquettes = zeros(prod(PDim),4); % Number of plaquettes depends on BCs.

for n = 1:N
    Vertices(n,1) = (n - (1-Bound(1))*(ceil(n/Dim(1))-1))*((Bound(1)+mod(n,Dim(1)))>0);
    Vertices(n,2) = max([0,find(Bonds(:,1)==n)]) - (1-Bound(1))*(ceil(n/Dim(1))-1)*sum(Bonds(:,1)==n); 
    Vertices(n,3) = (Ns1 + n)*((Bound(2)+mod(ceil(n/Dim(1)),Dim(2)))>0); 
    Vertices(n,4) = max([0,find(Bonds(:,2)==n)])+(Ns1*sum(Bonds(:,2)==n)); 
end

for x = 1:PDim(1)
    for y = 1:PDim(2)
        SA = x + (y-1)*PDim(1); SB = x + (y-1)*Dim(1);
        Plaquettes(SA,:) = [SA, LatA.Bonds(SA,2), SB+Ns1, LatB.Bonds(SB,1)+Ns1];
    end
end

DualObj.SLInds = [ones(Ns1,1);2*ones(Ns2,1)];
DualObj.ExtraLabels{1} = Vertices; DualObj.ExtraLabels{2} = Plaquettes;
end