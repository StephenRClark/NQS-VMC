function [BondMapR] = RotateBonds(GraphObj)

N = size(GraphObj.BondMap{1},1); Dimension = numel(GraphObj.Dim);
Nb = numel(GraphObj.BondMap); L = GraphObj.Dim(1);
if Dimension == 2 && (GraphObj.Dim(1) == GraphObj.Dim(2))
    BondMapR = cell(1,4*Nb);
    for b = 1:Nb
        BondMapR{b} = GraphObj.BondMap{b};
        BondMapR{b+Nb} = zeros(N,1);
        BondMapR{b+2*Nb} = zeros(N,1);
        BondMapR{b+3*Nb} = zeros(N,1);
        for n = 1:N
            n0 = GraphObj.BondMap{b}(n);
            i = 1 + mod(n0-1,L); j = ceil(n0/L);
            n90 = 1 + mod(-j,L) + (i-1)*L;
            n180 = 1 + mod(-i,L) + mod(-j,L)*L;
            n270 = j + mod(-i,L)*L;
            BondMapR{b+Nb}(n) = n90;
            BondMapR{b+2*Nb}(n) = n180;
            BondMapR{b+3*Nb}(n) = n270;
        end
    end
end