% --- Graph connectivity mapping function ---

function AllBonds = ReverseBond(Bonds)
% Given a particular lookup table Bonds, this function will return both the
% original lookup table and the lookup table for the reversed bonds in an
% interleaved fashion ( B1 | -B1 | B2 | -B2 | ... ).

N = size(Bonds,1); Nlist = size(Bonds,2);
AllBonds = zeros(N, Nlist);

for l = 1:Nlist
    AllBonds(:,2*l-1) = Bonds(:,l);
    for n = 1:N
        Ind = find(Bonds(:,l)==n);
        if isempty(Ind)
            AllBonds(n,2*l) = 0;
        else
            AllBonds(n,2*l) = Ind;
        end
    end
end

end