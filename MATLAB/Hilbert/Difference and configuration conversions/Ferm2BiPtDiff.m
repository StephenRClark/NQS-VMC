% --- Difference structure conversion function ---

function [Diff] = Ferm2BiPtDiff(Diff,N)
% This function converts a Diff struct compatible with fermionic moves to a
% Diff struct compatible with a bipartite Modifier.

% Assumption is that Diff will represent either no moves, single fermion
% movements or fermion site swaps.

if size(Diff.val,1) == 1 % Single fermion move.
    if Diff.type == 1 % Up fermion move - no adjustment to Diff.pos required
        Diff.num = 2; % Two local occupations modified by a fermion moving in BiPt scheme.
    elseif Diff.type == -1 % Down fermion move.
        Diff.pos = Diff.pos + N; % Move takes place in later half of expanded lattice.
        Diff.num = 2; % Two local occupations modified by a fermion moving in BiPt scheme.
    end
elseif size(Diff.val,1) == 2 % Two fermions move.
    Diff.val = [Diff.val(1,:), Diff.val(2,:)];
    Diff.num = 4; % Four local occupations modified by two fermions moving in BiPt scheme.
    Diff.pos = [Diff.pos, Diff.pos + N];
end