% --- Fermionic sign calculation function ---

function [Diff] = FermSignCfg(Diff,Cfg)
% This function adds fermionic anticommutator sign information to the Diff
% struct using Cfg.up and Cfg.dn.

% Assumption is that Diff will represent either no moves, single fermion
% movements or fermion site swaps.

if Diff.num == 0 % No change.
    Diff.sign = 1;
elseif Diff.num == 1 % Single fermion change.
    if Diff.type == 1 % Up fermion is moving.
        Diff.sign = (-1)^sum((Cfg.up<min(Diff.pos)) .* (Cfg.up<max(Diff.pos)));
        % Sign depends on number of up fermions between start and end.
    elseif Diff.type == -1
        Diff.sign = (-1)^sum((Cfg.dn<min(Diff.pos)) .* (Cfg.dn<max(Diff.pos)));
        % Sign depends on number of down fermions between start and end.
    end
elseif Diff.num == 2 % Fermion pair location swaps or pair moves.
    if Diff.type == 0 % Pair location swap.
        Diff.sign = (-1) * (-1)^(sum(Cfg.up<Diff.pos(1)) + sum(Cfg.up<Diff.pos(2))) * ...
            (-1)^(sum(Cfg.dn<Diff.pos(2)) + sum(Cfg.dn<Diff.pos(1)));
        % Extra sign arises from index ordering of a+(m)a-(n)b+(n)b-(m).
    elseif Diff.type == 2 % Pair move.
        Diff.sign = Diff.sign * (-1)^sum((Cfg.up<min(Diff.pos)) .* (Cfg.up<max(Diff.pos)))...
            * (-1)^sum((Cfg.dn<min(Diff.pos)) .* (Cfg.dn<max(Diff.pos)));
        % Sign change is product of sign changes from individual moves
        % (assuming C*{up}C{up}C*{dn}C{dn} ordering). Any sign changes from
        % rearranging to this form should be pre-assigned.
    end
else
    Diff.sign = 1;
end