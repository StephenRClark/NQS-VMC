% --- Difference to configuration conversion function ---

function [CfgP] = Diff2CfgFerm(Diff,Cfg)
% This function converts a Diff struct compatible with fermion moves to a
% CfgP struct when provided a starting Cfg.

% Assumption is that Diff will represent either no moves, single site
% modifications or two-site alterations. If multiple Diffs are passed,
% assume the Diffs are in a Nd x 1 struct array.

Nd = numel(Diff); % Function should be able to handle multiple Diffs being passed.

% Initialise structure array CfgP with last element first.
CfgP(Nd) = Cfg; CfgP = CfgP.'; % Transpose to output Nd x 1 struct array.

for n = 1:Nd
    CfgP(n) = Cfg;
    % Alter each CfgP using associated Diff. In fermion case, Diff.val
    % records local density changes. Moves are classified in Diff.type.
    for d = 1:Diff(n).num % Diff.num refers to number of moving particles.
        if Diff(n).type == 1 % Up fermion move.
            % Diff.pos for single moves are formatted as [Start Dest]
            % normally, but will use Diff.val to check.
            Start = Diff(n).pos(Diff(n).val<0);
            Dest = Diff(n).pos(Diff(n).val>0);
            CfgP(n).up(CfgP(n).up==Start) = [];
            CfgP(n).up = [CfgP(n).up(CfgP(n).up<Dest) ...
                Dest CfgP(n).up(CfgP(n).up>Dest)];
        elseif Diff(n).type == -1 % Down fermion move.
            Start = Diff(n).pos(Diff(n).val<0);
            Dest = Diff(n).pos(Diff(n).val>0);
            CfgP(n).dn(CfgP(n).dn==Start) = [];
            CfgP(n).dn = [CfgP(n).dn(CfgP(n).dn<Dest) ...
                Dest CfgP(n).up(CfgP(n).up>Dest)];
        elseif Diff(n).type == 0 % Fermion pair swap.
            % Two fermions move - convention is to handle up fermion first.
            % First listed site will be initial up fermion site.
            Start = Diff(n).pos(1+mod(d-1,2));
            Dest = Diff(n).pos(1+mod(d,2));
            if mod(d,2) == 1
                CfgP(n).up(CfgP(n).up==Start) = [];
                CfgP(n).up = [CfgP(n).up(CfgP(n).up<Dest) ...
                    Dest CfgP(n).up(CfgP(n).up>Dest)];
            else
                CfgP(n).dn(CfgP(n).dn==Start) = [];
                CfgP(n).dn = [CfgP(n).dn(CfgP(n).dn<Dest) ...
                    Dest CfgP(n).up(CfgP(n).up>Dest)];
            end
        elseif Diff(n).type == 2 % Fermion pair move.
            % Two fermions move - convention is to handle up fermion first.
            % First listed site will be initial doublon site.
            Start = Diff(n).pos(2*(ceil(d/2))-1);
            Dest = Diff(n).pos(2*(ceil(d/2)));
            if mod(d,2) == 1
                CfgP(n).up(CfgP(n).up==Start) = [];
                CfgP(n).up = [CfgP(n).up(CfgP(n).up<Dest) ...
                    Dest CfgP(n).up(CfgP(n).up>Dest)];
            else
                CfgP(n).dn(CfgP(n).dn==Start) = [];
                CfgP(n).dn = [CfgP(n).dn(CfgP(n).dn<Dest) ...
                    Dest CfgP(n).up(CfgP(n).up>Dest)];
            end
        end
    end
end

end