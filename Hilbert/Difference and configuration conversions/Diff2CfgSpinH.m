% --- Difference to configuration conversion function ---

function [CfgP] = Diff2CfgSpinH(Diff,Cfg)
% This function converts a Diff struct compatible with spin-1/2 moves to a
% CfgP struct when provided a starting Cfg.

% Assumption is that Diff will represent either no moves, single flips or
% spin swaps. If multiple Diffs are passed, assume the Diffs are in a
% Nd x 1 struct array.

Nd = numel(Diff); % Function should be able to handle multiple Diffs being passed.

% Initialise structure array CfgP with last element first.
CfgP(Nd) = Cfg; CfgP = CfgP.'; % Transpose to output Nd x 1 struct array.

for n = 1:Nd
    CfgP(n) = Cfg;
    % Alter each CfgP using associated Diff. In spin-1/2 case, Diff.val
    % records local Sz projection differences. This will require adjustment
    % of spin up/down lists.
    for d = 1:Diff(n).num
        if Diff(n).val(d) > 0 % Spin raised.
            % Remove site from down list.
            CfgP(n).dn(CfgP(n).dn==Diff(n).pos(d)) = [];
            % Add site to up list.
            CfgP(n).up = [CfgP(n).up(CfgP(n).up<Diff(n).pos(d)) ...
                Diff(n).pos(d) CfgP(n).up(CfgP(n).up>Diff(n).pos(d))];
        elseif Diff(n).val(d) < 0 % Spin lowered.
            % Remove site from up list.
            CfgP(n).up(CfgP(n).up==Diff(n).pos(d)) = [];
            % Add site to down list.
            CfgP(n).dn = [CfgP(n).dn(CfgP(n).dn<Diff(n).pos(d)) ...
                Diff(n).pos(d) CfgP(n).dn(CfgP(n).dn>Diff(n).pos(d))];
        end
    end
end

end