% --- Difference to configuration conversion function ---

function [CfgP] = Diff2CfgBose(Diff,Cfg)
% This function converts a Diff struct compatible with bosonic moves to a
% CfgP struct when provided a starting Cfg.

% Assumption is that Diff will represent either no moves, boson moves or
% boson additions. If multiple Diffs are passed, assume the Diffs are in a
% Nd x 1 struct array.

Nd = numel(Diff); % Function should be able to handle multiple Diffs being passed.

% Initialise structure array CfgP with last element first.
CfgP(Nd) = Cfg; CfgP = CfgP.'; % Transpose to output Nd x 1 struct array.

for n = 1:Nd
    CfgP(n) = Cfg;
    % Alter each CfgP using associated Diff. In bosonic case, Diff.val
    % records site occupation differences.
    for d = 1:Diff(n).num
        CfgP(n).occ(Diff(n).pos(d)) = Cfg.occ(Diff(n).pos(d)) + Diff(n).val(d);
    end
end

end