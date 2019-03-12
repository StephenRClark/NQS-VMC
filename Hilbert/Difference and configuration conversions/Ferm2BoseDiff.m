% --- Difference structure conversion function ---

function [Diff] = Ferm2BoseDiff(Diff)
% This function converts a Diff struct compatible with fermionic moves to a
% Diff struct compatible with bosonic ansatz.

% Assumption is that Diff will represent either no moves, single fermion
% movements or fermion site swaps.

if Diff.num~=0
    if abs(Diff.type) == 1 % Fermion has moved.
        Diff.val = [-1, 1]; Diff.num = 2; % Two occupation numbers change.
        % Single fermion moves always have Diff format [Start Destination]. 
    elseif Diff.type == 0 % Fermion / spin exchange.
        Diff.num = 0; % No difference in configuration for boson swaps.
    elseif Diff.type == 2 % Fermion pair move.
        Diff.val = [-2, 2]; Diff.num = 2; % Two occupation numbers change.
    end
end