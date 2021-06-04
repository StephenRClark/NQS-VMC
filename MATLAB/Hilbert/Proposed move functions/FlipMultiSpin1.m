% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipMultiSpin1(Cfg) 
% Randomly selects one site to change its spin projection in the Sz basis.
% This move does not preserve the total spin projection.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.SzT = total configuration Sz - if set to empty, will be treated as
% indefinite, and requires other configuration move function.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

% Select between 1 and N/e sites to flip.

Nflip = randi(round(Cfg.N/exp(1))); Diff.num = Nflip; vals = [-1 0 1];
Diff.pos = zeros(1,Nflip); Diff.val = zeros(1,Nflip);
for f = 1:Nflip
    ind_flip = randi(Cfg.N); % Pick a random site.
    while sum(Diff.pos == ind_flip)>0
        ind_flip = randi(Cfg.N);
    end
    Diff.pos(f) = ind_flip;
end
CfgP = Cfg; % Proposed configuration.
for f = 1:Nflip
    % Find out if site is up or down.
    val = sum(Cfg.up==Diff.pos(f)) - sum(Cfg.dn==Diff.pos(f));
    val_p = vals(randi(3));
    while val_p == val
        val_p = vals(randi(3));
    end
    Diff.val(f) = val_p - val;
    if val == 1
        CfgP.up(CfgP.up==Diff.pos(f)) = [];
    elseif val == -1
        CfgP.dn(CfgP.dn==Diff.pos(f)) = [];
    end
    if val_p == 1
        CfgP.up = [CfgP.up(CfgP.up<Diff.pos(f)) Diff.pos(f) CfgP.up(CfgP.up>Diff.pos(f))];
    elseif val_p == -1
        CfgP.dn = [CfgP.dn(CfgP.dn<Diff.pos(f)) Diff.pos(f) CfgP.dn(CfgP.dn>Diff.pos(f))];
    end
end
Diff.sign = 1; % For compatibility with fermionic implementations.
end