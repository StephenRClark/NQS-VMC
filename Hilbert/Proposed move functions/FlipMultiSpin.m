% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipMultiSpin(Cfg)
% Randomly picks a small number of sites in the configuration and flips
% their spins to propose a new configuration state of the system.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

% Select between 1 and N/e sites to flip.

Nflip = randi(round(Cfg.N/exp(1))); Diff.num = Nflip;
Diff.pos = zeros(1,Nflip); Diff.val = zeros(1,Nflip);

for f = 1:Nflip
    ind_flip = randi(Cfg.N); % Pick a random site.
    while sum(Diff.pos == ind_flip)>0
        ind_flip = randi(Cfg.N);
    end
    Diff.pos(f) = ind_flip;
end

CfgP = Cfg; % The proposed configuration state.
for f = 1:Nflip
    % Find out if site is up or down.
    val = sum(Cfg.up==Diff.pos(f)) - sum(Cfg.dn==Diff.pos(f));
    Diff.val(f) = -2*val;
    if val == 1
        CfgP.up(CfgP.up==Diff.pos(f)) = [];
        CfgP.dn = [CfgP.dn(CfgP.dn<Diff.pos(f)) Diff.pos(f) CfgP.dn(CfgP.dn>Diff.pos(f))];
    elseif val == -1
        CfgP.dn(CfgP.dn==Diff.pos(f)) = [];
        CfgP.up = [CfgP.up(CfgP.up<Diff.pos(f)) Diff.pos(f) CfgP.up(CfgP.up>Diff.pos(f))];
    end
end
Diff.sign = 1; % For compatibility with fermionic implementations.
end