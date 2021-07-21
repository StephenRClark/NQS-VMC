% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipMultiSpinGS(Cfg)
% Randomly picks a small number of sites in the configuration and flips
% their spins to propose a new configuration state of the system.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

% Select 2 or 3 sites to flip.

Nflip = 1 + randi(2); % 2 or 3 spin flips.
Diff.num = Nflip;
if Nflip == 2 % Two sites separated by 3
    Diff.pos = [randi(Cfg.N) 0];
    Diff.pos(2) = 1+mod(Diff.pos(1)+2,Cfg.N);
elseif Nflip == 3 % Consecutive trio of sites.
    Diff.pos = [randi(Cfg.N) 0 0];
    Diff.pos(2) = 1+mod(Diff.pos(1),Cfg.N);
    Diff.pos(3) = 1+mod(Diff.pos(1)+1,Cfg.N);
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
Diff.Tfac = 1; % Trial probability in both directions is equal.
end