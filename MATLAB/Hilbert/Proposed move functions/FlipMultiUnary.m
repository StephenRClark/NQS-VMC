% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipMultiUnary(Cfg)
% Randomly picks a random number of sites in the configuration and flips
% their spins to propose a new configuration state of the system.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% - Cfg.N0 = original number of sites.
% - Cfg.d = original site dimension.
% ---------------------------------

Nflip = randi(max([2,round(Cfg.N/exp(1))])); Diff.num = 2*Nflip;
Diff.pos = zeros(1,2*Nflip); Diff.val = zeros(1,2*Nflip); CfgP = Cfg;
s0list = zeros(1,Nflip);
for f = 1:Nflip
    ind_flip = randi(Cfg.N0); % Pick a random site.
    while sum(s0list == ind_flip)>0
        ind_flip = randi(Cfg.N0);
    end
    s0list(f) = ind_flip;
    ind0 = sum(Cfg.up <= ind_flip*Cfg.d); v0 = 1+mod(Cfg.up(ind0)-1,Cfg.d);
    vp = randi(Cfg.d); % Randomly select new value.
    while vp == v0
        vp = randi(Cfg.d);
    end
    su = v0 + (s0-1)*Cfg.d; sd = vp + (s0-1)*Cfg.d;
    % Set up Diff struct.
    Diff.val((1:2)+2*(f-1)) = [-2 2];
    Diff.pos((1:2)+2*(f-1)) = [su sd];
    % Adjust up and down lists in Cfg.
    CfgP.up(ind0) = sd; CfgP.dn(CfgP.dn==sd) = [];
    CfgP.dn = [CfgP.dn(CfgP.dn<su) su CfgP.dn(CfgP.dn>su)];
end
Diff.Tfac = 1; % Trial probability in both directions is equal.
end