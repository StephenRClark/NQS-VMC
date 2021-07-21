% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipUnaryCfg(Cfg) 
% Randomly picks one site in the configuration and flips its spin to 
% propose a new configuration state of the system.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% - Cfg.N0 = original number of sites.
% - Cfg.d = original site dimension.
% ---------------------------------

s0 = randi(Cfg.N0); % Randomly select original site index.
ind0 = sum(Cfg.up <= s0*Cfg.d); v0 = 1+mod(Cfg.up(ind0)-1,Cfg.d);
vp = randi(Cfg.d); % Randomly select new value.
while vp == v0
    vp = randi(Cfg.d);
end
su = v0 + (s0-1)*Cfg.d; sd = vp + (s0-1)*Cfg.d;
% Set up Diff struct.
Diff.num = 2; Diff.sign = 1; Diff.val = [-2 2];
Diff.pos = [su sd];
% Adjust up and down lists in Cfg.
CfgP = Cfg; CfgP.up(ind0) = sd; CfgP.dn(CfgP.dn==sd) = [];
CfgP.dn = [CfgP.dn(CfgP.dn<su) su CfgP.dn(CfgP.dn>su)];
Diff.Tfac = 1; % Trial probability in both direction is equal - 1/N*(d-1)
end