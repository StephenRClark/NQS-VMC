% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = MoveFermCfg(Cfg) 
% Randomly picks a fermion and moves it to any permissible location to
% propose a new configuration state of the system. This move proposal will 
% preserve the number of up and down fermions for any given input 
% configuration state.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------
% Fermion indices <= N_up correspond to up.
% Fermion indices > N_up correspond to down. 
% Convention - first position is site being vacated, second is destination.
% Difference values are in reference to differences in FullFermCfg.

Nf = numel(Cfg.up) + numel(Cfg.dn); Sites = 1:Cfg.N;
Ferm = randi(Nf); % Pick a random fermion.

CfgP = Cfg; 

if Ferm <= numel(Cfg.up)
    Start = Cfg.up(Ferm); Diff.type = 1;
    % Move to a site with no up fermions.
    Sites(Cfg.up) = []; Dest = Sites(randi(numel(Sites))); CfgP.up(Ferm) = [];
    CfgP.up = [CfgP.up(CfgP.up<Dest) Dest CfgP.up(CfgP.up>Dest)];    
    Diff.val = [-1 1];
    Diff.sign = (-1)^(sum( (Cfg.up > min(Start,Dest)) .* (Cfg.up < max(Start,Dest)) ));
else
    Start = Cfg.dn(Ferm-numel(Cfg.up)); Diff.type = -1;
    % Move to a site with no down fermions.
    Sites(Cfg.dn) = []; Dest = Sites(randi(numel(Sites))); CfgP.dn(Ferm-numel(Cfg.up)) = [];
    CfgP.dn = [CfgP.dn(CfgP.dn<Dest) Dest CfgP.dn(CfgP.dn>Dest)];
    Diff.val = [-1 1];
    Diff.sign = (-1)^(sum( (Cfg.dn > min(Start,Dest)) .* (Cfg.dn < max(Start,Dest)) ));
end
Diff.pos = [Start; Dest];
Diff.num = 1; % Only one fermion moves.