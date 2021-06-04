% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipSpin1Cfg(Cfg) 
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

CfgP = Cfg; % Proposed configuration.
N = Cfg.N; 
% All moves change two sites, type of change depends on site values.
Diff.num = 1; Diff.sign = 1; Vals = [-1 0 1];
Site = randi(N); Val1 = sum(Cfg.up == Site) - sum(Cfg.dn == Site);
Vals(Vals == Val1) = []; Val1P = Vals(randi(numel(Vals)));
Diff.pos = Site; Diff.val = Val1P - Val1;
% Adjust lists in CfgP.
if Val1 == 1
    CfgP.up(CfgP.up==Site) = [];
elseif Val1 == -1
    CfgP.dn(CfgP.dn==Site) = [];
end
if Val1P == 1
    CfgP.up = [CfgP.up(CfgP.up<Site), Site, CfgP.up(CfgP.up>Site)]; 
elseif Val1P == -1
    CfgP.dn = [CfgP.dn(CfgP.dn<Site), Site, CfgP.dn(CfgP.dn>Site)]; 
end
end