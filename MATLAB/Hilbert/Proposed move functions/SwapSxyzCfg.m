% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = SwapSxyzCfg(Cfg)
% Randomly selects one of the potential S = 1 configuration Markov chain
% moves in the xyz basis, depending on which are viable. All these moves
% preserve the total spin parity of the configuration.
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
Diff.num = 2; Diff.sign = 1;
Site1 = randi(N);
Site2 = randi(N);
while Site2 == Site1
    Site2 = randi(N);
end
Diff.pos = [Site1, Site2];
% Determine what the values of each site are.
Val1 = sum(Cfg.up == Site1) - sum(Cfg.dn == Site1);
Val2 = sum(Cfg.up == Site2) - sum(Cfg.dn == Site2);
if (Val1 == Val2) % Same values - change both to same end value.
    Vals = [-1, 0, 1]; Vals(Vals == Val1) = []; % Must change value.
    Endval = Vals(randi(2));
    Diff.val = [Endval - Val1, Endval - Val1];
    Val1P = Endval; Val2P = Endval;
else % Dissimilar - swap their values.
    Diff.val = [(Val2 - Val1), (Val1 - Val2)];
    Val1P = Val2; Val2P = Val1;
end
% Adjust lists in CfgP.
if Val1 == 1
    CfgP.up(CfgP.up==Site1) = [];
elseif Val1 == -1
    CfgP.dn(CfgP.dn==Site1) = [];
end
if Val2 == 1
    CfgP.up(CfgP.up==Site2) = [];
elseif Val2 == -1
    CfgP.dn(CfgP.dn==Site2) = [];
end
if Val1P == 1
    CfgP.up = [CfgP.up(CfgP.up<Site1), Site1, CfgP.up(CfgP.up>Site1)]; 
elseif Val1P == -1
    CfgP.dn = [CfgP.dn(CfgP.dn<Site1), Site1, CfgP.dn(CfgP.dn>Site1)]; 
end
if Val2P == 1
    CfgP.up = [CfgP.up(CfgP.up<Site2), Site2, CfgP.up(CfgP.up>Site2)]; 
elseif Val2P == -1
    CfgP.dn = [CfgP.dn(CfgP.dn<Site2), Site2, CfgP.dn(CfgP.dn>Site2)]; 
end
Diff.Tfac = 1; % Trial probability in both directions is equal.
end