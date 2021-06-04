% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = MoveSpin1Cfg(Cfg) 
% Randomly selects one of four potential S = 1 configuration Markov chain
% moves, depending on which are viable. All these moves preserve the total
% spin projection of the configuration.
% Type 1 - PP - Pair promotion - two Sz = 0 sites are given +1 and -1.
% Type 2 - PA - Pair annihilation - two Sz = +/- 1 sites cancel to Sz = 0.
% Type 3 - SM - Spin move - a Sz = +/- 1 site swaps with a Sz = 0 site.
% Type 4 - SS - Spin swap - two Sz = +/- 1 sites swap Sz values.
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
Site1 = randi(N); Val1 = sum(Cfg.up == Site1) - sum(Cfg.dn == Site1);
Site2 = randi(N); Val2 = sum(Cfg.up == Site2) - sum(Cfg.dn == Site2);
% Select a pair of sites that either differ in spin or sum to zero.
while (Site2 == Site1) || ((Val2 == Val1) && ((Val1 + Val2) ~= 0))
    Site2 = randi(N); Val2 = sum(Cfg.up == Site2) - sum(Cfg.dn == Site2);
end
Diff.pos = [Site1, Site2];
if (Val1 == Val2) && (Val1 == 0) % Both zero - pair promotion.
    Val1P = (-1)^(randi(2)); Val2P = -Val1P;
elseif (Val1 == -Val2) && (Val1 ~= 0) % +/- pair;
    Val1P = - Val1 * randi(2); Val2P = -Val1P; % Either swaps or zeros the pair. 
else % Simply swap the pair.
    Val1P = Val2; Val2P = Val1;
end
Diff.val = [(Val1P - Val1), (Val2P - Val2)];
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
end