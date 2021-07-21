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
N = Cfg.N; N_up = numel(Cfg.up); N_dn = numel(Cfg.dn);
% All moves change two sites, type of change depends on site values.
Diff.num = 2; Diff.sign = 1; Vals = [-1 0 1];
Site1 = randi(N); Val1 = sum(Cfg.up == Site1) - sum(Cfg.dn == Site1);
Vals(Vals==Val1) = []; Val1P = Vals(randi(2)); Delta = Val1P - Val1;
switch Delta
    case -2
        Sites = Cfg.dn;
        if isempty(Sites)
            Delta = -1; Sites = 1:N; Sites(Cfg.up) = []; Sites(Sites==Site1) = [];
        end
    case 2
        Sites = Cfg.up;
        if isempty(Sites)
            Delta = 1; Sites = 1:N; Sites(Cfg.dn) = []; Sites(Sites==Site1) = [];
        end
    case -1
        Sites = 1:N; Sites(Cfg.up) = []; Sites(Sites==Site1) = [];
    case 1
        Sites = 1:N; Sites(Cfg.dn) = []; Sites(Sites==Site1) = [];
end
Site2 = Sites(randi(numel(Sites)));
Val2 = sum(Cfg.up == Site2) - sum(Cfg.dn == Site2); Val2P = Val2 - Delta;
Diff.pos = [Site1, Site2]; Diff.val = [Delta, -Delta];
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
N_upP = numel(CfgP.up); N_dnP = numel(CfgP.dn);
% Trial probability ratio depends on number of up / dn sites.
Diff.Tfac = ((N-N_upP)*(N-N_dnP)*(2*N-N_up-N_dn)) / ...
    ((N-N_up)*(N-N_dn)*(2*N-N_upP-N_dnP));
end