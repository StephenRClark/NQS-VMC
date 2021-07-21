% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = MoveBoseCfg(Cfg)
% Randomly picks an occupied site and moves a boson from it to any
% permissible location to propose a new configuration state of the system.
% This move proposal will preserve the number of bosons for any given input
% bosonic configuration state.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-0 bosons here.
% - Cfg.N = total number of sites in the system.
% - Cfg.Nb = total number of bosons in the system.
% - Cfg.occ = (N x 1) vector - boson occupation numbers by site.
% - Cfg.Nmax = maximum number of bosons on a single site.
% ---------------------------------

% Convention - first position is site being vacated, second is destination.
% Difference values are in reference to differences in Cfg.occ.

Cfg_vec = Cfg.occ; CfgP = Cfg; Nb = Cfg.Nb; N = Cfg.N; Nmax = Cfg.Nmax;

% Pick a start site that has bosons already.
SiteStart = find(Cfg_vec>0); Start = SiteStart(randi(numel(SiteStart)));
% Randomly select any site with under Nmax bosons.
SiteDest = find(Cfg_vec<Nmax); Dest = SiteDest(randi(numel(SiteDest)));
while Dest == Start
    Dest = SiteDest(randi(numel(SiteDest))); % Reroll if Start is selected.
end
Diff.num = 2; % Two occupation numbers change.
Diff.val = [-1, 1]; Diff.pos = [Start, Dest];
% Calculate trial probability ratio factor.
Nmax = sum(Cfg_vec == Nmax); Nocc = sum(Cfg_vec>0);
dNmax = (Cfg_vec(Dest) == (Nmax-1)) - (Cfg_vec(Start) == Nmax);
dNocc = (Cfg_vec(Dest) == 0) - (Cfg_vec(Start) == 1);
% Trial probability ratio depends on number of occupied sites.
Diff.Tfac = (Nocc+dNocc)/Nocc; 
for d = 1:2
    CfgP.occ(Diff.pos(d)) = CfgP.occ(Diff.pos(d)) + Diff.val(d);
end
end