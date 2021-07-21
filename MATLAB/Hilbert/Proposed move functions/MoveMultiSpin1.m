% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = MoveMultiSpin1(Cfg)
% Randomly select a subset of sites and change their spin values, keeping
% total spin projection.
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
Cfg_vec = zeros(Cfg.N,1); Cfg_vec(Cfg.up) = 1; Cfg_vec(Cfg.dn) = -1;
N = Cfg.N; N_up = numel(Cfg.up); N_dn = numel(Cfg.dn);
% Can change either two or three sites.
Diff.sign = 1; Vals = [-1 0 1];
Site1 = randi(N); Val1 = sum(Cfg.up == Site1) - sum(Cfg.dn == Site1);
Vals(Vals==Val1) = []; Val1P = Vals(randi(2)); Delta = Val1P - Val1;
switch Delta
    case -2
        Sites = Cfg.dn;
        if isempty(Sites)
            Diff.num = 3;
        elseif (N-N_up)<2
            Diff.num = 2;
        else
            Diff.num = 1+randi(2);
        end
        if Diff.num == 3
            Sites = 1:N; Sites(Cfg.up) = []; Sites(Sites==Site1) = [];
        end
    case 2
        Sites = Cfg.up;
        if isempty(Sites)
            Diff.num = 3;
        elseif (N-N_dn)<2
            Diff.num = 2;
        else
            Diff.num = 1+randi(2);
        end
        if Diff.num == 3
            Sites = 1:N; Sites(Cfg.dn) = []; Sites(Sites==Site1) = [];
        end
    case -1
        Diff.num = 2; Sites = 1:N; Sites(Cfg.up) = []; Sites(Sites==Site1) = [];
    case 1
        Diff.num = 2; Sites = 1:N; Sites(Cfg.dn) = []; Sites(Sites==Site1) = [];
end
Diff.pos = zeros(1,Diff.num); Diff.pos(1) = Site1;
for d = 2:Diff.num
    Diff.pos(d) = Sites(randi(numel(Sites))); Sites(Sites==Diff.pos(d)) = [];
end
if Diff.num == 3
    Diff.val = [Delta -sign(Delta) -sign(Delta)];
else
    Diff.val = [Delta -Delta];
end
Vals0 = zeros(1,Diff.num);
for d = 1:Diff.num
    Vals0(d) = Cfg_vec(Diff.pos(d));
end
ValsP = Vals0 + Diff.val;
% Adjust lists in CfgP.
for d = 1:Diff.num
    if Vals0(d) == 1
        CfgP.up(CfgP.up==Diff.pos(d)) = [];
    elseif Vals0(d) == -1
        CfgP.dn(CfgP.dn==Diff.pos(d)) = [];
    end
    if ValsP(d) == 1
        CfgP.up = [CfgP.up(CfgP.up<Diff.pos(d)), Diff,pos(d), CfgP.up(CfgP.up>Diff.pos(d))];
    elseif ValsP(d) == -1
        CfgP.dn = [CfgP.dn(CfgP.dn<Diff.pos(d)), Diff.pos(d), CfgP.dn(CfgP.dn>Diff.pos(d))];
    end
end
N_upP = numel(CfgP.up); N_dnP = numel(CfgP.dn);
switch Diff.num % Trial probabilities treated differently depending on number.
    case 3
        if Delta < 0
            N_ex = N_up; N_exP = N_dnP;
        else
            N_ex = N_dn; N_exP = N_upP;
        end
        Diff.Tfac = ((N-N_exP)*(N-N_exP-1))/((N-N_ex)*(N-N_ex-1));
    case 2
        Diff.Tfac = ((N-N_upP)*(N-N_dnP)*(2*N-N_up-N_dn)) / ...
            ((N-N_up)*(N-N_dn)*(2*N-N_upP-N_dnP));
end
end