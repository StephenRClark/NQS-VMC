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

CfgP = Cfg; N = Cfg.N; ZeroSites = (1:N).'; Cfg_vec = zeros(N,1);
Cfg_vec(Cfg.up) = 1; Cfg_vec(Cfg.dn) = -1;
ZeroSites(Cfg.up) = 0; ZeroSites(Cfg.dn) = 0; ZeroSites(ZeroSites==0) = [];
PlusSites = Cfg.up.'; MinusSites = Cfg.dn.';
Nsite = 1+randi(floor((N-1)/2)); 
if sum(abs(Cfg_vec))==0
    Nsite = Nsite - mod(Nsite,2);
end
if (sum(abs(Cfg_vec)) == 1) || ((abs(sum(Cfg_vec)) == sum(abs(Cfg_vec))) && (sum(abs(Cfg_vec))>0))
    Sign = -sign(sum(Cfg_vec));
else
    Sign = (-1)^(randi(2));
end
DiffSites = zeros(1,Nsite); Shifts = zeros(1,Nsite); Vals0 = zeros(1,Nsite);
if mod(Nsite,2)
    Shifts(1) = 2 * Sign; Shifts([2 3]) = -Sign;
    for s = 4:Nsite
        Shifts(s) = (-1)^(s);
    end
else
    Shifts(1) = Sign; Shifts(2) = -Sign;
    for s = 3:Nsite
        Shifts(s) = (-1)^(s);
    end
end
for s = 1:Nsite
    if abs(Shifts(s)) == 2
        if Sign > 0
            Ind = randi(numel(MinusSites));
            DiffSites(s) = MinusSites(Ind);
            MinusSites(Ind) = []; Vals0(s) = -1;
        else
            Ind = randi(numel(PlusSites));
            DiffSites(s) = PlusSites(Ind);
            PlusSites(Ind) = []; Vals0(s) = 1;
        end
    else
        if sign(Shifts(s)) > 0
            SiteList = [ZeroSites; MinusSites];
            Ind = randi(numel(SiteList)); N0 = numel(ZeroSites);
            DiffSites(s) = SiteList(Ind); Vals0(s) = -(Ind>N0);
            if Ind > N0
                MinusSites(Ind-N0) = [];
            else
                ZeroSites(Ind) = [];
            end
        else
            SiteList = [ZeroSites; PlusSites];
            Ind = randi(numel(SiteList)); N0 = numel(ZeroSites);
            DiffSites(s) = SiteList(Ind); Vals0(s) = (Ind>N0);
            if Ind > N0
                PlusSites(Ind-N0) = [];
            else
                ZeroSites(Ind) = [];
            end
        end
    end
end
ValsP = Vals0 + Shifts;
% Adjust lists in CfgP.
Diff.pos = DiffSites; Diff.num = Nsite; Diff.val = Shifts; Diff.sign = 1;
for s = 1:Nsite
    if Vals0(s) == 1
        CfgP.up(CfgP.up==DiffSites(s)) = [];
    elseif Vals0(s) == -1
        CfgP.dn(CfgP.dn==DiffSites(s)) = [];
    end
    if ValsP(s) == 1
        CfgP.up = [CfgP.up(CfgP.up<DiffSites(s)), DiffSites(s), CfgP.up(CfgP.up>DiffSites(s))];
    elseif ValsP(s) == -1
        CfgP.dn = [CfgP.dn(CfgP.dn<DiffSites(s)), DiffSites(s), CfgP.dn(CfgP.dn>DiffSites(s))];
    end
end
end