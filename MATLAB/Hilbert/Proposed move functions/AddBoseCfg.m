% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = AddBoseCfg(Cfg)
% Randomly picks any site and changes its occupation within the range
% permitted by Cfg.Nmax. This move does not preserve the overall number of
% bosons in the system, and it is assumed that the Monte Carlo calculation
% is working in the grand canonical ensemble.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-0 bosons here.
% - Cfg.N = total number of sites in the system.
% - Cfg.Nb = total number of bosons in the system.
% - Cfg.occ = (N x 1) vector - boson occupation numbers by site.
% - Cfg.Nmax = maximum number of bosons on a single site.
% ---------------------------------

Cfg_vec = Cfg.occ; CfgP = Cfg; N = Cfg.N; Nmax = Cfg.Nmax;

Start = randi(N); % Pick a random site.
OldNum = Cfg_vec(Start); NewNum = OldNum;
while NewNum == OldNum
    NewNum = randi([0,Nmax]); % Ensure new proposed number is not the same.
end
Diff.num = 1; % One occupation number changes.
Diff.val = NewNum - OldNum; Diff.pos = Start;
Diff.Tfac = 1; % Trial probability in both directions is equal - 1/N*(Nmax+1).
CfgP.occ(Start) = NewNum;
end