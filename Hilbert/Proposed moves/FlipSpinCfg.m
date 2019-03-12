% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = FlipSpinCfg(Cfg) 
% Randomly picks one site in the configuration and flips its spin to 
% propose a new configuration state of the system.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

ind_flip = randi(Cfg.N); % Pick a random site.
% Retain the difference information (CfgP - Cfg) for later use:
Diff.pos = ind_flip;

CfgP = Cfg; % The proposed configuration state.
% Reassign the indexed site to either CfgP.up or CfgP.dn:
if isempty(Cfg.up) % If start configuration is a vector of ones
    CfgP.dn(CfgP.dn==ind_flip) = []; % Remove ind_flip site from CfgP.dn
    CfgP.up = ind_flip; % CfgP.up now has an element - ind_flip
    Diff.val = 2; % Sz(i) goes from -1 to 1
elseif isempty(Cfg.dn) % If start configuration is a vector of zeros
    CfgP.up(CfgP.up==ind_flip) = []; % Remove ind_flip site from CfgP.up
    CfgP.dn = ind_flip; % CfgP.dn now has an element - ind_flip
    Diff.val = -2; % Sz(i) goes from 1 to -1
else
    % Find out if ind_flip is an up site or down site
    check = CfgP.up==ind_flip;
    if sum(check)==0
        CfgP.dn(CfgP.dn==ind_flip) = []; % Remove ind_flip from down site vector
        CfgP.up = [CfgP.up ind_flip]; % Add ind_flip to up site vector
        Diff.val = 2;
    else
        CfgP.up(check) = []; % Remove ind_flip from up site vector
        CfgP.dn = [CfgP.dn ind_flip]; % Add ind_flip to down site vector
        Diff.val = -2;
    end
end 
Diff.num = 1; % Only one difference between configurations
Diff.sign = 1; % For compatibility with fermionic implementations.