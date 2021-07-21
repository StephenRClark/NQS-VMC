% --- Monte Carlo move proposal function ---

function [Diff,CfgP] = SwapSpinCfg(Cfg) 
% Randomly picks two sites with different spin states and swaps their
% states to propose a new configuration state of the system. This move
% proposal will preserve the number of any given spin state in the
% input configuration state.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

ind_up = randi(numel(Cfg.up)); % Pick a random up spin site.
site_up = Cfg.up(ind_up);
ind_dn = randi(numel(Cfg.dn)); % Pick a random down spin site.
site_dn = Cfg.dn(ind_dn);

CfgP = Cfg; % The proposed configuration state.

% Exchange the sites in the lists for the proposed configuration, then 
% reorder by lattice position for fermionic use (does not affect others):
UpList = Cfg.up; DnList = Cfg.dn; % Temporary listings of indices
UpList(UpList==site_up)=[]; DnList(DnList==site_dn)=[]; % Remove sites being changed
CfgP.up = [UpList(UpList<site_dn) site_dn UpList(UpList>site_dn)];
CfgP.dn = [DnList(DnList<site_up) site_up DnList(DnList>site_up)];
% Retain the difference information (CfgP - Cfg) for later use:
Diff.pos = [site_up,site_dn];
Diff.val = [-2,+2]; 
Diff.num = 2; % Swap moves only ever create two differences.
% If fermionic, include a sign calculation by processing Diff with a
% separate function to alter the sign field.
Diff.sign = 1;
Diff.type = 0; % For compatibility with fermionic reference states.
Diff.Tfac = 1; % Trial probability in both directions is equal - 1/N_up*N_dn.
end