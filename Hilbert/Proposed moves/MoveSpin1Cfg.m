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

Type = [1 2 3 4];

if N_up == 0 && N_dn == 0 % All spins in Sz = 0.
    Type = 1; % Only pair promotion yields any change.
else
    if N_up == 0 || N_dn == 0
        Type = [1 3]; % If all Sz in one direction, cannot use PA or SS.
    else
        if N_up + N_dn == N % All spins in Sz =/= 0.
            Type = [2 4]; % No Sz = 0 sites, cannot use PP or SM.
        elseif N - N_up - N_dn == 1 % One Sz = 0 site.
            Type = [2 3 4]; % Not enough Sz = 0 sites for PP.
        end
    end
end

Type = Type(randi(numel(Type))); % Randomly select from viable moves.

Diff.type = Type; Diff.sign = 1; % Placeholder for fermionic treatments.

Diff.num = 2; % All moves change two Sz values.

if Type == 1 % Pair promotion.
    sites = 1:N; sites([Cfg.up Cfg.dn]) = []; % Remove non-zero sites.
    
    SPlus = sites(randi(numel(sites))); % First site selected is given +1.
    sites(sites==SPlus) = []; % Remove chosen site.
    SMinus = sites(randi(numel(sites))); % Second site selected is given -1.
    
    Diff.val = [1, -1]; Diff.pos = [SPlus, SMinus]; 
    
    % Add new sites to CfgP lists.
    CfgP.up = [CfgP.up(CfgP.up<SPlus), SPlus, CfgP.up(CfgP.up>SPlus)]; 
    CfgP.dn = [CfgP.dn(CfgP.dn<SMinus), SMinus, CfgP.dn(CfgP.dn>SMinus)];
    
elseif Type == 2 % Pair annihilation.
    site_up = Cfg.up(randi(N_up)); site_dn = Cfg.dn(randi(N_dn)); % Select sites.
    
    Diff.val = [-1, 1]; Diff.pos = [site_up, site_dn]; 
    
    % Remove sites from CfgP lists.
    CfgP.up(CfgP.up == site_up) = []; CfgP.dn(CfgP.dn == site_dn) = [];
    
elseif Type == 3 % Spin move.
    sites_1 = [CfgP.up, CfgP.dn]; sites_0 = 1:N; sites_0(sites_1) = [];
    
    S1 = sites_1(randi(numel(sites_1))); S0 = sites_0(randi(numel(sites_0)));
    
    SM = (-1)^(sum(CfgP.up==S1)-1); % Determine if selected moving spin is +/-
    
    Diff.pos = [-1, 1] * SM; Diff.pos = [S1, S0];
    
    if SM == 1 % Up spin moves.
        CfgP.up(CfgP.up==S1) = [];
        CfgP.up = [CfgP.up(CfgP.up<S0), S0, CfgP.up(CfgP.up>S0)];
        Diff.val = [-1, 1];
    elseif SM == -1 % Down spin moves.
        CfgP.dn(CfgP.dn==S1) = [];
        CfgP.dn = [CfgP.dn(CfgP.dn<S0), S0, CfgP.dn(CfgP.dn>S0)];
        Diff.val = [1, -1];
    end

elseif Type == 4 % Spin swap.
    site_up = CfgP.up(randi(numel(CfgP.up))); site_dn = CfgP.dn(randi(numel(CfgP.dn))); 
    
    % Remove swapped sites from CfgP lists.
    CfgP.up(CfgP.up==site_up) = []; CfgP.dn(CfgP.dn==site_dn) = [];
    
    Diff.val = [-2, 2]; Diff.pos = [site_up, site_dn];
    
    % Add swapped sites to new CfgP lists.
    CfgP.up = [CfgP.up(CfgP.up<site_dn), site_dn, CfgP.up(CfgP.up>site_dn)];
    CfgP.dn = [CfgP.dn(CfgP.dn<site_up), site_up, CfgP.dn(CfgP.dn>site_up)];
end