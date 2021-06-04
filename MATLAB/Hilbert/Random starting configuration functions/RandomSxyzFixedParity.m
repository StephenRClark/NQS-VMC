function Cfg = RandomSxyzFixedParity(CParams)
% Creates a random configuration state for a system of N spin-1s in the xyz
% basis with a fixed parity set by CParams.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.SzT = total configuration Sz.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

N = CParams.N; Parity = CParams.Sector;
if (mod(N+Parity,2) ~= 0)
    error('Specified configuration parity not achievable with the given site number.');
end
Cfg_vec = zeros(N,1); SiteList = 1:N; % Knockout list.
if (mod(N,2) == 0) % Even case - ensure even number of each spin alignment.    
    while numel(SiteList) ~= 0
        Site1 = SiteList(randi(numel(SiteList))); 
        SiteList(SiteList == Site1) = []; % Remove site 1 from knockout list.
        if (numel(SiteList) > 1)
            Site2 = SiteList(randi(numel(SiteList)));
        else
            Site2 = SiteList;
        end
        SiteList(SiteList == Site2) = []; % Remove site 2 from knockout list.
        Val = randi(3) - 2; Cfg_vec([Site1,Site2]) = Val;
    end
else
    while numel(SiteList) > 3 % Add spin pairs until the last three.
        Site1 = SiteList(randi(numel(SiteList))); 
        SiteList(SiteList == Site1) = []; % Remove site 1 from knockout list.
        Site2 = SiteList(randi(numel(SiteList)));
        SiteList(SiteList == Site2) = []; % Remove site 2 from knockout list.
        Val = randi(3) - 2; Cfg_vec([Site1,Site2]) = Val;
    end
    Cfg_vec(SiteList) = randperm(3) - 2;
end
Cfg.up = find(Cfg_vec==1).'; Cfg.dn = find(Cfg_vec==-1).';
Cfg.N = N;
end