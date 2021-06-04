function Cfg = RandomSpin1FixedMag(CParams)
% Creates a random configuration state for a system of N spin-1s 
% with zero magnetisation.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.SzT = total configuration Sz.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

N = CParams.N; SzT = CParams.SzT;

if abs(SzT) >= N
    error('Requested SzT gives fully aligned state - not much point in that, is there?')
end

Cfg.N = N; Cfg.SzT = SzT;

if SzT ~= 0 % Non-zero total Sz projection.
    Cfg_vec = ones(N,1) * sign(SzT); Sd = N - abs(SzT); % Sd - spin difference.
    while Sd > 0 
        site = randi(N); % Choose a site randomly
        if abs(Cfg_vec(site)- sign(SzT))<=1 % Avoid unphysical spin on site.
            Cfg_vec(site) = Cfg_vec(site) - sign(SzT);
            Sd = Sd - 1;
        end
    end
    Cfg.up = find(Cfg_vec==1).'; Cfg.dn = find(Cfg_vec==-1).';    
else % Zero total Sz projection.
    Np = randi(floor(N/2)); % Choose some number of +/- pairs to make.
    sites = randperm(N,2*Np); % Randomly permute all the site numbers.
    Cfg.up = sort(sites(1:(Np))); % Assign first half to up spins.
    Cfg.dn = sort(sites((Np+1):(2*Np))); % Assign second half to down spins.
end

end