function Cfg = RandomBoseUnifFill(CParams)
% Creates a random configuration state for a system of N spinless bosons.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-0 bosons here.
% - Cfg.N = total number of sites in the system.
% - Cfg.Nb = total number of bosons in the system.
% - Cfg.occ = (N x 1) vector - boson occupation numbers by site.
% - Cfg.Nmax = maximum number of bosons on a single site.
% ---------------------------------

N = CParams.N; Cfg.N = N; Nb = CParams.Nb; Cfg.Nb = Nb; Cfg.Nmax = CParams.Nmax;

Cfg.occ = floor(Nb/N)*ones(N,1); Nrem = mod(Nb-N,N);

b = 0;
while b < Nrem
    site = randi(N);
    if Cfg.occ(site) < Cfg.Nmax
        Cfg.occ(site) = Cfg.occ(site) + 1;
        b = b + 1;
    end
end