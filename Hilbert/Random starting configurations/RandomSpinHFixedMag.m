function Cfg = RandomSpinHFixedMag(CParams)
% Creates a random configuration state for a system of N spin-1/2's 
% with zero magnetisation.
% ---------------------------------
% Format for configuration states:
% - Cfg.type = identifier for the type of states, assumed spin-1/2 here.
% - Cfg.N = total number of sites in the system.
% - Cfg.up = (Nup x 1) vector of sites where the spin is up.
% - Cfg.dn = (Ndn x 1) vector of sites where the spin is down.
% ---------------------------------

N = CParams.N; SzT = CParams.SzT;
N_up = round((N + SzT)/2); 
if mod(N,2)==0
  Cfg.N = N;
  sites = randperm(N); % Randomly permute all the site numbers.
  Cfg.up = sort(sites(1:N_up)); % Assign first half to up spins.
  Cfg.dn = sort(sites((N_up+1):N)); % Assign second half to down spins.
else
  error('RandomSpinFixedMag(N): Even number of sites N required.')
end 

