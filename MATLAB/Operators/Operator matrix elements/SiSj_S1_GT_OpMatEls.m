% --- Two site spin correlation matrix element function ---

function [Diff, OpMatEls] = SiSj_S1_GT_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a spin-1 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|S{i}.S{j}|CfgP>. Only the
% limited number of "reachable" configurations, where the matrix element is
% non-zero, are returned.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds; SLInds = GraphObj.SLInds; 
Nbond = numel(Bonds); % Count the number of bonds provided.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% There will be at most 2*Nbond+1 reachable configurations.
OpMatEls = zeros(2*Nbond+1,1);
% Initialise a structure array by populating the final element first:
Diff(2*Nbond+1).pos = 1; % To be safe put one site in..
Diff(2*Nbond+1).val = 0; %
Diff(2*Nbond+1).num = 0; % No differences.
Diff(2*Nbond+1).sign = 1; % For compatibility with fermionic implementations.
Diff(2*Nbond+1).type = 0; % Also for compatibility with fermionic implementations.
% Now initialise the rest:
for n=1:2*Nbond
    Diff(n).pos = zeros(1,2); % Array of positions of configuration differences.
    Diff(n).val = zeros(1,2); % Array of difference values at those positions.
    Diff(n).num = 2; % Heisenberg terms create at most 2 two differences.
    Diff(n).sign = 1;
    Diff(n).type = 0;
end

% Deal with each XX + YY contributions first:
% Compute all configurations obtained from Cfg by a pair of
% nearest-neighbour spin flips.
% Loop over sites and note the difference to the configuration:
for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); DInd = CoOrd*(n-1) + 2*c - 1; % Difference index.
        % Compute the differences CfgP - Cfg:
        if m ~= 0
            Diff(DInd).pos(1) = n;
            Diff(DInd).val(1) = -1; % Lower the spin at site n.
            Diff(DInd).pos(2) = m;
            Diff(DInd).val(2) = 1; % Raise the spin at site m.
            OpMatEls(DInd) = (-1)^(SLInds(n) - SLInds(m)) * ...
                (Cfg_vec(n)>=0) * (Cfg_vec(m)<1); % Compute matrix element.
        end
        
        % Conjugate moves:
        if m ~= 0
            Diff(DInd+1).pos(1) = m;
            Diff(DInd+1).val(1) = -1; % Lower the spin at site m.
            Diff(DInd+1).pos(2) = n;
            Diff(DInd+1).val(2) = 1; % Raise the spin at site n.
            OpMatEls(DInd+1) = (-1)^(SLInds(n) - SLInds(m)) * ...
                (Cfg_vec(m)>=0) * (Cfg_vec(n)<1); % Compute matrix element.
        end
    end
end
% Deal with the ZZ contribution last.
% This is a diagonal matrix element for H so there is no configuration difference
% in this case for any of the N terms. We already defined Diff(N+1) above.
% Loop over sites to add up the single contribution:
for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); % Nearest-neighbour site to m.
        if m ~= 0
            OpMatEls(2*Nbond+1) = OpMatEls(2*Nbond+1) + Cfg_vec(n)*Cfg_vec(m); % Compute matrix element.
        end
    end
end
% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.
end