% --- Two site fermion correlation matrix element function ---

function [Diff, OpMatEls] = PpPm_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a fermion configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements
% <Cfg|C*{i,dn}C*{i,up}C{j,up}C{j,dn}|CfgP>. Only the limited number of
% "reachable" configurations, where the matrix element is non-zero, are
% returned.

Cfg_vec = HilbertObj.FullCfg(Cfg);  % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds;
Nbond = numel(Bonds); % Count the number of bonds provided.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% There will be at most 2*Nbond reachable configurations.
OpMatEls = zeros(2*Nbond,1);

% Initialise a structure array by populating the final element first:
Diff(2*Nbond).pos = 1; % To be safe put one site in..
Diff(2*Nbond).val = 0;
Diff(2*Nbond).num = 0; % No differences.
Diff(2*Nbond).sign = 1; % For later assignment.
Diff(2*Nbond).type = 2; % Type 2 refers to a pair move.
% Now initialise the rest:
for n=1:2*Nbond
    Diff(n).pos = zeros(1,2); % Array of positions of configuration differences.
    Diff(n).val = zeros(2,2); % Array of difference values at those positions.
    Diff(n).num = 2; % Pair hopping terms count as two differences.
    Diff(n).sign = 1;
    Diff(n).type = 2;
end
for c = 1:CoOrd
    for n = 1:N
        m = Bonds(n,c); DInd = 2*CoOrd*(n-1) + 2*c - 1; % Difference index.
        if m ~= 0 && m ~= n
            Diff(DInd).pos(1) = n;
            Diff(DInd).val(:,1) = [-1; -1]; % Remove fermion pair from site n.
            Diff(DInd).pos(2) = m;
            Diff(DInd).val(:,2) = [1; 1]; % Add fermion pair to site m.
            OpMatEls(DInd) = prod(Cfg_vec(n,:)) * (sum(Cfg_vec(m,:))==0);
            % Require pair at site n and empty site m.
            
            % Hermitian conjugate / reverse move.
            Diff(DInd+1).type = 1; 
            Diff(DInd+1).pos(1) = m;
            Diff(DInd+1).val(:,1) = [-1; 1]; % Remove fermion pair from site m.
            Diff(DInd+1).pos(2) = n;
            Diff(DInd+1).val(:,2) = [1; 1]; % Add fermion pair to site n.
            OpMatEls(DInd+1) = prod(Cfg_vec(m,:)) * (sum(Cfg_vec(n,:))==0);
            % Require pair at site m and empty site n.
            
        elseif m == n % No pairs move, becomes total double occupancy counter.
            Diff(1).type = 1;
            Diff(1).pos = 1;
            Diff(1).val = 0;
            Diff(1).num = 0;
            OpMatEls(1) = sum(prod(Cfg_vec,2));
            break
        end
    end
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.