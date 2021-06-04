% --- Fermionic hopping matrix element function ---

function [Diff, OpMatEls] = CpCm_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a spin-1/2 fermionic configuration Cfg this function computes the
% list of configurations CfgP's and matrix elements <Cfg|C*{i}C{j}|CfgP>.
% Only the limited number of "reachable" configurations, where the matrix
% element is non-zero, are returned. UD signifies

% As operator acts to left, the conjugate is applied to the supplied
% configuration.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds;
Nbond = numel(Bonds); CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% There will be at most 4*Nbond reachable configurations.
OpMatEls = zeros(4*Nbond,1);
% Initialise a structure array by populating the final element first:
Diff(4*Nbond).pos = 1; % To be safe put one site in..
Diff(4*Nbond).val = 0; %
Diff(4*Nbond).num = 0; % No differences.
Diff(4*Nbond).sign = 1; % For later assignment.
Diff(4*Nbond).type = 0; % Irrelevant for zero difference, but needed for other elements.
% Now initialise the rest:
for n=1:4*Nbond
    Diff(n).pos = zeros(1,2); % Array of positions of configuration differences.
    Diff(n).val = zeros(1,2); % Array of difference values at those positions.
    Diff(n).num = 1; % Fermion moves count as single differences.
    Diff(n).sign = 1;
    Diff(n).type = 0; % For later assignment.
end

for c = 1:CoOrd
    for n = 1:N
        m = Bonds(n,c); DInd = 2*CoOrd*(n-1) + 2*c - 1; % Difference index.
        if m ~= 0 && m ~= n
            
            % Up fermion moves.
            
            Diff(DInd).type = 1; % Signifies up fermion.
            Diff(DInd).pos(1) = m;
            Diff(DInd).val(1) = -1; % Remove up fermion from site m.
            Diff(DInd).pos(2) = n;
            Diff(DInd).val(2) = 1; % Add up fermion to site n.
            OpMatEls(DInd) = (1 - Cfg_vec(m,1)) * Cfg_vec(n,1);
            % Zero if site m already occupied by up fermion or site n empty.
            
            % Hermitian conjugate / reverse move.
            Diff(DInd+1).type = 1; % Signifies up fermion.
            Diff(DInd+1).pos(1) = n;
            Diff(DInd+1).val(1) = -1; % Remove up fermion from site n.
            Diff(DInd+1).pos(2) = m;
            Diff(DInd+1).val(2) = 1; % Add up fermion to site m.
            OpMatEls(DInd+1) = (1 - Cfg_vec(n,1)) * Cfg_vec(m,1);
            % Zero if site n already occupied by up fermion or site m empty.
            
            % Down fermion moves.
            
            DInd = 2*Nbond + 2*CoOrd*(n-1) + 2*c - 1; % Difference index.
            
            Diff(DInd).type = -1; % Signifies down fermion.
            Diff(DInd).pos(1) = m;
            Diff(DInd).val(1) = -1; % Remove down fermion from site n.
            Diff(DInd).pos(2) = n;
            Diff(DInd).val(2) = 1; % Add down fermion to site m.
            OpMatEls(DInd) = (1 - Cfg_vec(m,2)) * Cfg_vec(n,2);
            % Zero if site m already occupied by down fermion or site n empty.
            
            % Hermitian conjugate / reverse move.
            Diff(DInd+1).type = -1; % Signifies down fermion.
            Diff(DInd+1).pos(1) = n;
            Diff(DInd+1).val(1) = -1; % Remove down fermion from site n.
            Diff(DInd+1).pos(2) = m;
            Diff(DInd+1).val(2) = 1; % Add down fermion to site m.
            OpMatEls(DInd+1) = (1 - Cfg_vec(n,2)) * Cfg_vec(m,2);
            % Zero if site n already occupied by down fermion or site m empty.
            
        elseif m == n % No fermions move, becomes number counter.
            Diff(DInd).type = 1;
            Diff(DInd).pos = m;
            Diff(DInd).val = 0;
            Diff(DInd).num = 0;
            OpMatEls(DInd) = Cfg_vec(m,1);
            
            Diff(DInd+1).type = -1;
            Diff(DInd+1).pos = 1;
            Diff(DInd+1).val = 0;
            Diff(DInd+1).num = 0;
            OpMatEls(DInd+1) = Cfg_vec(m,2);
        end
    end
end
% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.

