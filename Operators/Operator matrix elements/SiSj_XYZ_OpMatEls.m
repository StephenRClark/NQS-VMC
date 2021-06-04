% --- Two site spin correlation matrix element function ---

function [Diff, OpMatEls] = SiSj_XYZ_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a spin-1 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|S{i}.S{j}|CfgP>. Only the
% limited number of "reachable" configurations, where the matrix element is
% non-zero, are returned.

% This operator acts in the xyz basis, with the numerical mapping:
% x -> +1, y -> 0, z -> -1.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds;
Nbond = numel(Bonds); % Count the number of bonds provided.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.
ValTable = [1, 2; 1, -1; -2, -1]; % Value table to help with pair promotion moves.
% Table format:
% [z->y z->x]
% [y->x y->z]
% [x->z x->y]

% There will be at most 2*Nbond reachable configurations.
OpMatEls = zeros(2*Nbond,1);
% Initialise a structure array by populating the final element first:
Diff(2*Nbond).pos = 1; % To be safe put one site in..
Diff(2*Nbond).val = 0; %
Diff(2*Nbond).num = 0; % No differences.
Diff(2*Nbond).sign = 1; % For compatibility with fermionic implementations.
Diff(2*Nbond).type = 0; % Also for compatibility with fermionic implementations.
% Now initialise the rest:
for n=1:2*Nbond
    Diff(n).pos = zeros(1,2); % Array of positions of configuration differences.
    Diff(n).val = zeros(1,2); % Array of difference values at those positions.
    Diff(n).num = 2; % Operator terms create 2 two differences.
    Diff(n).sign = 1;
    Diff(n).type = 0;
end

% Deal with each XX + YY contributions first:
% Compute all configurations obtained from Cfg by a pair of
% nearest-neighbour spin flips.
% Loop over sites and note the difference to the configuration:
for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); DInd = 2*CoOrd*(n-1) + 2*c - 1; % Difference index.
        % Compute the differences CfgP - Cfg:
        if m ~= 0
            % In xyz basis, this operator either swaps the spins if
            % dissimilar or changes both to the same spin if spins start
            % similar.            
            if (Cfg_vec(n) ~= Cfg_vec(m)) % Swap spin values.
                Diff(DInd).pos(1) = n;
                Diff(DInd).val(1) = Cfg_vec(m) - Cfg_vec(n);
                Diff(DInd).pos(2) = m;
                Diff(DInd).val(2) = Cfg_vec(n) - Cfg_vec(m);
                OpMatEls(DInd) = 1; % Compute matrix element.
            else
                TabInd = Cfg_vec(n) + 2;
                % First pair change - Diff value found in ValTable.
                Diff(DInd).pos(1) = n;
                Diff(DInd).val(1) = ValTable(TabInd,1);
                Diff(DInd).pos(2) = m;
                Diff(DInd).val(2) = ValTable(TabInd,1);
                OpMatEls(DInd) = -1; % Compute matrix element.
                % Second pair change - Diff value found in ValTable.
                Diff(DInd+1).pos(1) = n;
                Diff(DInd+1).val(1) = ValTable(TabInd,2);
                Diff(DInd+1).pos(2) = m;
                Diff(DInd+1).val(2) = ValTable(TabInd,2);
                OpMatEls(DInd+1) = -1; % Compute matrix element.
            end
        end
    end
end
% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.

end