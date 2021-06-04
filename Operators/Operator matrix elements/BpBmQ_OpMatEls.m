% --- Boson hopping matrix element function ---

function [Diff, OpMatEls] = BpBmQ_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a boson configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements
% <Cfg|exp{1i*qj}B+{i}*B{j}|CfgP>. Only the limited number of "reachable"
% configurations, where the matrix element is non-zero, are returned.

% As operator acts to left, the conjugate is applied to the supplied
% configuration.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds;
Nbond = numel(Bonds); CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% There will be at most Nbond reachable configurations.
OpMatEls = zeros(Nbond,N); q = ((1:N)-1)*2*pi/N;
% Each new configuration has N associated values (from N values of
% momentum q).

% Initialise a structure array by populating the final element first:
Diff(Nbond).pos = 1; % To be safe put one site in..
Diff(Nbond).val = 0; %
Diff(Nbond).num = 0; % No differences.
Diff(Nbond).type = 0; % For compatibility.
Diff(Nbond).sign = 1; % For compatibility.

% Now initialise the rest:
for n=1:Nbond
    Diff(n).pos = zeros(1,2); % Array of positions of configuration differences.
    Diff(n).val = zeros(1,2); % Array of difference values at those positions.
    Diff(n).num = 2; % Boson hopping terms alter two occupation numbers.
    Diff(n).type = 0; % For compatibility.
    Diff(n).sign = 1; % For compatibility.
end

for c = 1:CoOrd
    for n = 1:N
        m = Bonds(n,c); DInd = CoOrd*(n-1) + c; % Difference index.
        if m ~= 0 && m ~= n
            Diff(DInd).pos(1) = n;
            Diff(DInd).val(1) = 1; % Add boson to site n.
            Diff(DInd).pos(2) = m;
            Diff(DInd).val(2) = -1; % Remove boson from site m.
            OpMatEls(DInd,:) = sqrt(Cfg_vec(m)*(Cfg_vec(n)+1)) * ...
                ((Cfg_vec(n)+1) <= Cfg.Nmax) * exp(1i*n*q)/sqrt(N);
            % Second Boolean term comes into play for hardcore boson cases.
            
        elseif m == n % No bosons move, becomes total number counter.
            Diff(DInd).type = 1;
            Diff(DInd).pos = 1;
            Diff(DInd).val = 0;
            Diff(DInd).num = 0;
            OpMatEls(DInd,:) = sum(exp(1i*reshape(Cfg_vec,N,1)*q),1)/sqrt(N);
            break
        end
    end
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls(:,1)~=0;
OpMatEls = OpMatEls(ind,:);
Diff = Diff(ind); % Trim the structure array down.