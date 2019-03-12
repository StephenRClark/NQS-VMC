% Script for sampling pre-optimised Ansatz objects with the Fermi Hubbard
% model and any specified Operators.

%% Specify the basic information for the lattice.
N = 20; Dim = N;

%% Specify Graph object.
LVecs = 1;
GraphObj = HypCub(Dim,'PBC',LVecs);

%% Specify Hilbert object.
N_up = N/2; N_dn = N/2;
HilbertObj = Ferm(N,N_up,N_dn);

%% Specify initial Hamiltonian object.
dRT = 1; dRU = 0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT); HypCub(Dim,'PBC',dRU)];
% Specify the Hamiltonian energy terms e.g. hopping.
t = -1; U = [0 1/4 1/3 1/2 2/3 1 3/2 2 5/2 3 7/2 4 9/2 5 6 8];
HamiltonianObj = Hubbard(HilbertObj,HGraphs,t,U(1));

%% Initialise Sampler.

% Construct the Operator objects for sampling and place in a cell array. In
% this case, operators are double occupancy number, two-site density-density
% correlations and rigidity.
SampOperators = {OperatorDg(HilbertObj,HGraphs(2),@DbOc_Ferm_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@NiNj_Ferm_CfgVal),...
    OperatorDg(HilbertObj,HGraphs(2),@Rgdy_Ferm_CfgVal)};

SamplerObj = Sampler(HilbertObj,GraphObj,HamiltonianObj,SampOperators);

SamplerObj = SetNsamp(SamplerObj,10000);

% Initialise storage for the sampled quantities.
EneGS = zeros(numel(U),1); VarE = zeros(numel(U),1);
DbOc = zeros(numel(U),1); NiNj = zeros(numel(U),N);
Rgdy = zeros(numel(U),1);

% If sampling over many different Ansatz types, specify their names here. 
AnsStr = {'SDet-JastTI'};%{'SDet-Gutz','SDet-JastTI','SDet-NQSTI Alpha 1','SDet-NQSTISSJ Alpha 1',...
    % 'SDet-NQSTIDDJ Alpha 1','SDet-NQSTIDDJ-FW Alpha 2','SDet-NQSTISS Alpha 1','SDet-NQSTISS Alpha 2'}; % Ansatz identifier string.

% Set of line formats for plotting purposes.
LineFormat = {'-*'};%{'-^','-v','-x','-o','-*','-d','-s','-p'};

for a = 1:numel(AnsStr)
    for u = 1:numel(U)
        if u > 1
            % Create new Hamiltonian and cycle out old one.
            NewH = Hubbard(HilbertObj,HGraphs,t,U(u));
            SamplerObj = SetHamiltonian(SamplerObj,NewH);
        end
        tic;
        load(['FHM U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr{a} ' Logs.mat'],'AnsatzObj','EnIter','RunTime');
        AnsatzObj = PsiUpdate(AnsatzObj,dP);
        [EnAvg,EnSamp,EvalAvg,~] = EvalSample(SamplerObj,AnsatzObj);
        % Store sampled quantities.
        EneGS(u) = EnAvg/N; DbOc(u) = EvalAvg{1}; NiNj(u,:) = EvalAvg{2};
        Rgdy(u) = EvalAvg{3};
        VarE(u) = mean(((EnSamp - EnAvg)/EnAvg).^2); % Normalised local energy variance.
        disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
    end
    
    figure(1); hold on; plot(U,real(EneGS),LineFormat{a});
    
    figure(2); hold on; plot(U,DbOc,LineFormat{a});
    
    figure(3); hold on; plot(U,abs(Rgdy),LineFormat{a});
    
    save(['FHM N ' num2str(N) ' U Scan ' AnsStr{a} ' expectation values.mat'],'EneGS','DbOc','NiNj','VarE','Rgdy');
end

figure(1); legend(AnsStr);

figure(2); legend(AnsStr);

figure(3); legend(AnsStr);