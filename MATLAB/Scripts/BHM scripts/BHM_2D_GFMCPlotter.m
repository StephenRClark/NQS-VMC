U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

L = 4; N = L^2; AnsStr = 'BECR-Gutz-NNMB'; EneGSEx = zeros(numel(U),1);

EneGSU = zeros(numel(U),1); VarNU = zeros(numel(U),1); VarNFW = zeros(numel(U),1);

BiBjU = cell(numel(U),1); BiBjFW = cell(numel(U),1); NiNjU = cell(numel(U),1); NiNjFW = cell(numel(U),1);

OcFrU = zeros(numel(U),5); OcFrFW = zeros(numel(U),5); DbHlU = cell(numel(U),1); DbHlFW = cell(numel(U),1);

CondFracU = zeros(numel(U),1); CondFracFW = zeros(numel(U),1);

for u = 1:numel(U)
    load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' GFMC expectation values.mat']);
    EneGSU(u) = EneGS/N; VarNU(u) = EvalAvgMA{1}; VarNFW(u) = EvalAvgFW{1};
    % SPDM has double counted in the off-diagonal elements - rectifying
    % here.
    SPDMM = (EvalAvgMA{5} + diag(diag(EvalAvgMA{5})))/2;
    SPDMF = (EvalAvgFW{5} + diag(diag(EvalAvgFW{5})))/2;
    SPDMMQ = fft2(SPDMM); SPDMFQ = fft2(SPDMF);
    CondFracU(u) = SPDMMQ(1)/(N^2); CondFracFW(u) = SPDMFQ(1)/(N^2);
    BiBjU{u} = SPDMM; BiBjFW{u} = SPDMF;
    % Want to plot NiNj - 1 with VarN as diagonal elements. Existing
    % correlations are scaled by factor N, so making corrections here.
    NiNjM = EvalAvgMA{2}/N; NiNjF = EvalAvgFW{2}/N; 
    NiNjM(2:end) = NiNjM(2:end)-1; NiNjF(2:end) = NiNjF(2:end)-1;
    NiNjU{u} = reshape(NiNjM,L,L); NiNjFW{u} = reshape(NiNjF,L,L);
    % Need to normalise doublon-holon correlations according to the
    % doublon and holon densities, which is recorded in the occupation
    % number distribution OcFr. Properly scaled, should have values around
    % 1 for mostly uncorrelated sites.
    OcFrU(u,:) = EvalAvgMA{4}.'; OcFrFW(u,:) = EvalAvgFW{4}.';
    DHM = OcFrU(u,1)*OcFrU(u,3); DHF = OcFrFW(u,1)*OcFrFW(u,3);
    DbHlU{u} = reshape(EvalAvgMA{3}/DHM,L,L); DbHlFW{u} = reshape(EvalAvgFW{3}/DHF,L,L);
    % Exact diagonalisation results for 3x3.
%     load(['BHM 2D N ' num2str(N) ' U ' num2str(U(u)) ' exact ground state.mat']);
%     EneGSEx(u) = gs_en/N;
end

figure(1); plot(U,EneGSU); hold on; % plot(U,EneGSEx,'-x');

EneErr = abs((EneGSU - EneGSEx)./EneGSEx);

% figure(2); semilogy(U,EneErr);

save(['BHM 2D U Scan N ' num2str(N) ' ' AnsStr ' GFMC expectation values.mat'],...
    'EneGSU','VarNU','VarNFW','BiBjU','BiBjFW','NiNjU','NiNjFW','OcFrU','OcFrFW',...
    'DbHlU','DbHlFW','CondFracU','CondFracFW');