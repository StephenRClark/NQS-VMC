N = 10; U = [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 ...
    25/4 20/3 15/2 35/4 10 12 16 20 30 50 100]; dx = (1:N) - 1;

FidUFull = zeros(numel(U),1); FidUSub = zeros(numel(U),1); EnePsiU = zeros(numel(U),1); EneGSEx = zeros(numel(U),1);

EneError = zeros(numel(U),1); load('dbhl_op_n10_nb11.mat'); load('ssf_op_n10_nb11.mat');

DbHlPsi = zeros(numel(U),N); DbHlGS = zeros(numel(U),N); SSFPsi = zeros(numel(U),N); SSFGS = zeros(numel(U),N);

DbDenPsi = zeros(numel(U),1); DbDenGS = zeros(numel(U),1); HlDenPsi = zeros(numel(U),1); HlDenGS = zeros(numel(U),1);

Q = (0:1:(N-1))'*2*pi/N;

AnsStr = 'BECR-NQSTI-I Alpha 2'; ULegend = cell(numel(U),1); NStr = ' N+1 ';

for u = 1:numel(U)
    ULegend{u} = ['U = ' num2str(U(u))];
    load(['BHM 1D' NStr num2str(N) ' U ' num2str(U(u)) ' exact ground state.mat']);
    load(['BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat']);
    Psi = AnsatzObj.Reference.PsiGenerate(basis); PsiMod = cell(numel(AnsatzObj.Modifier),1);
    for m = 1:numel(PsiMod)
        PsiMod{m} = AnsatzObj.Modifier{m}.PsiGenerate(basis);
        if sum(abs(PsiMod{m}))>0
            Psi = Psi .* PsiMod{m};
        else
            disp(['Error in calculating exact wavefunction for Modifier ' num2str(m) ', U = ' num2str(U(u)) '.']);
        end
    end
    % Calculate fidelity with ground state.
    Psi = Psi / sqrt(sum(abs(Psi).^2)); FidUFull(u) = abs(gs'*Psi)^2;
    % Divide by phase of highest weighted state.
    [~,I] = max(abs(Psi.^2));
    Psi = Psi/sign(Psi(I)); gs = gs * -sign(gs(I));
    % Find Nmax = 2 subspace and find fidelity within subspace.
    SubInds = find(sum(basis>2,2)==0);
    SubGS = gs(SubInds); SubGS = SubGS / sqrt(sum(abs(SubGS).^2));
    SubPsi = Psi(SubInds); SubPsi = SubPsi / sqrt(sum(abs(SubPsi).^2));
    FidUSub(u) = abs(SubGS'*SubPsi)^2;
    % Calculate full energy.
    EnePsi = Psi'*H*Psi; EnePsiU(u) = EnePsi; EneError(u) = abs((EnePsi - gs_en)/gs_en); EneGSEx(u) = gs_en;
    % Calculate doublon and holon densities.
    DbDenPsi(u) = real(Psi'*DbDen*Psi)/N; DbDenGS(u) = gs'*DbDen*gs/N;
    HlDenPsi(u) = real(Psi'*HlDen*Psi)/N; HlDenGS(u) = gs'*HlDen*gs/N;
    % Evaluate static structure factors and doublon-holon correlations
    for n = 1:N
        DbHlPsi(u,n) = real(Psi'*dbhl{n,1}*Psi) / (DbDenPsi(u)*HlDenPsi(u)); 
        DbHlGS(u,n) = gs'*dbhl{n,1}*gs  / (DbDenGS(u)*HlDenGS(u));
        SSFPsi(u,n) = real(Psi'*Nq{n}*Psi); SSFGS(u,n) = gs'*Nq{n}*gs;
    end
    % Plot wavefunctions and expectation values.
%     figure(4); subplot(5,5,u);
%     plot(Q,real(SSFPsi(u,:)).'./Q,'-o'); hold on; plot(Q,SSFGS(u,:).'./Q,'-x');
%     title(['U = ' num2str(U(u)) ' SSF']);
%     legend({'Ansatz','Exact'});
end

figure(1); plot(U,real(EnePsiU),'-o'); hold on; plot(U,EneGSEx,'-x');
title(['Variational versus exact ground state energy for ' AnsStr]);

figure(2); plot(U,FidUFull,'-o'); hold on; plot(U,FidUSub,'-x');
title(['Fidelity with ground state for ' AnsStr]);

figure(3); plot(U,EneError); hold on;
title(['Relative energy error for ' AnsStr]);
% 
% figure(5); 
% subplot(1,2,1); plot(dx,DbHlGS,'-x'); title('Doublon holon correlations, exact'); 
% xlabel('dx'); ylabel('DiHJ/DH'); legend(ULegend); % axis([0 (N-1) 0 0.55]);
% subplot(1,2,2); plot(dx,DbHlPsi,'-o'); title('Doublon holon correlations, ansatz'); 
% xlabel('dx'); ylabel('DiHJ/DH'); legend(ULegend); % axis([0 (N-1) 0 0.55]);

save(['BHM 1D' NStr num2str(N) ' ' AnsStr ' exact expectation values.mat'],...
    'FidUFull','FidUSub','EnePsiU','EneError','DbHlPsi','SSFPsi','DbDenPsi','DbDenGS','HlDenPsi','HlDenGS');