L = 3; N = L^2; U = [1 16 23 32]; % [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48]; % 
SFlag = 0; % (numel(U)==18); 

NStr = ' N '; 
GSStr = ' exact ground state.mat'; Nmax = 3;
% GSStr = ' Nmax 2 exact ground state.mat'; Nmax = 2;

FidUFull = zeros(numel(U),1); FidUSub = zeros(numel(U),1); EnePsiU = zeros(numel(U),1); EneExact = zeros(numel(U),1);

DbHlPsi = zeros(N,N,numel(U)); DbHlGS = zeros(N,N,numel(U)); SSFPsi = zeros(L,L,numel(U)); SSFGS = zeros(L,L,numel(U));

HiHjPsi = zeros(N,N,numel(U)); HiHjGS = zeros(N,N,numel(U)); DiDjPsi = zeros(N,N,numel(U)); DiDjGS = zeros(N,N,numel(U));

DbDenPsi = zeros(numel(U),1); HlDenPsi = zeros(numel(U),1); DbDenGS = zeros(numel(U),1); HlDenGS = zeros(numel(U),1);

DenPsi = zeros(numel(U),1); DenSqPsi = zeros(numel(U),1); DenGS = zeros(numel(U),1); DenSqGS = zeros(numel(U),1);

NiNjPsi = zeros(L,L,numel(U)); NiNjGS = zeros(L,L,numel(U)); EneError = zeros(numel(U),1);  

OccFracGS = zeros(numel(U),4); OccFracPsi = zeros(numel(U),4);

SPDMPsi = cell(numel(U),1); SPDMGS = cell(numel(U),1); SPDMQPsi = cell(numel(U),1); SPDMQGS = cell(numel(U),1); 

SFFracPsi = zeros(numel(U),1); SFFracGS = zeros(numel(U),1);

Q = (0:1:(L-1))'*2*pi/L;

% AnsStr = 'BECR-Jast FT';
% AnsStr = 'BECR-Jast-NNMB FT';
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 2'; 
% AnsStr = 'BECR-NQSNHTI-HD4-NNMB Alpha 2';
% AnsStr = 'BECR-NQSNHTI-HD4 Alpha 2';
% AnsStr = 'BECR-NQSMHTI-HD4-LR Alpha 1';
% AnsStr = 'BECR-Gutz-NNMB FT';
% AnsStr = 'BECR-NQS-HD4-ABNH-DHSR Alpha 2';
% AnsStr = 'BECR-NQSTI-I Alpha 1';
% AnsStr = 'BECR-NQSOHTI-VD4 Alpha 1';
% AnsStr = 'BECR-CPSATI-V4H3 Alpha 1 HTest FT';
% AnsStr = 'BECR-CPSXTI-V4H2 Alpha 1 FT';
% AnsStr = 'BECR-NQSS1TI-V4-S2 Alpha 1 FT'; 
% AnsStr = 'BECR-NQSSXTI-V4-S2 Alpha 1 FT';

% AnsStr = 'BECR-Jast-NNMB-NQSTI-I Alpha 1 FT';
% AnsStr = 'BECR-Jast-NNMB-NQSTI Alpha 1 FT';
% AnsStr = 'BECR-Jast-NNMB-NQSMHTI-HD4 Alpha 1 FT';
% AnsStr = 'BECR-Jast-NNMB-NQSNHTI-HD4 Alpha 1 FT';
% AnsStr = 'BECR-Jast-NNMB-NQSOHTI-VD4 Alpha 1 FT';

for u = 1:numel(U)
    load(['BHM 2D' NStr num2str(N) ' U ' num2str(U(u)) GSStr]);
    load(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj','EnIter','RunTime');
    Psi = AnsatzObj.Reference.PsiGenerate(basis); PsiMod = cell(numel(AnsatzObj.Modifier),1);
    for m = 1:numel(PsiMod)
        PsiMod{m} = AnsatzObj.Modifier{m}.PsiGenerate(basis);
        if sum(abs(PsiMod{m}.^2))>0
            Psi = Psi .* PsiMod{m};
        else
            disp(['Error in calculating exact wavefunction for Modifier ' ...
                num2str(m) ', U = ' num2str(U(u)) '.']);
        end
    end
    % Calculate fidelity with ground state.
    Psi = Psi / sqrt(sum(abs(Psi).^2)); FidUFull(u) = abs(gs'*Psi)^2;
    % Find Nmax = 2 subspace and find fidelity within subspace.
    SubInds = find(sum(basis>2,2)==0);
    SubGS = gs(SubInds); SubGS = SubGS / sqrt(sum(abs(SubGS).^2));
    SubPsi = Psi(SubInds); SubPsi = SubPsi / sqrt(sum(abs(SubPsi).^2));
    FidUSub(u) = abs(SubGS'*SubPsi)^2;
    % Calculate full energy.
    EnePsi = Psi'*H*Psi; EnePsiU(u) = EnePsi; EneError(u) = abs((EnePsi - gs_en)/gs_en); EneExact(u) = gs_en;
    % Calculate the part of Psi orthogonal to ground state.
    PsiRemnant = Psi - (gs'*Psi)*gs; SubPsiRemnant = SubPsi - (SubGS'*SubPsi)*SubGS;
    % Evaluate static structure factors and doublon-holon correlations
    for n = 0:Nmax
        OccFracGS(u,n+1) = gs'*ocfr{n+1}*gs;
        OccFracPsi(u,n+1) = Psi'*ocfr{n+1}*Psi;
    end    
    DbDenPsi(u) = OccFracPsi(u,3); HlDenPsi(u) = OccFracPsi(u,1); 
    DbDenGS(u) = OccFracGS(u,3); HlDenGS(u) = OccFracGS(u,1);
    for i = 1:L
        for j = 1:L
            SSFPsi(i,j,u) = Psi'*Nq{i+(j-1)*L}*Psi; SSFGS(i,j,u) = gs'*Nq{i+(j-1)*L}*gs;
            for i2 = 1:L
                for j2 = 1:L
                    n = i+(j-1)*L; m = i2+(j2-1)*L; I = 1+mod(i2-i,L); J = 1+mod(j2-j,L);
                    NiNjPsi(I,J,u) = NiNjPsi(I,J,u) + Psi'*ninj{n,m}*Psi/N; 
                    NiNjGS(I,J,u) = NiNjGS(I,J,u) + gs'*ninj{n,m}*gs/N;
                end
            end
        end
    end
    SPDMGS{u} = zeros(N); SPDMPsi{u} = zeros(N);
    for n = 1:N
        DenPsi(u) = DenPsi(u) + Psi'*num{n}*Psi/N; DenGS(u) = DenGS(u) + gs'*num{n}*gs/N;
        DenSqPsi(u) = DenSqPsi(u) + Psi'*(num{n}^2)*Psi/N; DenSqGS(u) = DenSqGS(u) + gs'*(num{n}^2)*gs/N;
        for m = 1:N
            DbHlPsi(n,m,u) = Psi'*dbhl{n,m}*Psi; DbHlGS(n,m,u) = gs'*dbhl{n,m}*gs;
            HiHjPsi(n,m,u) = Psi'*hihj{n,m}*Psi; HiHjGS(n,m,u) = gs'*hihj{n,m}*gs;
            DiDjPsi(n,m,u) = Psi'*didj{n,m}*Psi; DiDjGS(n,m,u) = gs'*didj{n,m}*gs;
            SPDMPsi{u}(n,m) = Psi'*hop{m+(n-1)*N}*Psi; SPDMGS{u}(n,m) = gs'*hop{m+(n-1)*N}*gs;
        end
    end
    SPDMQPsi{u} = fft2(SPDMPsi{u}); SPDMQGS{u} = fft2(SPDMGS{u});
    SFFracPsi(u) = real(SPDMQPsi{u}(1,1))/(sum(den)*N); SFFracGS(u) = real(SPDMQGS{u}(1,1))/(sum(den)*N);
    % Add calculated quantities to Logs.
    BiBj = SPDMPsi{u}; BpBq = SPDMQPsi{u}; ninj = NiNjPsi(:,:,u); DbHl = DbHlPsi(:,:,u);
    VarN = DenSqPsi(u) - DenPsi(u)^2; SSF = SSFPsi(:,:,u); Fid = FidUFull(u); 
    Infid = 1- Fid; OccDist = OccFracPsi(u,:); HiHj = HiHjPsi(:,:,u); DiDj = DiDjPsi(:,:,u);
    HlDen = HlDenPsi(u); DbDen = DbDenPsi(u); SFFrac = SFFracPsi(u); Params = AnsatzObj.ParamList;
    save(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','RunTime','EnIter','Psi','EnePsi','VarN','BiBj','BpBq','ninj',...
        'DbDen','HlDen','DbHl','SSF','SFFrac','Fid','Infid','Params','OccDist');
end

VarNPsi = DenSqPsi - (DenPsi.^2); VarNGS = DenSqGS - (DenGS.^2);

figure(1); 
plot(U,EneExact,'-x'); hold on; plot(U,EnePsiU,'-o');
title(['Variational versus exact ground state energy for ' AnsStr]);
xlabel('U/t'); ylabel('E/N'); 

figure(2); 
plot(U,FidUFull,'-o'); hold on;
title(['Fidelity with ground state for ' AnsStr]);
xlabel('U/t'); ylabel('F');

figure(3); 
plot(U,EneError,'-o'); hold on;
title(['Relative energy error for ' AnsStr]);
xlabel('U/t'); ylabel('(E_{\Psi} - E_{0} )/ E_{0}');

figure(4); 
plot(U,sqrt(VarNGS),'-x'); hold on; plot(U,sqrt(VarNPsi),'-o');
title(['Deviation in occupation for ' AnsStr]);
xlabel('U/t'); ylabel('\sigma');

figure(5);
plot(U,SFFracGS,'-x'); hold on; plot(U,SFFracPsi,'-o');
title(['Condensate fraction for ' AnsStr]);
xlabel('U/t'); ylabel('\rho_{0} / N'); 

if SFlag == 1
save(['BHM 2D' NStr num2str(N) ' ' AnsStr ' exact expectation values.mat'],...
    'FidUFull','FidUSub','EnePsiU','EneExact','EneError','VarNPsi','VarNGS','SFFracPsi','SFFracGS',...
    'NiNjPsi','NiNjGS','DenPsi','DenGS','DenSqPsi','DenSqGS','SPDMPsi','SPDMQPsi','SPDMGS','SPDMQGS',...
    'DbHlPsi','SSFPsi','DbHlGS','SSFGS','DbDenPsi','DbDenGS','HlDenPsi','HlDenGS','OccFracGS','OccFracPsi');
end