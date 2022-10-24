U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

Nvec = [16 36 64 100]; Nmax = 4; % 

for n = 1:numel(Nvec)
    
    N = Nvec(n); 
    
%     AnsStr = 'BECR-Gutz'; AnsMark = 'c-o';
%     AnsStr = 'BECR-Gutz-NNMB'; AnsMark = 'b-s';
%     AnsStr = 'BECR-Jast'; AnsMark = 'r-o';
%     AnsStr = 'BECR-Jast-NNMB'; AnsMark = 'r-s';
%     AnsStr = 'BECR-NQSTI-I Alpha 1'; AnsMark = 'b-^';
%     AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 1'; AnsMark = 'r-^';
%     AnsStr = 'BECR-NQSSHTI-HD5-bB Alpha 1'; AnsMark = 'g-^';
%     AnsStr = 'BECR-NQSNHTI-HD5-bB Alpha 1'; AnsMark = 'c-^';
%     AnsStr = 'BECR-NQSTI-I Alpha 2'; AnsMark = 'b-d';
%     AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 2'; AnsMark = 'r-d';
%     AnsStr = 'BECR-NQSSHTI-HD5 Alpha 2'; AnsMark = 'g-d';
%     AnsStr = 'BECR-NQSNHTI-HD5 Alpha 2'; AnsMark = 'c-d';
%     AnsStr = 'BECR-NQSMHTI-HD5-LR Alpha 1'; AnsMark = 'k-v';
    AnsStr = 'BECR-NQSNHTI-HD5-NNMB Alpha 2'; AnsMark = 'g-s';
%     AnsStr = 'BECR-NQSOHTI-VD5 Alpha 1'; AnsMark = 'k-*';
%     AnsStr = 'BECR-CPSATI-V5H2 Alpha 1'; AnsMark = 'g-*';
%     AnsStr = 'BECR-CPSXTI-V5H2 Alpha 1'; AnsMark = 'c-*';
    
    % Initialise storage for the sampled quantities.
    EneGSU = zeros(numel(U),1); VarEU = zeros(numel(U),1);
    VarNU = zeros(numel(U),1); NiNjU = cell(numel(U),1);
    DbHlU = cell(numel(U),1); OcFrU = zeros(numel(U),Nmax+1);
    BiBjU = cell(numel(U),1);
    
    for u = 1:numel(U)
        load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
        EneGSU(u) = EneGS; VarNU(u) = VarN; VarEU(u) = VarE; NiNjU{u} = NiNj;
        DbHlU{u} = DbHl; OcFrU(u,:) = OcFr.'; BiBjU{u} = BiBj;
    end
    
    figure(1); hold on; plot(U,real(EneGSU),AnsMark);
    
    RunDate = date;
    
    save(['BHM 2D U Scan N ' num2str(N) ' ' AnsStr ' MMC expectation values.mat'],...
        'EneGSU','VarNU','NiNjU','VarEU','DbHlU','OcFrU','BiBjU','RunDate');
end