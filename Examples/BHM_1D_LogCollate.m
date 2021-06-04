% Script for collating logs for energy comparisons and plotting.
N = 20; U = [2 4 6 8 10]; Nmax = 4; SysStr = 'BHM 1D';
% Specify the Ansatz type to collate logs for. Ensure folder with results
% is added to MATLAB path.
AnsStr = 'BECR-Jast';

% Initialise storage for the sampled quantities.
EneGSU = zeros(numel(U),1); VarEU = zeros(numel(U),1);
VarNU = zeros(numel(U),1); NiNjU = cell(numel(U),1);
DbHlU = cell(numel(U),1); OcFrU = zeros(numel(U),Nmax+1);
BiBjU = cell(numel(U),1);

for u = 1:numel(U)
    load([SysStr ' U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    EneGSU(u) = EneGS; VarNU(u) = VarN; VarEU(u) = VarE; NiNjU{u} = NiNj;
    DbHlU{u} = DbHl; OcFrU(u,:) = OcFr.'; BiBjU{u} = BiBj;
end

RunDate = date;
% Save collated observables in separate file.
save([SysStr ' U Scan N ' num2str(N) ' ' AnsStr ' expectation values.mat'],...
    'EneGSU','VarNU','NiNjU','VarEU','DbHlU','OcFrU','BiBjU','RunDate');

% Plot sampled energy.
figure(1); hold on; plot(U,real(EneGSU),'-o');