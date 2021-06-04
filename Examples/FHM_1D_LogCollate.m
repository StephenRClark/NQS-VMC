% Script for collating logs for energy comparisons and plotting.
N = 20; U = [2 4 6 8 10]; Nmax = 4; SysStr = 'FHM 1D';
% Specify the Ansatz type to collate logs for. Ensure folder with results
% is added to MATLAB path.
AnsStr = 'BECR-Jast';

% Initialise storage for the sampled quantities.
EneGSU = zeros(numel(U),1); VarEU = zeros(numel(U),1);
DbOcU = cell(numel(U),1); NiNjU = cell(numel(U),1);
RgdyU = cell(numel(U),1);

for u = 1:numel(U)
    load([SysStr ' U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    EneGSU(u) = EneGS; VarEU(u) = VarE; NiNjU{u} = NiNj;
    DbOcU{u} = DbOc; RgdyU{u} = Rgdy;
end

RunDate = date;
% Save collated observables in separate file.
save([SysStr ' U Scan N ' num2str(N) ' ' AnsStr ' expectation values.mat'],...
    'EneGSU','VarEU','NiNjU','DbOcU','RgdyU','RunDate');

% Plot sampled energy.
figure(1); hold on; plot(U,real(EneGSU),'-o');