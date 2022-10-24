% Script for collating logs for energy comparisons and plotting.
N = 20; JzV = [0 0.25 0.5 1 2 4]; SysStr = 'XXZ 1D';
% Specify the Ansatz type to collate logs for. Ensure folder with results
% is added to MATLAB path.
AnsStr = 'Plus-NQSTI Alpha 1';

% Initialise storage for the sampled quantities.
EneGSJ = zeros(numel(JzV),1); VarEJ = zeros(numel(JzV),1);
SzSzJ = cell(numel(JzV),1); SiSjJ = cell(numel(JzV),1);

for j = 1:numel(JzV)
    load([SysStr ' Jz ' num2str(JzV(j)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    EneGSJ(j) = EneGS; VarEJ(j) = VarE; SiSjJ{j} = SiSj; SzSzJ{j} = SzSz;
end

RunDate = date;
% Save collated observables in separate file.
save([SysStr ' Jz Scan N ' num2str(N) ' ' AnsStr ' expectation values.mat'],...
    'EneGSJ','VarEJ','SiSjJ','SzSzJ','RunDate');

% Plot sampled energy.
figure(1); hold on; plot(JzV,real(EneGSJ),'-o');