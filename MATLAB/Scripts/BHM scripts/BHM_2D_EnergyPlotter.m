L = 6; N = L^2; U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

FigNum = 1; urange = 5:13;

AnsStr{1} = 'BECR-Gutz'; AnsStr{2} = 'BECR-Jast'; 
AnsStr{3} = 'BECR-Gutz-NNMB'; AnsStr{4} = 'BECR-Jast-NNMB';
LineSpec{1} = 'g-o'; LineSpec{2} = 'b-o'; LineSpec{3} = 'g-d'; LineSpec{4} = 'b-d';
LegEntry{1} = 'Gutzwiller'; LegEntry{2} = 'Jastrow'; 
LegEntry{3} = 'Gutz w/ MBC'; LegEntry{4} = 'Jast w/ MBC';

% AnsStr{1} = 'BECR-Jast'; AnsStr{2} = 'BECR-Jast-NNMB'; 
% AnsStr{3} = 'BECR-NQSSHTI-HD2-AW Alpha 1'; AnsStr{4} = 'BECR-NQSSHTI-HD5-bB Alpha 1';
% LineSpec{1} = 'b-o'; LineSpec{2} = 'b-d'; LineSpec{3} = 'k-+'; LineSpec{4} = 'g-s';
% LegEntry{1} = 'Jastrow'; LegEntry{2} = 'Jast w/ MBC'; 
% LegEntry{3} = 'N_{A} d_{h} = 2'; LegEntry{4} = 'N_{AB} d_{h} = d_{v}';

% AnsStr{1} = 'BECR-Jast-NNMB'; AnsStr{2} = 'BECR-NQSSHTI-HD5-bB Alpha 1'; 
% AnsStr{3} = 'BECR-NQSMHV-HD5 Alpha 1'; AnsStr{4} = 'BECR-NQSMHTI-HD5-LR Alpha 1';
% LineSpec{1} = 'b-d'; LineSpec{2} = 'g-s'; LineSpec{3} = 'm-x'; LineSpec{4} = 'm-d';
% LegEntry{1} = 'Jast w/ MBC'; LegEntry{2} = 'N_{AB} d_{h} = d_{v}'; 
% LegEntry{3} = 'N_{DH} hybrid'; LegEntry{4} = 'N_{DH} long range';

% AnsStr{1} = 'BECR-Jast-NNMB'; AnsStr{2} = 'BECR-NQSSHTI-HD5-bB Alpha 1'; 
% AnsStr{3} = 'BECR-NQSSHTI-HD2-AW Alpha 2'; AnsStr{4} = 'BECR-NQSSHTI-HD5 Alpha 2';
% LineSpec{1} = 'b-d'; LineSpec{2} = 'g-s'; LineSpec{3} = 'k-x'; LineSpec{4} = 'm-d';
% LegEntry{1} = 'Jast w/ MBC'; LegEntry{2} = 'N_{AB} d_{h} = d_{v} \alpha = 1'; 
% LegEntry{3} = 'N_{AB} d_{h} = 2 \alpha = 2'; LegEntry{4} = 'N_{AB} d_{h} = d_{v} \alpha = 2';;

for a = 1:numel(AnsStr)
    load(['BHM 2D U Scan N ' num2str(N) ' ' AnsStr{a} ' MMC expectation values.mat']);
    figure(FigNum); 
    subplot(1,2,1); hold on; plot([0 U],[-4; EneGSU],LineSpec{a});
    subplot(1,2,2); hold on; plot(U(urange),EneGSU(urange),LineSpec{a});
end

subplot(1,2,1); xlabel('U /t'); ylabel('E / N'); legend(LegEntry); title('Whole interaction range');
subplot(1,2,2); xlabel('U /t'); ylabel('E / N'); legend(LegEntry); title('Critical region');