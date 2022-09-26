N = 10; U = [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 25/4 20/3 15/2 35/4 10 12 16 20 30 50 100];

AnsStr = {'BECR-Jast-NQSNHTI-HD5 Alpha 1 WLR','BECR-Jast-NQSNHTI-HD5 Alpha 1 WSR','BECR-Jast-NNMB SeqOpt'};

LineSpec = {'-o','-s','-d'}; NStr = ' N-1 ';

for u = 1:numel(U)
    for a = 1:numel(AnsStr)
        figure(ceil(u/4)); subplot(2,2,1+mod(u-1,4)); hold on;
        load(['BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr{a} ' Logs.mat']);
        plot(real(EnIter));
        title(['U = ' num2str(U(u))]); 
    end
    legend(AnsStr);
end