U = [1 16 23 32];

u_plot = [1 2 3 4];

N = 100; Nmax = 4; npref = 600;

np_jast = N; np_j_hdc = N+1;
np_nqs_i_a1 = N+1; np_nqs_a_a1 = N+2; np_nqs_ab_a1 = N+3;
np_nqs_i_a2 = 2*N+2; np_nqs_a_a2 = 2*N+3; np_nqs_ab_a2 = 2*N+5;
np_nqs_hd = 2*N+3; np_nqs_ab_hdc = N+4; np_nqs_u = 1+Nmax+N*(Nmax+1);

np_ans_bec = [np_jast; np_j_hdc; np_nqs_i_a1; np_nqs_a_a1; np_nqs_ab_a1; np_nqs_i_a2;...
    np_nqs_a_a2; np_nqs_ab_a2; np_nqs_hd; np_nqs_ab_hdc; np_nqs_u];

AnsStrBEC ={'BECR-Jast';'BECR-Jast-HDC';'BECR-NQS-I Alpha 1';'BECR-NQS-A Alpha 1';...
    'BECR-NQS-AB-HDim5 Alpha 1';'BECR-NQS-I Alpha 2';'BECR-NQS-A Alpha 2';...
    'BECR-NQS-AB-HDim5 Alpha 2';'BECR-NQS-HD-HDim5 Alpha 1';'BECR-NQS-AB-HDim5-HDC';...
    'BECR-NQS-U-VDim5 Alpha 1';};

LegStrBEC = {'Jastrow','Jast-HDC','I \alpha = 1','A \alpha = 1','AB \alpha = 1',...
    'I \alpha = 2','A \alpha = 2','AB \alpha = 2','HD \alpha = 1','AB-HDC','U \alpha = 1'};

AnsMarkBEC = {'bs','rs','bv','cv','gs','bd','cd','gd','kv','ks','k*'};

Ind_J = 1; Ind_JHD = 2; Ind_U = 11;

n_ans_bec = numel(np_ans_bec); err_ans_bec = cell(n_ans_bec,1); ene_ans_bec = cell(n_ans_bec,1);

for n = 1:n_ans_bec
    ene_ans_bec{n} = zeros(numel(u_plot),1);
    for u = 1:numel(u_plot)
        load(['BHM 2D U ' num2str(U(u_plot(u))) ' N ' num2str(N) ' ' AnsStrBEC{n} ' Logs.mat']);
        ene_ans_bec{n}(u) = EneGS;
        err_ans_bec{n}(u) = sqrt(VarE);
    end
end
save(['BHM 2D N ' num2str(N) ' BEC reference ansatz comparisons.mat'],'n_ans_bec',...
    'np_ans_bec','ene_ans_bec','err_ans_bec','AnsStrBEC','LegStrBEC','AnsMarkBEC','U');

AnsStrJHD = {'BECR-JHD-NQS-I';'BECR-JHD-NQS-R';'BECR-JHD-NQS-AB-HDim5';...
    'BECR-JHD-NQS-HD-HDim5';'BECR-JHD-NQS-U-VDim5'};

LegStrJHD = {'JHDC'; 'JHDC-NQS-R'; 'JHDC-NQS-I'; 'JHDC-NQS-AB';'JHDC-NQS-HD'; 'JHDC-NQS-U';'Reference'};

np_ans_jhd = [np_nqs_i_a1; np_nqs_i_a1; np_nqs_ab_a1; np_nqs_hd; np_nqs_u];

n_ans_jhd = numel(np_ans_jhd); err_ans_jhd = cell(n_ans_jhd,1); ene_ans_jhd = cell(n_ans_jhd,1);

AnsMarkJHD = {'rs','c^','bv','gs','kd','k*'};

for n = 1:n_ans_jhd
    ene_ans_jhd{n} = zeros(numel(u_plot),1);
    for u = 1:numel(u_plot)
        load(['BHM 2D U ' num2str(U(u_plot(u))) ' N ' num2str(N) ' ' AnsStrJHD{n} ' Logs.mat']);
        ene_ans_jhd{n}(u) = EneGS;
        err_ans_jhd{n}(u) = sqrt(VarE);
    end
end
save(['BHM 2D N ' num2str(N) ' JHD reference ansatz comparisons.mat'],'n_ans_jhd',...
    'np_ans_jhd','ene_ans_jhd','err_ans_jhd','AnsStrJHD','LegStrJHD','AnsMarkJHD','U');

for p = 1:numel(u_plot)
    u = u_plot(p); figure; hold on;
    for a = 1:n_ans_bec
        errorbar(np_ans_bec(a),ene_ans_bec{a}(p),err_ans_bec{a}(p),AnsMarkBEC{a});
    end
    % Title and axis labelling.
    legend(LegStrJHD); xlabel('Number of parameters'); ylabel('Energy per site');
    title(['U = ' num2str(U(u)) ' energy vs parameters']);
    % Energy reference lines
    plot([0 npref],ene_ans_bec{Ind_J}(p)*ones(1,2),'b--');
    plot([0 npref],ene_ans_bec{Ind_JHD}(p)*ones(1,2),'r--');
    plot([0 npref],ene_ans_bec{Ind_U}(p)*ones(1,2),'k--');
end

for p = 1:numel(u_plot)
    u = u_plot(p); figure; hold on;
    for a = 1:n_ans
        errorbar(np_ans_jhd(a),ene_ans_jhd{a}(p),err_ans_jhd{a}(p),AnsMarkJHD{a});
    end
    % Energy reference lines
    plot([0 np_ref],ene_ans_bec{Ind_JHD}(p)*ones(1,2),'r--');
    % Title and axis labelling.
    legend(LegStrJHD); xlabel('Number of parameters'); ylabel('Energy per site');
    title(['U = ' num2str(U(u)) ' energy vs parameters']);
end