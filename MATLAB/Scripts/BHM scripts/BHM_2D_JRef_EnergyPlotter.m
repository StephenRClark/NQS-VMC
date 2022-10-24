U = [1 16 23 32]; u_plot = [1 2 3 4]; % Alter u_plot if you only want to see specific interaction strengths

N = 9; % Choose from [9 16 36 64 100].

load(['BHM 2D N ' num2str(N) ' ansatz comparisons.mat']);

for p = 1:numel(u_plot)
    u = u_plot(p); figure(p); hold on;
    if N == 9
        for a = 1:n_ans
            scatter(p_nqs_ans(a),ene_ans{a}(p),AnsMark{a});
        end
        % Energy reference lines.
        plot([0 p_nqs_oh],ene_ans{1}(p)*ones(1,2),'r--');
        % Title and axis labelling.
        legend(LegStr); xlabel('Number of parameters'); ylabel('Energy per site');
        title(['U = ' num2str(U(u)) ' energy vs parameters']);
        % Infidelity plot.
        figure(p+numel(u_plot)); ax = gca; ax.YScale = 'log'; hold on;
        for a = 1:n_ans
            scatter(p_nqs_ans(a),infid_ans{a}(p),AnsMark{a});
        end       
        % Infidelity reference lines
        plot([0 p_nqs_oh],infid_ans{1}(p)*ones(1,2),'r--');
        % Title and axis labelling. 
        title(['U = ' num2str(U(u)) ' infidelity vs parameters']);
        legend(LegStr); xlabel('Number of parameters'); ylabel('Infidelity');
    else
        for a = 1:n_ans
            errorbar(p_nqs_ans(a),ene_ans{a}(p),err_ans{a}(p),AnsMark{a});
        end        
        % Energy reference lines
        plot([0 p_nqs_oh],ene_ans{1}(p)*ones(1,2),'r--');
        % Title and axis labelling.
        legend(LegStr); xlabel('Number of parameters'); ylabel('Energy per site');
        title(['U = ' num2str(U(u)) ' energy vs parameters']);
    end
end