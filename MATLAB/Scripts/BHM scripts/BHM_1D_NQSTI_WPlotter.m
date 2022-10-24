% Script for plotting the NQS weights of pre-optimised Ansatz objects
% saved in Logs. This only handles post-optimisation analysis - no Monte
% Carlo is performed here.

% Lattice and Hamiltonian parameters used for BHM 1D.

N = 8; k = -(1-2/N):2/N:1; dx = -(N/2-1):1:(N/2); Alpha = 2;

U = [0 1 2 3 4 6 8 12 16 24];

% Due to the number of NQS variants, one needs to be careful when
% specifying the Ansatz objects being loaded.

AnsStr = 'BECR-NQSTI-C Alpha 2';

% Sequentially load and plot NQS W couplings in real space and momentum
% space for the Hamiltonian parameters specified.
for u = 1:numel(U)
    load(['BHM 1D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    
    % Figure 1 for energy optimisation curves to ensure that each Ansatz has been optimised properly. 
    figure(1); subplot(2,5,u); hold on;
    plot(real(EnIter(EnIter<0))); title(['U = ' num2str(U(u))]);
    for a = 1:Alpha
        Wv = AnsatzObj.Modifier{1}.W(1+N*(a-1),:);
        WvQ = circshift(fft(Wv),(N/2)-1);
        
        % Even numbered figures for real space W couplings.
        figure(a+1); subplot(2,5,u);
        [~,I] = max(abs(Wv));
        plot(dx,circshift(real(Wv),N/2-I),'-+'); hold on; plot(dx,circshift(imag(Wv),N/2-I),'-x');
        title(['U = ' num2str(U(u)) ', Alpha = ' num2str(a)]); xlabel('dx'); ylabel('W(x)');
        
%         % Odd numbered figures for k-space W couplings.
%         figure(2*a+1); subplot(2,5,u);
%         plot(k,real(WvQ),'-o'); hold on; plot(k,imag(WvQ),'-^');
%         title(['U = ' num2str(U(u)) ', Alpha = ' num2str(a)]); xlabel('q'); ylabel('W(q)');
    end
end