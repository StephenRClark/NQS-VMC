% Script for plotting the NQS weights of pre-optimised Ansatz objects
% saved in Logs. This only handles post-optimisation analysis - no Monte
% Carlo is performed here.

% Lattice and Hamiltonian parameters used for FHM 1D.

N = 20; k = -(1-2/N):2/N:1; dx = -(N/2 - 1):1:(N/2); 

U = [0 1/4 1/3 1/2 2/3 1 3/2 2 5/2 3 7/2 4 9/2 5 6 8];

% Due to the number of NQS variants, one needs to be careful when
% specifying the Ansatz objects being loaded.

Alpha = 2; AnsStr = ['SDet-NQSTISS Alpha ' num2str(Alpha)];

% Sequentially load and plot NQS W couplings in real space and momentum
% space for the Hamiltonian parameters specified.

for u = 1:numel(U)
    load(['FHM U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);

    % Figure 1 for energy optimisation curves to ensure that each Ansatz has been optimised properly. 
    figure(1); subplot(4,5,u); hold on;
    plot(real(EnIter)); title(['U = ' num2str(U(u))]);
    
    for a = 1:Alpha
        % Separate W couplings acting on up/down sites and individually
        % Fourier transform.
        WvU = AnsatzObj.Modifier{1}.W(1+(a-1)*N,1:N); 
        WvD = AnsatzObj.Modifier{1}.W(1+(a-1)*N,N+(1:N)); 
        WvUQ = circshift(real(fft(WvU)),(N/2)-1);
        WvDQ = circshift(real(fft(WvD)),(N/2)-1);
        [~,I] = max(abs(WvU)); % Locate largest element and center plot on that site.
        
        % Even numbered figures for real space W couplings.
        figure(2*a); subplot(4,4,u);        
        plot(dx,circshift(real(WvU),N/2-I),'-+'); hold on; plot(dx,circshift(real(WvD),N/2-I),'-x');
        title(['U = ' num2str(U(u)) ', Alpha = ' num2str(a)]); xlabel('dx'); ylabel('W(x)');
        
        % Odd numbered figures for k-space W couplings.
        figure(2*a+1); subplot(4,4,u);
        plot(k,WvUQ,'-o'); hold on;  plot(k,WvDQ,'-^');
        title(['U = ' num2str(U(u)) ', Alpha = ' num2str(a)]); xlabel('q'); ylabel('W(q)');
    end
end