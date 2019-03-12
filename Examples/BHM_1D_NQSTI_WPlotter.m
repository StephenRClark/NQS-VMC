% Script for plotting the NQS weights of pre-optimised Ansatz objects
% saved in Logs. This only handles post-optimisation analysis - no Monte
% Carlo is performed here.

% Lattice and Hamiltonian parameters used for BHM 1D.

N = 20; k = -(1-2/N):2/N:1; dx = -(N/2-1):1:(N/2); Alpha = 1;

U = [2 2.1 2.2 2.25 2.3 2.35 2.375 2.4 2.425 2.45 2.475 2.5 2.525 2.55 2.575 2.6 2.7 2.8 2.9 3];

% Due to the number of NQS variants, one needs to be careful when
% specifying the Ansatz objects being loaded.

AnsStr = 'BECR-NQSNHTI-CJ-HD5 Alpha 1';

% Sequentially load and plot NQS W couplings in real space and momentum
% space for the Hamiltonian parameters specified.
for u = 1:numel(U)
    load(['BHM U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    
    % Figure 1 for energy optimisation curves to ensure that each Ansatz has been optimised properly. 
    figure(1); subplot(4,5,u); hold on;
    plot(real(EnIter)); title(['U = ' num2str(U(u))]);
    for a = 1:Alpha
        Wv = AnsatzObj.Modifier{1}.Wv(a,:);
        WvQ = circshift(fft(Wv),(N/2)-1);
        
        % Even numbered figures for real space W couplings.
        figure(2*a); subplot(4,5,u);
        plot(real(Wv),'-+'); % hold on; plot(imag(Wv),'-x');
        title(['U = ' num2str(U(u)) ', Alpha = ' num2str(a)]); xlabel('dx'); ylabel('W(x)');
        
        % Odd numbered figures for k-space W couplings.
        figure(2*a+1); subplot(4,5,u);
        plot(k,real(WvQ),'-o'); % hold on; plot(k,imag(WvQ),'-^');
        title(['U = ' num2str(U(u)) ', Alpha = ' num2str(a)]); xlabel('q'); ylabel('W(q)');
    end
end