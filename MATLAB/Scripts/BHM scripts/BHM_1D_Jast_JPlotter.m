% Script for plotting the Jastrow factors of pre-optimised Ansatz objects
% saved in Logs. This only handles post-optimisation analysis - no Monte
% Carlo is performed here.

% Lattice and Hamiltonian parameters used for BHM 1D.

N = 80; q = -(1-2/N):2/N:1; dx = -(N/2-1):1:(N/2);

U = [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 25/4 20/3 15/2 35/4 10 12 16 20 30 50 100];

AnsStr = 'BECR-Jast FT'; NStr = ' N ';

% load(['BHM 1D' NStr num2str(N) ' U Scan ' AnsStr ' MMC expectation values.mat']);

% Sequentially load and plot Jastrow factors in real space and momentum
% space for the Hamiltonian parameters specified.
for u = 1:numel(U)
    load(['BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat']);
    Js = AnsatzObj.Modifier{1}.Js(:,1);
    JsQ = circshift(real(fft(Js)),(N/2)-1) .*(q.'.^2); % q^2 multiplication optional.
    
    % Figure 1 for energy optimisation curves to ensure that each Ansatz has been optimised properly. 
    figure(1); subplot(5,5,u);
    plot(real(EnIter)); title(['U = ' num2str(U(u))]);
    
    % Figure 2 for real space Jastrow.
    figure(2); subplot(5,5,u); 
    plot(dx,circshift(Js,(N/2)-1),'-+');    
    title(['U = ' num2str(U(u))]); xlabel('dx'); ylabel('J(dx)');
    
    % Figure 3 for k-space Jastrow multiplied by q^2.
    figure(3); subplot(5,5,u);
    plot(q(q>0),JsQ(q>0),'-o'); title(['U = ' num2str(U(u))]); xlabel('q'); ylabel('J(q)'); 
    
    % Figure 4 for k-space dependence of SSF Nq/q
    % figure(4); subplot(5,5,u);
    % plot(q(q>1e-5),StatSF(u,2:(round(N/2)+1))./q(q>1e-5),'-o'); title(['U = ' num2str(U(u))]); xlabel('q'); ylabel('N(q)'); 
end
