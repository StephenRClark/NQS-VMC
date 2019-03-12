% Script for plotting the Jastrow factors of pre-optimised Ansatz objects
% saved in Logs. This only handles post-optimisation analysis - no Monte
% Carlo is performed here.

% Lattice and Hamiltonian parameters used for BHM 1D.

N = 20; k = -(1-2/N):2/N:1; dx = -(N/2-1):1:(N/2);

U = [2 2.1 2.2 2.25 2.3 2.35 2.375 2.4 2.425 2.45 2.475 2.5 2.525 2.55 2.575 2.6 2.7 2.8 2.9 3];

% Sequentially load and plot Jastrow factors in real space and momentum
% space for the Hamiltonian parameters specified.
for u = 1:numel(U)
    load(['BHM U ' num2str(U(u)) ' N ' num2str(N) ' BECR-JastTI Logs.mat']);
    Js = AnsatzObj.Modifier{1}.Js(:,1);
    JsQ = circshift(real(fft(Js)),(N/2)-1);% .*(k.'.^2); % k^2 multiplication optional.
    
    % Figure 1 for energy optimisation curves to ensure that each Ansatz has been optimised properly. 
    figure(1); subplot(4,5,u);
    plot(real(EnIter)); title(['U = ' num2str(U(u))]);
    
    % Figure 2 for real space Jastrow.
    figure(2); subplot(4,5,u); 
    plot(dx,circshift(Js,(N/2)-1),'-+');    
    title(['U = ' num2str(U(u))]); xlabel('dx'); ylabel('J(dx)');
    
    % Figure 3 for k-space Jastrow multiplied by k^2.
    figure(3); subplot(4,5,u);
    plot(k(k>0),JsQ(k>0),'-o'); title(['U = ' num2str(U(u))]); xlabel('q'); ylabel('J(q)'); 
end
