% Script for plotting the Jastrow factors of pre-optimised Ansatz objects
% saved in Logs. This only handles post-optimisation analysis - no Monte
% Carlo is performed here.

% Lattice and Hamiltonian parameters used for FHM 1D.

N = 20; k = -(1-2/N):2/N:1; dx = -(N/2-1):1:(N/2);

U = [0 1/4 1/3 1/2 2/3 1 3/2 2 5/2 3 7/2 4 9/2 5 6 8];

% Sequentially load and plot Jastrow factors in real space and momentum
% space for the Hamiltonian parameters specified.
for u = 1:numel(U)
    load(['FHM U ' num2str(U(u)) ' N ' num2str(N) ' SDet-JastTI Logs.mat']);
    Js = AnsatzObj.Modifier{1}.Js(:,1);
    JsQ = circshift(real(fft(Js)),(N/2)-1); % .*(k.'.^2); % k^2 multiplication optional.
    
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
