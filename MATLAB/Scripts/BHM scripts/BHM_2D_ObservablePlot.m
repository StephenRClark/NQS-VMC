U = [0 2 4 6 8 9 10 11 12 14 16 18 20 22 24 28 32 40];

AnsStr = 'Jastrow + MBC';

figure(1); 
plot(U(2:18),EneExact(2:18),'--'); hold on; plot(U(2:18),EnePsiU(2:18),'-o');
title('Energy per site');
xlabel('U/t'); ylabel('E/N'); 
legend('Exact',AnsStr);

figure(2); 
semilogy(U(2:18),1-FidUFull(2:18),'-o'); hold on;
title(['Fidelity difference with ground state, ' AnsStr])
xlabel('U/t'); ylabel('1 - F');

figure(3); 
semilogy(U(2:18),EneError(2:18),'-o'); hold on; 
title(['Relative energy error, ' AnsStr])
xlabel('U/t'); ylabel('(E_{\Psi} - E_{0} )/ E_{0}');

figure(4); 
plot(U(1:18),sqrt(VarNGS(1:18)),'--'); hold on; plot(U(1:18),sqrt(VarNPsi(1:18)));
title('On-site occupation deviation');
xlabel('U/t'); ylabel('\sigma'); axis([0 Inf 0 1]);
legend('Exact',AnsStr);

figure(5);
plot(U(1:18),SFFracGS(1:18),'--'); hold on; plot(U(1:18),SFFracPsi(1:18),'-o');
title('Condensate fraction');
xlabel('U/t'); ylabel('\rho_{0} / N'); 
legend('Exact',AnsStr);