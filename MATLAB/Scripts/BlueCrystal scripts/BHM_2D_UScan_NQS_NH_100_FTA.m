% Script for constructing and optimising Ansatz objects with number-hidden
% NQS Modifiers and Bose condensate References for the Bose Hubbard model
% on a 2D lattice with periodic boundary conditions.

%% Add the folders containing class definitions to MatLab's path.
addpath(genpath('Variational Monte Carlo'));

%% Use setup script for all the necessary objects.
BHM_2D_FTSetup_100_Q;

%% Initialise starting Ansatz.

AnsStr0 = 'BECR-NQSNHTI-HD5-bB Alpha 1'; % Ansatz name string for easy identification later.
AnsStr = [AnsStr0 ' FT'];

urange = 1:numel(U);

DirStr = [AnsStr NStr num2str(L) 'x' num2str(L)]; mkdir(DirStr); 

for u = urange
    load(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr0 ' Logs.mat'],'AnsatzObj');
    BHM_2D_2Opt
end

exit