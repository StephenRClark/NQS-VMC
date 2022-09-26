% Script for constructing and optimising Ansatz objects with number-hidden
% NQS Modifiers and Bose condensate References for the Bose Hubbard model
% on a 2D lattice with periodic boundary conditions.

%% Add the folders containing class definitions to MatLab's path.
addpath(genpath('Variational Monte Carlo'));

%% Use setup script for all the necessary objects.
BHM_2D_OptSetup_100_Q;
% Hilbert, Hamiltonian, Sampler, StochasticReconfig and Graph objects created in BHM_2D_OptSetup_100_Q

%% Initialise starting Ansatz.

% Initialise Reference and Modifier(s) to slot into Ansatz.

% Chosen Reference: BECR
% Necessary fields: Hilbert, Graph, Params.SPH

% Can construct a single particle Hamiltonian (hopping only) using Graph2Array.
CArr = Graph2Array(GraphObj,0); RefParams.SPH = -CArr;

Ref = BECR(HilbertObj,GraphObj,RefParams);

% Noise terms - nmag is a small perturbation such that a = a_0 + nmag*rand
% nphs adds random phase such that a = a_0*exp(1i*nphs*rand)
ModParams.nmag = 0.01; ModParams.nphs = 0; 

ModParams.Nh = N; % System size N = 100 initialised in OptSetup
% Set base values of NQS parameters. Can fine-grain and specify A and B as well in NQSNH case
ModParams.a = -0.01; ModParams.b = 0.1; ModParams.W = -0.01;

% Chosen Modifier: translation invariant NQS with number-like hidden units and square biases, Alpha = 1.
NQSObj = NQSNHTI(HilbertObj,GraphObj,ModParams,1);

% Initialise Ansatz object. Mods need to be passed as a cell list.

AnsatzObj = Ansatz(Ref,{JastObj;NNMBObj},HilbertObj); 

AnsStr = 'BECR-Jast'; % Ansatz name string for easy identification later.

urange = 1:numel(U); % Values of U in a vector initialised in OptSetup.

DirStr = [AnsStr NStr num2str(L) 'x' num2str(L)]; mkdir(DirStr); 

for u = urange
    BHM_2D_3Opt % This chunk handles the 3 stage optimisation of the ansatz. Optimisation hyperparameters in OptSetup
end

exit