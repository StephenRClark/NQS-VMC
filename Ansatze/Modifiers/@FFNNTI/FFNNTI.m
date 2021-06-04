classdef FFNNTI < FFNN
    % FFNN - a Modifier subclass that modifies configuration amplitudes
    % using a feed forward neural network (FFNN) architecture of visible neurons
    % and layers of hidden spins.
    %   FFNNTI is a subclass of FFNN, Modifier is the overarching class.
    
    % ---------------------------------
    % Format for FFNN Modifier object:
    % - FFNN.Nv = number of "visible" spins.
    % - FFNN.Nh1 = number of "hidden" spins in first layer.
    % - FFNN.Nh2 = number of "hidden" spins in second layer.
    % - FFNN.Np = number of parameters in the ansatz = Alpha1 + Alpha1*Nv +
    % Alpha1*Nh2 + Nh2
    % - FFNN.b1 = (Nh1 x 1) vector - first hidden site bias.
    % - FFNN.b2 = (Nh2 x 1) vector - second hidden site bias.
    % - FFNN.W1 = (Nh1 x Nv) matrix - hidden1-visible coupling terms.
    % - FFNN.W2 = (Nh2 x Nh1) matrix - hidden2-hidden1 coupling terms.
    % - FFNN.W3 = (1 x Nh2) vector - average pooling layer matrix of ones.
    % - FFNN.Theta1 = (Nh1 x 1) vector - effective angles (hidden1).
    % - FFNN.Theta2 = (Nh2 x 1) vector - effective angles (hidden2).
    % - FFNN.h1 = (Nh1 x 1) vector - hidden unit values (first layer).
    % - FFNN.h2 = (Nh2 x 1) vector - hidden unit values (second layer).
    % Properties added with translation invariance:
    % - FFNN.Alpha1 = number of unique W1 parameter sets
    % - FFNN.b1r = (Alpha1 x 1) vector - unique 1st layer hidden biases.
    % - FFNN.W1r = (Alpha1 x Nv) matrix - unique 1st layer coupling terms.
    % - FFNN.W2r = (Nh2 x Alpha1) matrix - unique 1st-2nd layer coupling terms.
    % ---------------------------------
    % Format for Update is a vector of effective angles ThetaP
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Alpha1 x 1) for d/db1.
    % - (Nh2 x 1) for d/db2.
    % - (Alpha1*Nv x 1) for d/dW1.
    % - (Nh2*Alpha1 x 1) for d/dW2.
    % ---------------------------------
    
    
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier.
    end
    
    properties
        b1 = 0; % Hidden spin bias terms, Nh1 x 1 vector.
        b2 = 0; % Hidden spin bias terms, Nh2 x 1 vector.
        W1 = 0; % Hidden1-visible coupling terms, Nh1 x Nv matrix.
        W2 = 0; % Hidden2-hidden1 coupling terms, Nh2 x Nh1 matrix.
        W3 = 0; % Average pooling layer matrix, 1 x Nh2 vector.
        b1r = 0; % Reduced parameter set from imposing TI, Alpha1 x 1 vector.
        W1r = 0; % Reduced parameter set from imposing TI, Alpha1 x Nv matrix.
        W2r = 0; % Reduced parameter set from imposing TI, Nh2 x Alpha 1 matrix.
    end
    
    properties (SetAccess = protected)
        La = 2; % Number of layers
        Nv = 0; % Number of visible neurons.
        Nh1 = 0; % Number of hidden spins in the first layer.
        Nh2 = 0; % Number of hidden spins in the second layer.
        Np = 0; % Number of parameters.
        Alpha1 = 1; % Number of distinct 1st layer hidden unit sets.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        Theta1 = 0; % Layer 1 local effective angles, Nh1 x 1 vector.
        Theta2 = 0; % Layer 2 local effective angles, Nh2 x 1 vector.
        h1 = 0; % Layer 1 hidden unit values, Nh1 x 1 vector.
        h2 = 0; % Layer 2 hidden unit values, Nh2 x 1 vector.
        ParamCap = 5; % Parameter cap to mitigate effects of erroneous parameter changes.
        OptInds = zeros(162,1); % Individual parameter flags for variational purposes.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
        FFlag = 0; % Fermionic flag for incorporating Diff conversion in PsiRatio. Set to N if active.
        SFlag = 0; % Spin symmetry flag for fermionic case.
    end
    
    methods
        % Constructor for unconstrained FFNN with no symmetries:
        function obj = FFNNTI(Hilbert,Graph,Params,VFlag)
            % Graph necessary for FFNN subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified.
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            obj.FullCfg = @FullSpinCfg;
            obj.Nv = Hilbert.N;
            obj.Nh1 = Params.Nh1;
            obj.Nh2 = Params.Nh2;
            obj.Np = Params.Np;
            obj.Graph = Graph;
            obj = RandomInitPsiFFNNTI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateFFNNTI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from FFNN.
        
        % PrepPsi inherited from FFNN.
        
        % PsiGenerate inherited from FFNN.
        
        % PsiRatio inherited from FFNN.
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            [dLogp] = LogDerivFFNNTI(obj,Cfg);
        end
        
        % ParamList inherited from FFNN for the time being.
    end
end
