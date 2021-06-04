classdef TimeVariational < Optimiser
    % TimeVariational - an Operator subclass that is used to calculate the
    % parameter changes that give the trajectory closest to the true time
    % evolution of the wavefunction in Hilbert space using the time
    % dependent variational principle.
    %    Optimiser is overarching class. TimeVariational features single
    %    and multithread options, as well as a choice of integration
    %    scheme.
    
    properties
        Int % Integration scheme - defaults to naive step integration.
        TStep % Time step size, which will be normalised to the characteristic time scale of the system.
        TScale % Time step rescaling to be dimensionless.
    end
    
    properties (Hidden) % Generally do not need to modify these, as these are just numerical checks.
        STol = 1e-30; % Tolerance in overlap matrix element magnitude - anything below is ignored.
    end
    
    methods
        % General constructor for subclass - additional arguments in ETols and SRParams.
        function obj = TimeVariational(Npass,Ncore,TStep,TScale,Int)
            obj@Optimiser(Npass,Ncore);
            if nargin == 3
                error('No dimensionless time step scaling specified.')
            elseif nargin == 4 % Assume default step integration if not specified.
                obj.TStep = TStep; obj.TScale = TScale; obj.Int = 'Def';
            else
                obj.TStep = TStep; obj.TScale = TScale; obj.Int = Int;
            end            
            if strcmp(obj.Int,'Def') == 0 || strcmp(obj.Int,'RK4')
                error('Integration identifier unknown. Current options are "Def" (default) or "RK4" (Runge-Kutta, 4th order).')
            end
        end
        
        % Reassign time step and time scale if initial values are
        % unsuitable.
        function [obj] = SetTStep(obj,TStep_p,TScale_p)
            if nargin == 3 % Only change TScale if requested, otherwise leave alone.
                obj.TScale = TScale_p;
            end
            obj.TStep = TStep_p;
        end
        
        % Reassign integration style if initial choice is unsuitable.
        function [obj] = SetIntegration(obj,IntP)
            if strcmp(IntP,'Def') == 0 || strcmp(IntP,'RK4')
                error('Integration identifier unknown. Current options are "Def" (default) or "RK4" (Runge-Kutta, 4th order).')
            end
            obj.Int = IntP;
        end
    end
    
    methods
        function [Ansatz,EvalIter] = Optimise(obj,Sampler,Ansatz)
            if strcmp(obj.Int,'Def')
                [Ansatz,EvalIter] = TimeEvolveAnsatz(obj,Sampler,Ansatz);
            elseif strcmp(obj.Int,'RK4')
                [Ansatz,EvalIter] = TimeEvolveAnsatz_RK4(obj,Sampler,Ansatz);
            end
        end
    end
    
end

