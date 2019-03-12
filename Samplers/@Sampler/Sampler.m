classdef Sampler
    % Sampler - contains methods and parameters for a Markov chain Monte
    % Carlo sample.
    %   Sampler can be fed into an Optimiser object, or used on its own
    %   by invoking MCMCSample to output sampled quantities.
    
    properties
        N % Number of sites.
        Nequil = 5000; % Number of equilibration 'burn in' runs.
        Nblock % Number of Markov chain moves made before sampling.
        Nsamp = 10000; % Number of samples to take.
        Hamiltonian % Contains the Hamiltonian operator for energy calculations.
        Operators % Non-Hamiltonian operators for sampling.
        Hilbert % Contains configuration information for the Markov chain.
        Graph % Contains information about the lattice.
    end
    
    methods
        % Constructor for the Sampler object.
        function [obj] = Sampler(Hilbert,Graph,Hamiltonian,Operators)
            obj.Hilbert = Hilbert; obj.Graph = Graph;
            obj.N = Hilbert.N; obj.Nblock = Hilbert.N;
            obj.Hamiltonian = Hamiltonian; % Requires Hamiltonian object.
            if nargin == 3 % If no Operators specified, assume none.
                obj.Operators = {};
            else
                obj.Operators = Operators; % Requires cell array of Operator objects.
            end
        end
        
        % Reassign Nequil if default is unsuitable.
        function [obj] = SetNequil(obj,Nnew)
            obj.Nequil = Nnew;            
        end
        
        % Reassign Nsamp if default is unsuitable.
        function [obj] = SetNsamp(obj,Nnew)
            obj.Nsamp = Nnew;
        end
        
        % Reassign Nblock if default is unsuitable.
        function [obj] = SetNblock(obj,Nnew)
            obj.Nblock = Nnew;            
        end
        
        % Reassign Hamiltonian if initial input is unsuitable.
        function [obj] = SetHamiltonian(obj,NewH)
            obj.Hamiltonian = NewH;
        end
    end
    
    methods       
        % Markov Chain Monte Carlo sampling of an Ansatz using Sampler object.
        [EnAvg,dLogpAvg,EvalAvg,MRate] = MCMCSample(obj,Ansatz)
        
        % Markov Chain Monte Carlo sampling for final evaluation of
        % observables.
        [EnAvg,EnSamp,EvalAvg,EvalSamp] = EvalSample(obj,Ansatz);
    end
    
end
