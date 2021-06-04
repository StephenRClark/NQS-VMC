classdef Sampler
    % Sampler - contains methods and parameters for a Markov chain Monte
    % Carlo sample.
    %   Sampler can be fed into an Optimiser object, or used on its own
    %   by invoking MCMCSample to output sampled quantities.
    
    properties (SetAccess = protected)
        Nequil = 5000; % Number of equilibration 'burn in' runs.
        Nblock % Number of Markov chain moves made before sampling.
        Nsamp = 10000; % Number of samples to take.
        Hamiltonian % Contains the Hamiltonian operator for energy calculations.        
    end
    
    properties
        Operators % Non-Hamiltonian operators for sampling.
    end
    
    methods
        % Constructor for the Sampler object.
        function [obj] = Sampler(Hilbert,Hamiltonian,Operators)
            obj.Nblock = Hilbert.N; obj.Hamiltonian = Hamiltonian; % Requires Hamiltonian object.
            if nargin == 2 % If no Operators specified, assume none.
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
        
        % Reassign or add Operator of choice
        function [obj] = SetOperator(obj,NewOp,OpInd)
            if OpInd > numel(obj.Operators)
                disp(['New Operator index specified is greater than number of existing Operators - '...
                    'appending to current Operator list.'])
                obj.Operators{numel(obj.Operators)+1} = NewOp;
            else
                obj.Operators{OpInd} = NewOp;
            end
        end
    end
    
    methods       
        % MCMCSample: Markov Chain Monte Carlo sampling of an Ansatz using
        % Sampler object. Found in folder.
        [EnAvg,dLogpAvg,EvalAvg,MRate] = MCMCSample(obj,Ansatz)
        
        % EvalSample: Markov Chain Monte Carlo sampling for final
        % evaluation of observables. Found in folder.
        [EnAvg,EnSamp,EvalAvg,EvalSamp] = EvalSample(obj,Ansatz);
        
        % MultiChainSample: Multithreaded Markov Chain Monte Carlo sampling
        % for final evaluation of observables. Found in folder.
        [EnAvg,EnSamp,EvalAvg,EvalSamp] = MultiChainSample(obj,Ansatz,Ncore)
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Nsamp = obj.Nsamp; Properties.Nequil = obj.Nequil;
            Properties.Nblock = obj.Nblock; Properties.Hamiltonian = obj.Hamiltonian.PropertyList;
            Properties.Operator = cell(numel(obj.Operators),1);
            for o = 1:numel(obj.Operators)
                Properties.Operator{o} = obj.Operators{o}.PropertyList;
            end
        end
    end
    
end
