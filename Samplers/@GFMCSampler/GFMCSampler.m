classdef GFMCSampler
    % GFMCSampler - contains methods and parameters for a Green's function
    % continuous time Monte Carlo sample.
    %   GFMCSampler is used entirely to evaluate energies and operators,
    %   and does not contain provisions to be used with Optimisers.
    
    properties (SetAccess = protected)
        Nequil = 500; % Number of equilibration reconfigurations.        
        Nwalk = 10; % Number of walkers during GFMC.
        Tbranch = 0.2; % Time step length before branching.
        Tequil = 0.2; % Time step during equilibration phase.
        Nsamp = 100; % Number of samples to take.
        Pmax = 10; % Maximum projection of cumulant weights.
        Ncore = 1; % Number of cores / threads available.
        Hamiltonian % Contains the Hamiltonian operator for energy calculations.
    end
    
    properties
        Operators % Non-Hamiltonian operators for sampling.
    end
    
    methods
        % Constructor for the Sampler object.
        function [obj] = GFMCSampler(Ncore,HamiltonianObj,Operators)
            obj.Hamiltonian = HamiltonianObj; % Requires Hamiltonian object.
            obj.Ncore = Ncore; obj.Nwalk = Ncore; obj.Pmax = obj.Nwalk;
            % If multiple cores used, a 1 to 1 ratio of cores to walkers is assumed.
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
        
        % Reassign branching parameters if defaults are unsuitable.
        function [obj] = SetTbranch(obj,Tnew)
            obj.Tbranch = Tnew;
        end
        
        % Reassign equilbration time step if default is unsuitable.
        function [obj] = SetTequil(obj,Tnew)
            obj.Tequil = Tnew;
        end
              
        % Reassign maximum projection index if default is unsuitable.
        function [obj] = SetNbranch(obj,Nnew)
            obj.Nbranch = Nnew;
        end
        
        % Reassign maximum projection index if default is unsuitable.
        function [obj] = SetProjection(obj,Pnew)
            obj.Pmax = Pnew;
        end
        
        % Reassign walker and core numbers if defaults are unsuitable.
        function [obj] = SetNwalk(obj,NWnew,NCnew)
            if nargin == 3
                obj.Ncore = NCnew;
            end
            obj.Nwalk = NWnew;
        end
        
        % Reassign Hamiltonian if initial input is unsuitable.
        function [obj] = SetHamiltonian(obj,NewH)
            obj.Hamiltonian = NewH;
        end
        
        % Reassign or add Operator of choice
        function [obj] = SetOperator(obj,NewOp,OpInd)
            if OpInd > numel(obj.Operator)
                disp(['New Operator index specified is greater than number of existing Operators - '...
                    'appending to current Operator list.'])
                obj.Operator{numel(obj.Operator)+1} = NewOp;
            else
                obj.Operator{OpInd} = NewOp;
            end
        end
    end
    
    methods
        % GFMCSample: Continuous time Green's Function Monte Carlo sampling
        % of an Ansatz using GFMCSampler object. Found in folder.
        [EvalAvgP,WAvgP] = GFMCSampleCT(obj,AnsatzObj)
        % If no Operators being evaluated, EvalAvgP returns the energy.
        [EvalAvgP,WAvgP] = GFMCSampleDT(obj,AnsatzObj)
        
        % GFMCChain: CTGFMC sampling of an Ansatz, outputting the
        % configurations used and a record of reconfiguration changes.
        % Found in folder.
        [EnAvgP,WAvgP,CfgsP,RcfgIndsP] = GFMCChain(obj,AnsatzObj)
        
        % GFMCChainSample: Given a list of configurations and
        % reconfiguration changes, will sample the operators assigned with
        % both mixed averages and forward walking values.
        [EvalAvgMA,EvalAvgFW] = GFMCChainSample(obj,AnsatzObj,CfgsP,RcfgIndsP,WAvgP)
        
    end
    
end
