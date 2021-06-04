% --- Continuous time Green's Function Monte Carlo Markov chain generating function ---

function [EvalAvgMA,EvalAvgFW] = GFMCChainSample(GFMCObj,AnsatzObj,CfgsP,RcfgIndsP,WAvgP)

% Enacts importance sampling CTGFMC using provided Ansatz object for the
% guiding wavefunction.
Ncore = GFMCObj.Ncore; Nwalk = size(RcfgIndsP,2); Nsamp = GFMCObj.Nsamp;
Pmax = GFMCObj.Pmax; Operators = GFMCObj.Operators; OpNum = numel(Operators);

% Set up storage for evaluated quantities, giving both biased mixed
% averages and (mostly) unbiased forward walking values.
EvalAvgMA = cell(OpNum,1); EvalAvgFW = cell(OpNum,1);
for o = 1:OpNum
    Op0 = 0*Operators{o}.GraphSample(CfgsP{1},0,0,AnsatzObj);
    EvalAvgMA{o} = Op0; EvalAvgFW{o} = 0;
end
OpArray = cell(OpNum,Nwalk); % Set up array to get around parfor oddities.
for o = 1:OpNum
    for nw = 1:Nwalk
        OpArray{o,nw} = Operators{o};
    end
end

% Set up vectors of raw and projected weights.
RawWeights = WAvgP((1:Nsamp) + Pmax); ProjWeights = ones(Nsamp,1);
for p = 1:Pmax
    ProjWeights = ProjWeights .* WAvgP((1:Nsamp)+Pmax-p);
end
RawWeights = RawWeights/sum(RawWeights);
ProjWeights = ProjWeights/sum(ProjWeights);

% Clone Ansatz object for individual use by each worker.
AnsMA = cell(Nwalk,1);
for nw = 1:Nwalk
    AnsMA{nw} = AnsatzObj;
end
AnsFW = AnsMA;

% For the forward walking technique, need to trace which configurations
% appeared at the 'projection iteration' using RcfgIndsP.
PCfgs = cell(Nsamp,Nwalk);
% Need to step back by m = Pmax/2 iterations.
for n = 1:Nsamp
    Inds = RcfgIndsP(n+Pmax,:);
    % Trace back the indices step by step.
    for p = 1:(Pmax/2)
        Inds = RcfgIndsP(n+Pmax-p,Inds);
    end
    PCfgs(n,:) = CfgsP(n+(Pmax/2),Inds);
end

if Ncore > 1
    Pool = parpool(Ncore); % Save for deletion at the end.
end

for n = 1:Nsamp
    CfgMA = CfgsP(n+Pmax,:); CfgFW = PCfgs(n,:);    
    for nw = 1:Nwalk
        AnsMA{nw} = PrepPsi(AnsMA{nw},CfgMA{nw});
        AnsFW{nw} = PrepPsi(AnsFW{nw},CfgFW{nw});
    end
    for o = 1:OpNum
        OpP = OpArray(o,:);
        EvalMA = cell(Nwalk,1); EvalFW = cell(Nwalk,1);
        % parfor (nw = 1:Nwalk, Ncore)
        for nw = 1:Nwalk    
            % Mixed average naive sampling.            
            EvalMA{nw} = OpP{nw}.GraphSample(CfgMA{nw},0,0,AnsMA{nw})/Nwalk;
            % Forward walking projected sampling - accurate for diagonal
            % operators, discard for non-diagonal operators.
            EvalFW{nw} = OpP{nw}.GraphSample(CfgFW{nw},0,0,AnsFW{nw})/Nwalk;
        end
        for nw = 1:Nwalk
            EvalAvgMA{o} = EvalAvgMA{o} + EvalMA{nw}*RawWeights(n);
            EvalAvgFW{o} = EvalAvgFW{o} + EvalFW{nw}*ProjWeights(n);
        end        
    end
end

if Ncore > 1
    delete(Pool);
end

end
