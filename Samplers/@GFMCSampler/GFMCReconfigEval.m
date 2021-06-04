function [EnAvg,OpAvg,CfgNew,AnsNew,WAvg,Inds] = GFMCReconfigEval(EneLocVec,OpValCell,CfgW,AnsW,Weights)

WAvg = mean(Weights); EnAvg = sum(EneLocVec .* Weights)/sum(Weights); % Average energy after normalised weighting.
Nw = numel(Weights); OpAvg = cell(1,size(OpValCell,2)); Inds = zeros(Nw,1);
AnsNew = cell(Nw,1); CfgNew = cell(Nw,1);
for w = 1:Nw % Average and repackage any operator values.
    for o = 1:size(OpValCell,2)
        if w == 1
            OpAvg{o} = OpValCell{w,o} * Weights(w)/sum(Weights);
        else
            OpAvg{o} = OpAvg{o} + OpValCell{w,o} * Weights(w)/sum(Weights);
        end
    end
end
Weights = cumsum(Weights)/sum(Weights); % Conversion of weights to probability distribution.
RW = ((0:(Nw-1)) + rand)/Nw; % Correlated random numbers for reselection.
for w = 1:Nw % Reassign configurations and prepare Ansatz objects.
    Inds(w) = sum(Weights<RW(w))+1; CfgNew{w} = CfgW{Inds(w)};
    AnsNew{w} = PrepPsi(AnsW{w},CfgNew{w});
end
