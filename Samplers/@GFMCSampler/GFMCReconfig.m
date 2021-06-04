% --- Continuous time Green's Function Monte Carlo reconfiguration function ---

function [CfgNew,AnsNew] = GFMCReconfig(CfgW,AnsW,Weights)

Weights = cumsum(Weights)/sum(Weights); % Conversion of weights to probability distribution.
Nw = numel(CfgW); AnsNew = cell(Nw,1); CfgNew = cell(Nw,1);
RW = ((0:(Nw-1)) + rand)/Nw; Inds = zeros(Nw,1); % Correlated random numbers for reselection.
for w = 1:Nw % Reassign configurations and prepare Ansatz objects.
    Inds(w) = sum(Weights<RW(w))+1; 
    CfgNew{w} = CfgW{Inds(w)};
    AnsNew{w} = PrepPsi(AnsW{w},CfgNew{w});
end

end