% --- General fermionic Pfaffian wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioPfaf(PfafObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for fermionic Pfaffian Reference:
% - Pfaf.Nf = (1 x 2) vector - number of up/down fermions respectively.
% - Pfaf.PairMat = (2N x 2N) matrix - contains all pairing terms.
% - Pfaf.PfI = (Nf x Nf) matrix - inverse of reduced PfFull matrix.
% - Pfaf.PfG = (2N x Nf) matrix - matrix used for ratio calculations.
% - Pfaf.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
% - Pfaf.Np = number of variational parameters associated with Pfaf Reference.
% Pfaf properties used in variational version:
% - Pfaf.PfV = (2N x 2N) array - logs which variational parameters make up the elements of Pfaf.PairMat.
% - Pfaf.PfVR = (Nf x Nf) array - reduced matrix constructed from PfV.
% - Pfaf.PfVar = (Np x 1) vector - variational parameters in PfFull.
% ---------------------------------
% Format for Update is a struct containing updates for PfI, PfG and FermLoc.
% ---------------------------------

% Make local copies to reduce notation in code below.
N = numel(PfafObj.FermLoc)/2; % Number of sites.
Nf = sum(PfafObj.Nf); % Number of fermions.
% Compute the fermionic part of the ratio
FermLocP = PfafObj.FermLoc; PfMat = PfafObj.PairMat;
PfGUpdate = zeros(2*N,Nf); PfIUpdate = zeros(Nf); Ratio = 1;
if Diff.num ~= 0
    PfGMatP = PfafObj.PfG; PfIMatP = PfafObj.PfI;
    FermInd = zeros(Diff.num,1); % Storage for fermion indices that are being shifted
    FermDest = zeros(Diff.num,1); % Storage for fermion destinations (in FermLoc vector)
    
    for d = 1:Diff.num % Primitive implementation updates with each move
        if Diff.type == 0 % Type 0 - spin / fermion position exchange.
            % Diffs arising from a viable swap will list one pair of
            % positions, and Diff.val will be a [2x2] matrix, of which the
            % first row details the up fermion move.
            % First entry in Diff.pos will be initial site of up fermion.
            FermInd(d) = FermLocP(Diff.pos(d)+(mod(d-1,2)*N)); % Find index of fermion on site.
            FermDest(d) = Diff.pos(1+mod(d,Diff.num)) + (mod(d-1,2)*N); % Fermion destination
            FermLocP(FermLocP==FermInd(d)) = 0; % Remove fermion from starting point in FermLocP.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            
        elseif Diff.type == 1 % Type 1 - up fermion move.
            FermInd(d) = FermLocP(Diff.pos(2*d-1)); % Independent fermion moves have numel(Diff.pos) = 2*Diff.num.
            FermDest(d) = Diff.pos(2*d); % Destination of up fermion.
            FermLocP(FermLocP==FermInd(d)) = 0; % Remove fermion from starting point in FermLocP.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            
        elseif Diff.type == -1 % Type -1 - down fermion move.
            FermInd(d) = FermLocP(Diff.pos(2*d-1)+N); % Independent fermion moves have numel(Diff.pos) = 2*Diff.num.
            FermDest(d) = Diff.pos(2*d) + N; % Destination of down fermion.
            FermLocP(FermLocP==FermInd(d)) = 0; % Remove fermion from starting point in FermLocP.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            
        elseif Diff.type == 2 % Type 2 - pair move.
            % Diffs arising from a viable pair move will list one pair of
            % positions, and Diff.val will be a [2x2] matrix, of which the
            % first row details the up fermion move.
            % First entry in Diff.pos will be initial site of pair.
            FermInd(d) = FermLocP(Diff.pos(1) + N*mod(d-1,2)); % First site is site being vacated.
            FermDest(d) = Diff.pos(2) + N*mod(d-1,2); % Destination of pair.
            FermLocP(FermLocP==FermInd(d)) = 0; % Remove fermion from starting point in FermLocP.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
        end
        
        Ratio = Ratio * PfGMatP(FermDest(d),FermInd(d));
        
        % Update PfGMatP and PfIMatP
        
        PfInds = zeros(Nf,1);
        for f = 1:Nf
            PfInds(f) = find(FermLocP==f);
        end
        
        WVec = - PfGMatP(FermDest(d),:); WVec(FermInd(d)) = WVec(FermInd(d)) + 1;
        PfIUpdateP = ((PfIMatP(:,FermInd(d)) * WVec) + ...
            (WVec.' * PfIMatP(FermInd(d),:))) ...
            / PfGMatP(FermDest(d),FermInd(d));
        PfGUpdateP = ((PfGMatP(:,FermInd(d)) * WVec) + ...
            ((PfMat(:,FermDest(d)) + (PfGMatP * PfMat(FermDest(d),PfInds).')) * PfIMatP(FermInd(d),:)))...
            / PfGMatP(FermDest(d),FermInd(d));
        
        PfIUpdate = PfIUpdate + PfIUpdateP; PfGUpdate = PfGUpdate + PfGUpdateP;
        PfGMatP = PfGMatP + PfGUpdateP; PfIMatP = PfIMatP + PfIUpdateP;
        % Update working version of PfG and PfI within loop - overall
        % update gets stored in PfGUpdate / PfIUpdate
    end
end

% Store updates in struct for compatibility with current PsiRatio
% implementations
Update.FermLoc = FermLocP; Update.PfI = PfIUpdate; Update.PfG = PfGUpdate;