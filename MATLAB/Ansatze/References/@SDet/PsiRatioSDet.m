% --- General fermionic determinant wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioSDet(SDetObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for SDet Reference:
% - SDet.Nf = (1 x 2) vector - number of up/down fermions respectively.
% - SDet.Orbitals = (2N x 2N) matrix - contains all available single particle orbitals.
% - SDet.OrbMat = (2N x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
% - SDet.DetMat = (2N x Nf) matrix - elements are used for determinants in PsiRatio.
% - SDet.FermLoc = (2N x 1) vector - details locations of fermions by index for sign tracking purposes.
% - SDet.Np = number of variational parameters associated with SDet Reference.
% SDet properties used in variational version:
% - SDet.CArr = (2N x 2N x Np) array - connectivity array for the reference Hamiltonian.
% - SDet.WArr = (2N x 2N x Np) array - transformed connectivity array for the reference Hamiltonian.
% - SDet.EnFac = (2N x 2N) matrix - elements are used in LogDeriv function.
% - SDet.HVar = (Np x 1) vector - variational parameters in the reference Hamiltonian used.
% ---------------------------------
% Format for Update is a struct containing updates for DetMat and FermLoc.
% ---------------------------------

% Make local copies to reduce notation in code below.
N = numel(SDetObj.FermLoc)/2; % Number of sites.
Nf = sum(SDetObj.Nf); % Number of fermions.

% Compute the fermionic part of the ratio.
FermLocP = SDetObj.FermLoc;
DetMatUpdate = zeros(2*N,Nf);

% Diff.val is structured as [U1 U2; D1 D2] if both an up and down fermion
% move in the same configuration change, otherwise Diff.val is a [1x2]
% vector with up/down specified by Diff.type.

if Diff.num ~= 0
    BFe = zeros(Diff.num,Nf); WRed = zeros(Diff.num);
    FermInd = zeros(Diff.num,1); % Storage for fermion indices that are being shifted.
    FermDest = zeros(Diff.num,1); % Storage for fermion destinations (in FermLoc vector).
    for d = 1:Diff.num
        if Diff.type == 0 % Type 0 - spin / fermion position exchange.
            % Diffs arising from a viable swap will list one pair of
            % positions, and Diff.val will be a [2x2] matrix, of which the
            % first row details the up fermion move.
            % First entry in Diff.pos will be initial site of up fermion.
            FermInd(d) = FermLocP(Diff.pos(d)+(mod(d-1,2)*N)); % Find index of fermion on site.
            FermDest(d) = Diff.pos(1+mod(d,Diff.num)) + (mod(d,2)*N); % Fermion destination
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            FermLocP(Diff.pos(d)+(mod(d-1,2)*N)) = 0; % Remove fermion from starting point in FermLocP.
            
        elseif Diff.type == 1 % Type 1 - up fermion move.
            FermInd(d) = FermLocP(Diff.pos(2*d-1)); % Independent fermion moves have numel(Diff.pos) = 2*Diff.num.
            FermDest(d) = Diff.pos(2*d); % Destination of up fermion.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            FermLocP(Diff.pos(2*d-1)) = 0; % Remove fermion from starting point in FermLocP.
            
        elseif Diff.type == -1 % Type -1 - down fermion move.
            FermInd(d) = FermLocP(Diff.pos(2*d-1)+N); % Independent fermion moves have numel(Diff.pos) = 2*Diff.num.
            FermDest(d) = Diff.pos(2*d) + N; % Destination of down fermion.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            FermLocP(Diff.pos(2*d-1)+N) = 0; % Remove fermion from starting point in FermLocP.
            
        elseif Diff.type == 2 % Type 2 - pair move.
            % Diffs arising from a viable pair move will list one pair of
            % positions, and Diff.val will be a [2x2] matrix, of which the
            % first row details the up fermion move.
            % First entry in Diff.pos will be initial site of pair.
            FermInd(d) = FermLocP(Diff.pos(1) + N*mod(d-1,2)); % First site is site being vacated.
            FermDest(d) = Diff.pos(2) + N*mod(d-1,2); % Destination of pair.
            FermLocP(FermDest(d)) = FermInd(d); % Add fermion to destination in FermLocP.
            FermLocP(Diff.pos(1) + N*mod(d-1,2)) = 0; % Remove fermion from starting point in FermLocP.
        end
    end
    % Construct reduced W matrix
    for k = 1:Diff.num
        for l = 1:Diff.num
            WRed(k,l) = SDetObj.DetMat(FermDest(k),FermInd(l));
        end
    end
    WRedInv = WRed^(-1);
    
    % Include phase factor in Diff struct.
    Ratio = Diff.sign * det(WRed); % Ratio of fermionic wavefunction portion
    if isnan(Ratio) || isinf(Ratio)
        Ratio = 1; % Likely stuck in an unviable configuration and should exit as soon as possible.
        % Setting to 1 ensures escape is possible without potentially
        % ruining some expectation values.
    end
    
    % Prepare updates for WFe in case move is accepted
    for k = 1:Diff.num
        for l = 1:Diff.num
            for f = 1:Nf
                BFe(k,f) = BFe(k,f) - (WRedInv(k,l) * (SDetObj.DetMat(FermDest(l),f)-(f==FermInd(l))));
            end
            % BFe(k,:) = BFe(k,:) - (WRedInv(k,l) * Ansatz.WFe(FermDest(l),:));
            % BFe(k,FermInd(l)) = BFe(k,FermInd(l)) + WRedInv(k,l);
        end
        DetMatUpdate = DetMatUpdate + (SDetObj.DetMat(:,FermInd(k)) * BFe(k,:));
    end
    % Sanity check for any unreasonably large or small elements
    cap = 1e30; min = 1e-30;
    DetMatUpdate(abs(DetMatUpdate)>cap) = 0;
    DetMatUpdate(abs(DetMatUpdate)<min) = 0;
else
    Ratio = 1;
end

% Store updates in struct for compatibility with current PsiRatio
% implementations
Update.FermLoc = FermLocP; Update.DetMatD = DetMatUpdate;