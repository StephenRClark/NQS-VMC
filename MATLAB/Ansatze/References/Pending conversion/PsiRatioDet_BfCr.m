% --- General fermionic determinant wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioDet_BfCr(Ansatz,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% BfCr - backflow correlations, will alter determinant according to
% doublon-holon positions.
% ---------------------------------
% Format for fermionic determinant wavefunction addition:
% - Ansatz.Nf = number of fermions - set to Nv for spin models, specified at initialisation.
% - Ansatz.UFull = (2Nv x 2Nv) matrix - contains all available single particle orbitals.
% - Ansatz.UFe = (2Nv x Nf) matrix - obtained from diagonalisation of non-interacting terms of Hamiltonian.
% - Ansatz.WFe = (2Nv x Nf) matrix - elements are used for determinants in PsiRatio.
% - Ansatz.FermLoc = (2Nv x 1) vector - details locations of fermions by index for sign tracking purposes.
% Backflow correlation alterations:
% - Ansatz.UInv (Nf x 2Nv) replaces WFe as it is more convenient to work with.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = Ansatz.Nv; % Number of "visible" spins.
Nf = Ansatz.Nf; % Number of fermions.
UInv = Ansatz.UInv; BfP = Ansatz.BfP; UFe = Ansatz.UFe;

% Compute the fermionic part of the ratio.
FermLocP = Ansatz.FermLoc;

Cfg_vec = (FermLocP(1:Nv)>0) + (FermLocP(Nv+(1:Nv))>0); % Keeps track of density, though not spin.

if Diff.num ~= 0
    FermInd = zeros(Diff.num,1); FermStart = zeros(Diff.num,1); FermDest = zeros(Diff.num,1);
    QFe = eye(Nf); CFe = eye(Diff.num); UInvUpdateP = zeros(Nf,2*Nv); dURedP = zeros(Diff.num,Nf);
    for d = 1:Diff.num
        if Diff.val(d) == 1 % Moving fermion.
            if Diff.type == 1 % Up fermion move.
                FermInd(d) = FermLocP(Diff.start(d));
                FermStart(d) = Diff.start(d); FermDest(d) = Diff.dest(d);
                FermLocP(FermDest(d)) = FermInd; FermLocP(FermStart(d)) = 0;
            elseif Diff.type == -1 % Down fermion move.
                FermInd(d) = FermLocP(Diff.start(d)+Nv);
                FermStart(d) = Diff.start(d)+Nv; FermDest(d) = Diff.dest(d)+Nv;
                FermLocP(FermDest(d)) = FermInd; FermLocP(FermStart(d)) = 0;
            elseif Diff.type == 2 % Pair move - treat up fermion first.
                FermInd(d) = FermLocP(Diff.start(d)+(mod(d+1,2)*Nv));
                FermStart(d) = Diff.start(d)+(mod(d+1,2)*Nv); FermDest(d) = Diff.dest(d)+(mod(d+1,2)*Nv);
                FermLocP(FermDest(d)) = FermInd; FermLocP(FermStart(d)) = 0;
            end
            dURedP(d,:) = BfP(1) * (UFe(FermDest(d),:) - UFe(FermStart(d),:));
            if Cfg_vec(Diff.start(d)) == 2 % Contribution from doublons near start.
                for s = 1:numel(Diff.start_nn)
                    if Cfg_vec(Diff.start_nn(s)) == 0
                        dURedP(d,:) = dURedP(d,:) - BfP(2)*UFe(Diff.start_nn(s),:);
                    end
                end
            end
            if Cfg_vec(Diff.dest(d)) == 1 % Contribution from holons near destination.
                for t = 1:numel(Diff.dest_nn)
                    if Cfg_vec(Diff.dest_nn(t)) == 0
                        if Diff.type == 1 || (Diff.type == 2 && mod(d,2) == 1)
                            dURedP(d,:) = dURedP(d,:) + BfP(2)*UFe(Diff.dest_nn(t),:);
                        elseif Diff.type == -1 || (Diff.type == 2 && mod(d,2) == 0)
                            dURedP(d,:) = dURedP(d,:) + BfP(2)*UFe(Diff.dest_nn(t)+Nv,:);
                        end
                    end
                end
            end
        elseif Diff.val(d) == -2 % Fermion left behind on former doublon site loses backflow correlations.
            if Diff.type == 1 % Up fermion has left.
                FermInd = FermLocP(Diff.start(d)+Nv); FermStart(d) = Diff.start(d) + Nv;
                for s = 1:numel(Diff.start_nn)
                    if Cfg_vec(Diff.start_nn(s)) == 0
                        dURedP(d,:) = dURedP(d,:) - BfP(2)*UFe(Diff.start_nn(s)+Nv,:);
                    end
                end
            elseif Diff.type == -1 % Down fermion has left.
                FermInd(d) = FermLocP(Diff.start(d)); FermStart(d) = Diff.start(d);
                for s = 1:numel(Diff.start_nn)
                    if Cfg_vec(Diff.start_nn(s)) == 0
                        dURedP(d,:) = dURedP(d,:) - BfP(2)*UFe(Diff.start_nn(s)+Nv,:);
                    end
                end
            elseif Diff.val(d) == 2 % Doublons near vacated hole gain backflow correlation contributions.
                FermInd(d) = FermLocP(Diff.start(d) + (mod(d,2)*Nv)); % Both up and down fermions modified, order unimportant.
                FermStart(d) = Diff.start(d) + (mod(d,2)*Nv);
                HLoc = Diff.dest(d) + (mod(d,2)*Nv);
                dURedP(d,:) = - BfP(2) * UFe(HLoc,:);
            end
        end
    end
    
    for a = 1:Diff.num
        QFe = QFe + (UInv(:,FermStart(a)) * dURedP(a,:));
        for b = 1:Diff.num            
            CFe(a,b) = CFe(a,b) + dURedP(a,:)*UInv(:,FermStart(b));
            UInvUpdateP = UInvUpdateP - CFe(a,b)*(UInv(:,FermStart(a))*(dURedP(b,:).'*UInv));
        end
    end
    
    Ratio = Diff.sign * det(QFe);
    
    if isnan(Ratio) || isinf(Ratio)
        Ratio = 1; % Likely stuck in an unviable configuration and should exit as soon as possible.
        % Setting to 1 ensures escape is possible without potentially
        % ruining some expectation values.
    end  
else
    Ratio = 1;
end

% Store updates in struct for compatibility with current PsiRatio
% implementations
Update.FermLoc = FermLocP; Update.UInv = UInv;