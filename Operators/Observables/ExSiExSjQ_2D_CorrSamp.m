% --- Spin excitation wavefunction overlap evaluation function ---

function ExSiExSjQ = ExSiExSjQ_2D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of overlaps of excitations of
% particular momenta. The excitations are of the form of spin movements
% +/- S+/-(r+dR) S-/+(r) along with phase factors exp(iqr).

Nv = Ansatz.Nv; L1 = Ansatz.Dim(1); L2 = Ansatz.Dim(2);
Qx = (0:1:(L1-1))*2/L1;  Qy = (0:1:(L2-1))*2/L2; % Allowed momenta on lattice.

HParams = MCMC.HParams;

ExSiExSjQ = zeros(Nv,Nv,Nv);

for q1 = 1:L1
    for q2 = 1:L2
        AmpV = zeros(1,Nv);
        for dR1 = 1:L1
            for dR2 = 1:L2
                if dR1 == 1 && dR2 == 1
                    Cfg_mat = reshape(FullSpinCfg(Cfg),L1,L2);
                    Cfg_mat = - exp(1i*pi * ( (Qx(q1) * ((1:L1).'*ones(1,L2))) + ...
                        (Qy(q2) * (ones(L1,1)*(1:L2))) )) .* Cfg_mat;
                    AmpV(dR1+(dR2-1)*L1) = sum(Cfg_mat(:)) / 2;
                else
                    [DiffSpSm,SpSmMatEls] = SpSm_2D_PBC_GT_CorrMatEls(HParams,Cfg,[dR1-1, dR2-1]);
                    % Apply S+(r-dR)S-(r) operators to the right.
                    for d = 1:numel(SpSmMatEls)
                        PsiRatio = MCMC.PsiRatio(Ansatz,DiffSpSm(d));
                        Rx = 1+mod(DiffSpSm(d).pos(1)-1,L1); Ry = ceil(DiffSpSm(d).pos(1)/L1);
                        AmpV(dR1+(dR2-1)*L1) = AmpV(dR1+(dR2-1)*L1) + PsiRatio * SpSmMatEls(d)...
                            * exp(1i*pi*(Qx(q1)*Rx + Qy(q2)*Ry));
                    end
                    [DiffSmSp,SmSpMatEls] = SmSp_2D_PBC_GT_CorrMatEls(HParams,Cfg,[dR1-1, dR2-1]);
                    % Apply S-(r-dR)S+(r) operators to the right.
                    for d = 1:numel(SmSpMatEls)
                        PsiRatio = MCMC.PsiRatio(Ansatz,DiffSmSp(d));
                        Rx = 1+mod(DiffSmSp(d).pos(1)-1,L1); Ry = ceil(DiffSmSp(d).pos(1)/L1);
                        AmpV(dR1+(dR2-1)*L1) = AmpV(dR1+(dR2-1)*L1) - PsiRatio * SmSpMatEls(d)...
                            * exp(1i*pi*(Qx(q1)*Rx + Qy(q2)*Ry));
                    end
                end
            end
        end
        
        q = q1 + (q2-1)*L1;
        ExSiExSjQ(:,:,q) = (AmpV')*AmpV / sqrt(Nv);
        
    end
end