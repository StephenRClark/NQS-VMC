% --- Spin excitation wavefunction overlap evaluation function ---

function ExSiExSjQ = ExSiExSjQ_1D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of overlaps of excitations of
% particular momenta. The excitations are of the form of spin movements
% +/- S+/-(r+dR) S-/+(r) along with phase factors exp(iqr).

Nv = Ansatz.Nv; Q = (0:1:(Nv-1))*2/Nv; % Allowed momenta on lattice.

HParams = MCMC.HParams;

ExSiExSjQ = zeros(Nv,Nv,Nv);

for q = 1:Nv
    AmpV = zeros(1,Nv);
    for dR = 1:Nv
        if dR == 1
           Cfg_vec = FullSpinCfg(Cfg);
           AmpV(dR) = exp(1i*pi*Q(q)*(1:Nv)) * Cfg_vec / 2;
        else
        [DiffSpSm,SpSmMatEls] = SpSm_1D_PBC_GT_CorrMatEls(HParams,Cfg,dR-1);
        % Apply S+(r-dR)S-(r) operators to the right.
        for d = 1:numel(SpSmMatEls)
            PsiRatio = MCMC.PsiRatio(Ansatz,DiffSpSm(d));
            AmpV(dR) = AmpV(dR) - PsiRatio * SpSmMatEls(d) * exp(1i*pi*DiffSpSm(d).pos(1)*Q(q));
        end
        [DiffSmSp,SmSpMatEls] = SmSp_1D_PBC_GT_CorrMatEls(HParams,Cfg,dR-1);
        % Apply S-(r-dR)S+(r) operators to the right.
        for d = 1:numel(SmSpMatEls)
            PsiRatio = MCMC.PsiRatio(Ansatz,DiffSmSp(d));
            AmpV(dR) = AmpV(dR) + PsiRatio * SmSpMatEls(d) * exp(1i*pi*DiffSmSp(d).pos(1)*Q(q));
        end
        end
    end
    ExSiExSjQ(:,:,q) = (AmpV')*AmpV / sqrt(Nv); 
end
