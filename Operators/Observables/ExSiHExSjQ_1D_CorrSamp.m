% --- Spin excitation wavefunction overlap evaluation function ---

function ExSiHExSjQ = ExSiHExSjQ_1D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of Hamiltonian matrix
% elements of excitations of particular momenta. The excitations are of the
% form of spin movements +/- S+/-(r+dR) S-/+(r) along with phase factors
% exp(iqr).

Nv = Ansatz.Nv; Q = (0:1:(Nv-1))*2/Nv; % Allowed momenta on lattice.

HParams = MCMC.HParams;

ExSiHExSjQ = zeros(Nv,Nv,Nv);

for q = 1:Nv
    Amp = zeros(1,Nv);
    HAmp = zeros(1,Nv);
    for dR = 1:Nv
        if dR == 1
            Cfg_vec = FullSpinCfg(Cfg);
            Amp(dR) = exp(1i*pi*Q(q)*(1:Nv)) * Cfg_vec / (2*sqrt(Nv));
        else
            [DiffSpSm,SpSmMatEls] = SpSm_1D_PBC_GT_CorrMatEls(HParams,Cfg,dR-1);
            % Apply S+(r-dR)S-(r) operators to the right.
            for d = 1:numel(SpSmMatEls)
                PsiRatio = MCMC.PsiRatio(Ansatz,DiffSpSm(d));
                Amp(dR) = Amp(dR) - PsiRatio * SpSmMatEls(d) * exp(1i*pi*DiffSpSm(d).pos(1)*Q(q)) / sqrt(Nv);
            end
            [DiffSmSp,SmSpMatEls] = SmSp_1D_PBC_GT_CorrMatEls(HParams,Cfg,dR-1);
            % Apply S-(r-dR)S+(r) operators to the right.
            for d = 1:numel(SmSpMatEls)
                PsiRatio = MCMC.PsiRatio(Ansatz,DiffSmSp(d));
                Amp(dR) = Amp(dR) + PsiRatio * SmSpMatEls(d) * exp(1i*pi*DiffSmSp(d).pos(1)*Q(q)) / sqrt(Nv);
            end
        end
    end
    
    [DiffH,HMatEls] = MCMC.HamMatEls(HParams,Cfg);
    % Apply Hamiltonian to find linked configurations, then consider
    % overlap of those configurations with excited state.
    for h = 1:numel(HMatEls)
        CfgP = DiffCfgPSpin(Cfg,DiffH(h));
        [PsiRatioH,Update] = MCMC.PsiRatio(Ansatz,DiffH(h));
        AnsatzP = MCMC.PsiCfgUpdate(Update,Ansatz);
        for dR = 1:Nv
            if dR == 1
                Cfg_vec = FullSpinCfg(CfgP);
                HAmp(dR) = HAmp(dR) + PsiRatioH * HMatEls(h) * exp(1i*pi*Q(q)*(1:Nv)) * Cfg_vec / 2;
            else
                [DiffSpSm,SpSmMatEls] = SpSm_1D_PBC_GT_CorrMatEls(HParams,CfgP,dR-1);
                % Apply S+(r-dR)S-(r) operators to the right.
                for d = 1:numel(SpSmMatEls)
                    PsiRatio = MCMC.PsiRatio(AnsatzP,DiffSpSm(d));
                    HAmp(dR) = HAmp(dR) - PsiRatioH * HMatEls(h) * PsiRatio...
                        * SpSmMatEls(d) * exp(1i*pi*DiffSpSm(d).pos(1)*Q(q));
                end
                [DiffSmSp,SmSpMatEls] = SmSp_1D_PBC_GT_CorrMatEls(HParams,CfgP,dR-1);
                % Apply S-(r-dR)S+(r) operators to the right.
                for d = 1:numel(SmSpMatEls)
                    PsiRatio = MCMC.PsiRatio(AnsatzP,DiffSmSp(d));
                    HAmp(dR) = HAmp(dR) + PsiRatioH * HMatEls(h) * PsiRatio...
                        * SmSpMatEls(d) * exp(1i*pi*DiffSmSp(d).pos(1)*Q(q));
                end
            end
        end
    end
    ExSiHExSjQ(:,:,q) = (Amp')*HAmp / Nv;
end
