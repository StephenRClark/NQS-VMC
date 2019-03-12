% --- Time-dependent variational Monte Carlo function ---

function [AnsatzObj,EvalIter] = TimeEvolveAnsatz(TE,SamplerObj,AnsatzObj)
% MMC - Multi Markov Chain - parallelises the sampling across the TE.Ncore
% cores assigned in TE.

Nc = TE.Ncore; Nsamp0 = round(SamplerObj.Nsamp/Nc,-2); % Round to nearest hundred samples.

% Add the TE/SR correlations to the evaluation list:
id = numel(SamplerObj.Operators) + 1; % Index for TE correlations in the list.
SamplerObj.Operators = {SamplerObj.Operators{1:(id-1)},...
    OperatorDg(AnsatzObj.Hilbert,AnsatzObj.Graph,@SR_EnLogDerivCorr),...
    OperatorDg(AnsatzObj.Hilbert,AnsatzObj.Graph,@SR_LogDerivCorr),...
    OperatorDg(AnsatzObj.Hilbert,AnsatzObj.Graph,@TE_EnLocSqCorr)}; % Concatenate them into the array.
% EvalIter logs all of the quantities specified in MCMC.eval for each time step
EvalAvgIter = cell(TE.Npass,id+1); EvalIter = cell(id+1,1); p = 1; SampFlag = 0;
EnRollAvg = 0; 

if Nc > 1
    MMCpool = parpool(Nc); % Assigning object handle for later deletion / closure.
end

while p<=TE.Npass
    
    NsampP = Nsamp0 + SampFlag * round(TE.ExtraSamp/Nc,-2);
    
    [SamplerObj] = SetNsamp(SamplerObj,NsampP);
    
    % If time dependent Hamiltonian is specified, calculate values
    % according to the assigned TFuncs.
    time = TE.TStep * TE.TScale * p; SamplerObj.Hamiltonian.HParams = HParams0;
    SamplerObj.Hamiltonian = TimeEvolveH(SamplerObj.Hamiltonian,time);
    
    if Nc > 1
        % Setup for distributed multiple Markov chain sampling:
        Sampler_s = cell(Nc,1); Ansatz_s = cell(Nc,1);
        % Prepare storage for each chain / core:
        EnAvg_s = zeros(Nc,1); dLogpAvg_s = cell(Nc,1);
        EvalAvg_s = cell(Nc,1);
        for c = 1:Nc
            Sampler_s{c} = SamplerObj;
            Ansatz_s{c} = AnsatzObj;
        end
        
        % Perform MCMC sampling of the current state:
        parfor (c = 1:Nc, Nc)
            [EnAvg_s(c),dLogpAvg_s{c},EvalAvg_s{c},~] = MCMCSample(Sampler_s{c},Ansatz_s{c});
        end
        % Recombine sampled quantities, averaging over the number of chains.
        % If a chain outputs a ludicrous energy, that contribution is removed.
        EnAvgMed = median(real(EnAvg_s)); EnAvgMean = mean(real(EnAvg_s)); EnAvgMin = min(real(EnAvg_s));
        fprintf('Mean energy: %d. Median energy: %d. Minimum energy: %d \n',EnAvgMean,EnAvgMed,min(real(EnAvg_s)));
        EnAvgDev = (abs(real(EnAvg_s) - EnAvgMin) / abs(EnAvgMin));
        % Large |EnAvgDev| indicates significant deviation from the mean/median:
        DelIndC = find(abs(real(EnAvgDev))>1/Nc);
        NcEff = Nc - numel(DelIndC);
        % Since deleting entries from cells just makes them empty, still loop
        % over all Nc entries, with the empty ones contributing nothing.
        EnAvg = 0; dLogpAvg = zeros(AnsatzObj.Np,1); EvalAvg = cell(id+2,1);
        if NcEff > 0
            for c = 1:Nc
                if sum(DelIndC==c) == 0
                    EnAvg = EnAvg + (1/NcEff) * EnAvg_s(c);
                    dLogpAvg = dLogpAvg + (1/NcEff) * dLogpAvg_s{c};
                    for i = 1:(id+2)
                        if isempty(EvalAvg{i})
                            EvalAvg{i} = (1/NcEff) * EvalAvg_s{c}{i};
                        else
                            EvalAvg{i} = EvalAvg{i} + (1/NcEff) * EvalAvg_s{c}{i};
                        end
                    end
                end
            end
        else
            [EnAvg,dLogpAvg,EvalAvg,~] = MCMCSample(SamplerObj,AnsatzObj);
        end
        
        % Calculate differences in energy versus last step:
        if p == 1
            EnRef = EnAvg; AnsRef = AnsatzObj;
        else
            dEV = abs(real(EnAvg) - real(EnRef)) / AnsatzObj.Hilbert.N;
            dER = dEV / abs(real(EnRef));
        end
        % Unpack the TE correlations from the evaluation list:
        EnLogDerivCorr = EvalAvg{id}; % <Eloc O_k*>.
        LogDerivCorr = EvalAvg{id+1}; % <O_k* O_k'>.
        % Energy variance calculation - useful as a diagnostic.
        EnLocSqCorr = EvalAvg{id+2};
        VarE = EnLocSqCorr - abs(EnAvg)^2;
        
        % Form the covariance matrix S:
        S = LogDerivCorr - conj(dLogpAvg) * dLogpAvg.'; % S_kk' = <O_k* O_k'> - <O_k*><O_k'>.
        
        % N.B: No diagonal regularisation applied for time evolution
        
        %   S_pc = S./sqrt(diag(S)*diag(S)');
        %   S_reg = S_pc + SR.epsilon*eye(Ansatz.Np);
        
        % Form the force vector F:
        F = EnLogDerivCorr - EnAvg*conj(dLogpAvg); % F_k = <E_loc O_k'> - <E_loc><O_k*>.
        
        % Possible issue - some S and F values are zero - have to filter out:
        ind = find(sum(abs(S),2) > TE.STol);
        F_N0 = F(ind);
        S_N0 = S(ind,ind);
        
        dP_N0 = - 1i * TE.TStep * TE.TScale * (pinv(S_N0) * F_N0);
        
        % Sanity check that all parameter changes are sensible:
        dP_N0(isinf(dP_N0)) = 0;
        dP_N0(isnan(dP_N0)) = 0;
        dP_N0(abs(dP_N0)<TE.dPTol) = 0;
        
        dP = zeros(AnsatzObj.NpTotal,1);
        dP(ind) = dP_N0;
        if p == 1
            AnsatzObj = AnsatzObj.PsiUpdate(AnsatzObj,dP); % Update the ansatz with the new parameters.
        else
            if abs(dEV) < TE.dEVTol && abs(dER) < TE.dERTol
                % Energy remains within tolerances - save ansatz and energy as
                % reference for rollback in case of error.
                AnsatzObj = AnsatzObj.PsiUpdate(AnsatzObj,dP);
                SampFlag = 0; AnsRef = AnsatzObj;
            else
                disp('Energy fluctuation outside of acceptable bounds detected - resampling.')
                SampFlag = 1; AnsatzObj = AnsRef;
            end
        end
        
    else
        disp('Multi-chain Markov sampling has given vastly disparate energies - resampling.');
        SampFlag = 1; AnsatzObj = AnsRef;
    end
    
    if SampFlag == 0
        % Store evaluated quantities inside EvalIter - energy goes in first
        % slot by default, energy fluctuations in second, others specified by
        % MCMC.eval follow.
        EvalAvgIter{p,1} = EnAvg/AnsatzObj.Hilbert.N; EvalAvgIter{p,2} = VarE / AnsatzObj.Hilbert.N;
        for k = 1:(id-1)
            EvalAvgIter{p,k+2} = EvalAvg{k};
        end
        fprintf('Energy calculated per site from sampling: %6f\n',EvalAvgIter{p,1});
        EnRef = EnAvg; % Rolling mean energy calculation.
        EnRollAvg = (EnRollAvg * (p-1) + EnAvg)/p; p = p + 1;
        fprintf('Current energy per site averaged over iterations: %6f\n',EnRollAvg / AnsatzObj.Hilbert.N);
    end
end

for k = 1:id+1
    EvalIter{k} = cell2mat(EvalAvgIter(:,k));
end

if Nc > 1
    delete(MMCpool); % To avoid complications if running TE in serial.
end