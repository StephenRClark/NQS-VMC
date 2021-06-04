% --- Stochastic Reconfiguration function ---

function [AnsatzObj,EvalIter] = SROptimise(SR,SamplerObj,AnsatzObj)
% MMC - Multi Markov Chain - parallelises the sampling across the SR.Ncore
% cores assigned in SR. Object oriented variant.

Nc = SR.Ncore; Nsamp0 = round(SamplerObj.Nsamp/Nc,-1); % Round to nearest ten samples.

% Add the SR correlations to the evaluation list:
id = numel(SamplerObj.Operators) + 1; % Index for SR correlations in the list.
SamplerObj.Operators = {SamplerObj.Operators{1:id-1},...
    OperatorDg(AnsatzObj.Hilbert,HypCub(1,0,0,0),@SR_EnLogDerivCorr),...
    OperatorDg(AnsatzObj.Hilbert,HypCub(1,0,0,0),@SR_LogDerivCorr)}; % Concatenate them into the array.

EvalIter = cell(SR.Npass,id);

% Logs for Ansatz and EnLoc for error correction purposes.
EneLog = zeros(SR.PSave,1);
AnsLog = cell(SR.PSave,1);
for p = 1:SR.PSave
    AnsLog{p} = AnsatzObj; % Preparations for error mitigation procedure.
end
% Flags and counters for error mitigation procedure.
SampFlag = 0; itcount = 0; encount = 0; tcount = 0; pshift = 0;
% Energy change sensitivity and transition rate tolerances used to detect
% if the ansatz has become stuck in a configuration. Encapsulated in
% StochasticReconfig object.
Mcheckmin = max(round(SR.PSave/5),5); PSave = max(SR.PSave,5);
if Nc > 1
    MMCpool = parpool(Nc); % Assigning object handle for later deletion / closure.
end

for p=1:SR.Npass
    
    % Increase sampling as time goes on to improve accuracy.
    if p <= 0.2*SR.Npass
        NsampT = Nsamp0;
    elseif p > 0.2*SR.Npass && p <= 0.8*SR.Npass
        NsampT = 1.25 * Nsamp0;
    else
        NsampT = 1.5 * Nsamp0;
    end
    
    NsampP = NsampT + SampFlag*round(SR.ExtraSamp/Nc,-2);
    
    [SamplerObj] = SetNsamp(SamplerObj,NsampP);
    
    if Nc > 1
        % Setup for distributed multiple Markov chain sampling:
        Sampler_s = cell(Nc,1); Ansatz_s = cell(Nc,1);
        % Prepare storage for each chain / core:
        EnAvg_s = zeros(Nc,1); dLogpAvg_s = cell(Nc,1);
        EvalAvg_s = cell(Nc,1); MRate_s = zeros(Nc,1);
        for c = 1:Nc
            Sampler_s{c} = SamplerObj;
            Ansatz_s{c} = AnsatzObj;
        end
        
        % Perform MCMC sampling of the current state:
        parfor (c = 1:Nc, Nc)
            [EnAvg_s(c),dLogpAvg_s{c},EvalAvg_s{c},MRate_s(c)] = MCMCSample(Sampler_s{c},Ansatz_s{c});
        end
        % Recombine sampled quantities, averaging over the number of chains.
        % If a chain outputs a ludicrous energy, that contribution is removed.
        EnAvgMed = median(real(EnAvg_s)); EnAvgMean = mean(real(EnAvg_s)); EnAvgMin = min(real(EnAvg_s));
        fprintf('Mean energy: %4f. Median energy: %4f. Minimum energy: %4f \n',EnAvgMean,EnAvgMed,min(real(EnAvg_s)));
        EnAvg = 0; dLogpAvg = zeros(AnsatzObj.NpTotal,1); EvalAvg = cell(id+1,1); MRate = 0;
        for c = 1:Nc
            EnAvg = EnAvg + (1/Nc) * EnAvg_s(c);
            dLogpAvg = dLogpAvg + (1/Nc) * dLogpAvg_s{c};
            for i = 1:(id+1)
                if isempty(EvalAvg{i})
                    EvalAvg{i} = (1/Nc) * EvalAvg_s{c}{i};
                else
                    EvalAvg{i} = EvalAvg{i} + (1/Nc) * EvalAvg_s{c}{i};
                end
            end
            MRate = MRate + (1/Nc) * MRate_s(c);
        end
    else
        [EnAvg,dLogpAvg,EvalAvg,MRate] = MCMCSample(SamplerObj,AnsatzObj);
    end
    
    % Unpack the SR correlations from the evaluation list:
    EnLogDerivCorr = EvalAvg{id}; % <Eloc O_k*>.
    LogDerivCorr = EvalAvg{id+1}; % <O_k* O_k'>.
    
    % Form the covariance matrix S:
    S = LogDerivCorr - conj(dLogpAvg) * dLogpAvg.'; % S_kk' = <O_k* O_k'> - <O_k*><O_k'>.
    SDiag = diag(S);
    
    % Regularise the S matrix (only use lambda_min for this version)
    if numel(SR.lambda_min) == 1 % Same minimum regularisation for all parameters.
        lambda_reg = max(SR.lambda_min,SR.lambda0*(SR.b^(p-pshift)));
        S_reg = S + lambda_reg * diag(diag(S));
    else % If different parameters require different minimisations, SR.lambda_min should be a (N_min x 2) matrix
        S_reg = S; s_start = 0;
        for s = 1:size(SR.lambda_min,1)
            lambda_reg = SR.lambda_min(s,1);
            SInds = s_start + (1:SR.lambda_min(s,2)); % SR.lmabda_min(s,2) is the number of elements lambda_min(s,1) applies to.
            S_reg(SInds,SInds) = S_reg(SInds,SInds) + lambda_reg * diag(diag(S(SInds,SInds))) + ...
                (SR.DFac * max(diag(S(SInds,SInds))) * eye(SR.lambda_min(s,2)));
            s_start = s_start + SR.lambda_min(s,2);
        end
    end
    
    % Form the force vector F:
    F = EnLogDerivCorr - EnAvg*conj(dLogpAvg); % F_k = <E_loc O_k'> - <E_loc><O_k*>.
    
    % Possible issue - some S and F values are zero - have to filter out:
    ind = find(sum(abs(S_reg),2) > SR.STol);
    if numel(ind) > 0
        F_N0 = F(ind)./sqrt(SDiag(ind));
        S_N0 = S_reg(ind,ind) ./ sqrt(SDiag(ind) * SDiag(ind)');
        
        dP_N0 = - SR.LRate * minres(AnsatzObj.Hilbert.N * S_N0, F_N0) ./ sqrt(SDiag(ind));
        
        % Sanity check that all parameter changes are sensible:
        dP_N0(isinf(dP_N0)) = 0;
        dP_N0(isnan(dP_N0)) = 0;
        dP_N0(abs(dP_N0)<SR.dPTol) = 0;
        
        dP = zeros(AnsatzObj.NpTotal,1);
        dP(ind) = dP_N0;
    elseif numel(ind) == 0
        disp('Logarithmic derivative magnitude has fallen below threshold for all parameters - no parameter changes applied.');
        dP = zeros(AnsatzObj.NpTotal,1);
        encount = encount + 1;
    end
    
    EvalIter{p,1} = EnAvg/AnsatzObj.Hilbert.N;
    for i = 1:(id-1)
        EvalIter{p,i} = EvalAvg{i};
    end
    
    if abs(EnAvg/AnsatzObj.Hilbert.N) > SR.ETol || isnan(EnAvg) || isinf(EnAvg) % In case a nonsensical energy is outputted.
        disp('Unphysical energy calculated - proposed parameter changes rejected.')
        dP = zeros(AnsatzObj.NpTotal,1);
        EvalIter{p,1} = SR.ETol;
    end
    
    AnsatzObj = AnsatzObj.PsiUpdate(dP); % Update the ansatz with the new parameters.
    
    % Error mitigation procedures detailed in following section.
    
    % Store a certain number of energies and parameter sets defined by SR.PSave.
    if p <= PSave
        AnsLog{p} = AnsatzObj;
        EneLog(p) = EvalIter{p,1};
    else % Examine if energy changes are unfavourable.
        dE = real(EvalIter{p,1}) - real(EneLog(end)); dEr = dE / abs(real(EneLog(end)));
        EneLog = [EneLog(2:end); EvalIter{p,1}]; AnsLog = {AnsLog{2:end}, AnsatzObj};
        if SampFlag == 1 % Extra samples are turned on.
            itcount = itcount + 1;
            if itcount == SR.PSave
                SampFlag = 0; itcount = 0; % Switch off extra samples after SR.PSave sweeps.
            end
        end
        if abs(dEr) < SR.ESens && p > SR.Npass/5
            encount = encount + 1; % Add to counter if change in energy ratio falls below SR.ESens.
        else
            encount = 0;
        end
        % SR.ESens should be set much lower than the average relative fluctuations in energy ~ 1/NSamp.
        if MRate < SR.MRTol
            tcount = tcount + 1;
        end
        % Move acceptance rate may fall to a very low value even if there
        % is no error, so make SR.MTol quite small - MTol * NSweep ~ 5-10
        
        if SampFlag == 0 % If errors occur in normal sampling, add extra samples and re-iterate, setting up AnsRef as a backup.
            if dEr > SR.dERTol && dE > SR.dEVTol
                fprintf('Unfavourable energy change proposed - backtracking by %d iterations and increasing sampling.\n',PSave);
                pshift = pshift + SR.PShift; AnsatzObj = AnsLog{1}; SampFlag = 1; itcount = 1; AnsRef = AnsatzObj;
            elseif tcount >= Mcheckmin
                fprintf('Move acceptance rate has fallen below %f for %d runs - backtracking by %d iterations and increasing sampling.\n',SR.MRTol,Mcheckmin,PSave);
                pshift = pshift + SR.PShift; AnsatzObj = AnsLog{1}; SampFlag = 1; itcount = 1; AnsRef = AnsatzObj; tcount = 0;
            elseif abs(dEr) < SR.ESens && encount >= PSave
                fprintf('Solution energy variation has fallen below specified threshold for %d consecutive runs.\n Assuming convergence has been reached and terminating further runs.\n',encount,PSave);
                EvalIter = EvalIter(1:p,:); % Truncate EnIter.
                break
            end
        elseif SampFlag == 1 % If errors occur in extra sampling period, just roll back changes to AnsRef and try again.
            if dEr > SR.dERTol && dE > SR.dEVTol && itcount > 1
                fprintf('Unfavourable energy proposed during extra sampling period - backtracking by %d iterations.\n',itcount);
                pshift = pshift + itcount; AnsatzObj = AnsRef; itcount = 1;
            elseif abs(dEr) < SR.ESens && encount >= PSave
                fprintf('Solution has remained stuck for %d runs during extra sampling period - backtracking by %d iterations.\n',encount,itcount);
                pshift = pshift + itcount; AnsatzObj = AnsRef; itcount = 1; encount = 0;
            elseif tcount >= Mcheckmin
                fprintf('Move acceptance rate has fallen below %f for %d runs during extra sampling period - backtracking by %d iterations.\n',SR.MRTol,Mcheckmin,itcount);
                pshift = pshift + itcount; AnsatzObj = AnsRef; itcount = 1; tcount = 0;
            end
        end
        
    end
    fprintf('Energy calculated per site from sampling: %6f\n',EvalIter{p,1});
end

if id == 1 % Only evaluating energy by iteration - convert to vector.
    EvalIter = cell2mat(EvalIter);
end

if Nc > 1
    delete(MMCpool); % To avoid complications if running SR in serial.
end