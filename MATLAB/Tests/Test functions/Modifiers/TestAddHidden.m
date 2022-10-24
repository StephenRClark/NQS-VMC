function [Params_sum_cond,Params_len_cond] = TestAddHidden(Modifier,ModParams)
    % Test assumes new parameters are zeroes.
    % Extracting original parameters:
    Params = Modifier.ParamList();
    % Testing AddHidden works correctly:
    NewNQSA2 = Modifier.AddHidden(ModParams);
    Params_read = NewNQSA2.ParamList();
    Params_sum_cond = (abs(sum(Params_read)-sum(Params))<1e-12);
    Params_len_cond = (numel(Params_read)>numel(Params));
    if ~Params_sum_cond || ~Params_len_cond
        disp(['Length of new Params: ' num2str(length(Params_read))]);
        disp(['New number of parameters: ' num2str(NewNQSA2.Np)]);
        disp(['Alpha of new NQS: ' num2str(NewNQSA2.Alpha)]);
        disp(['Size of OptInds: ' num2str(size(NewNQSA2.OptInds,1)) ' '...
            num2str(size(NewNQSA2.OptInds,2))]);
    end
end