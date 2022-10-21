function [Params_list_cond] = TestParams_Load_List(Modifier,Params)
% Testing ParamLoad and ParamList:
Modifier = Modifier.ParamLoad(Params(:));
Params_read = Modifier.ParamList();
dParams = Params - Params_read;
Params_list_cond = (sum(abs(dParams))<1e-10);
if ~Params_list_cond
    disp(['Original parameters: ' num2str(Params(:).')]);
    disp(['Listed parameters: ' num2str(Params_read(:).')]);
end
end