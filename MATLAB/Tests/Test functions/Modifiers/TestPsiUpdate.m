function [Psi_update_cond] = TestPsiUpdate(Modifier)
% Testing PsiUpdate:
Params = Modifier.ParamList();
Modifier = Modifier.PsiUpdate(-Params);
Params_read = Modifier.ParamList();
Psi_update_cond = (sum(abs(Params_read)))<1e-15;
if ~Psi_update_cond
    disp(['Resulting parameters: ' num2str(Params_read).']);
end
end