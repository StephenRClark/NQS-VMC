function [Psi_len_cond,Psi_norm_cond] = TestPsiGenerate(Modifier,TestCfgMat)
% Testing PsiGenerate on medium set:
Psi_gen_m = Modifier.PsiGenerate(TestCfgMat);
Psi_norm_cond = (1-(Psi_gen_m'*Psi_gen_m))<1e-10;
Psi_len_cond = (size(Psi_gen_m,1)==size(TestCfgMat,1));
if ~Psi_len_cond || ~Psi_norm_cond
    disp(['Length of Psi: ' num2str(size(Psi_gen_m,1))]);
    disp(['Norm of Psi: ' num2str(Psi_gen_m'*Psi_gen_m)]);
end
end