function [Psi_rat_cond] = TestPsiRatio(Modifier,Cfg_test,Diff,RatioTester)
% Test that PrepPsi runs:
Modifier = Modifier.PrepPsi(Cfg_test);
% Testing PsiRatio runs:
[Ratio_df, ~] = Modifier.PsiRatio(Diff);
% Testing PsiGenerate on small set:
Psi_gen = Modifier.PsiGenerate(RatioTester);
Ratio_wf = Psi_gen(2)/Psi_gen(1);
dRatio = Ratio_df - Ratio_wf;
Psi_rat_cond = abs(dRatio/Ratio_wf) < 1e-10;
if ~Psi_rat_cond
    disp(['Ratio from PsiGenerate: ' num2str(Ratio_wf)]);
    disp(['Ratio from PsiRatio: ' num2str(Ratio_df)]);
end
end