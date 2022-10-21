function [dLogp_cond] = TestdLogp(Modifier,Cfg_test)
% Testing dLogp gives correct number of values:
dLogp = Modifier.LogDeriv(Cfg_test);
dLogp_cond = (length(dLogp) == Modifier.Np);
if ~dLogp_cond
    disp(['Length of dLogp: ' num2str(length(dLogp))]);
    disp(['Number of parameters: ' num2str(Modifier.Np)]);
end
end