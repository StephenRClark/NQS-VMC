function [Diff_cond] = TestPropMove(Hilbert)
% This function tests self consistency of Diff2Cfg and PropMove.
% Generate starting configuration:
Cfg = Hilbert.RandomCfg();
[Diff, CfgP] = Hilbert.PropMove(Cfg);
% Use Diff2Cfg to test for consistency:
CfgD = Hilbert.Diff2Cfg(Diff,Cfg);
CfgP_vec = Hilbert.FullCfg(CfgP);
CfgD_vec = Hilbert.FullCfg(CfgD);
dCfg = CfgD_vec - CfgP_vec;
Diff_cond = (sum(abs(dCfg))<1e-10);
if ~Diff_cond
    disp(['Configuration from PropMove: ' num2str(CfgP_vec(:).')]);
    disp(['Configuration from Diff2Cfg: ' num2str(CfgD_vec(:).')]);
    disp(['Sum of absolute configuration differences: ' num2str(sum(abs(dCfg)))])
end
end