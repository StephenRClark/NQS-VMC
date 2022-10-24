function [Cfg_cond] = TestCfgLength(Hilbert)
% This function tests self consistency between RandomCfg, Hilbert.N and
% FullCfg.
Cfg = Hilbert.RandomCfg();
Cfg_vec = Hilbert.FullCfg(Cfg);
Cfg_cond = (numel(Cfg_vec)==Hilbert.N);
if ~Cfg_cond
    disp(['Length of configuration (expect ' num2str(Hilbert.N) '): ' ...
        num2str(numel(Cfg_vec))]);
end
end