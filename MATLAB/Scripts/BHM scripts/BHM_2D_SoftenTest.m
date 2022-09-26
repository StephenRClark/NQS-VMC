U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

AnsStr = 'BECR-Jast'; Nv = 9; Nh = 9; Nmax = 3; urange = [1 7 16]; 

load(['BHM 2D N ' num2str(Nv) ' U 1 exact ground state.mat'],'basis');

S_vals = 1:10; S_fid = zeros(numel(S_vals),numel(urange)); psi_s = cell(numel(S_vals),numel(urange)); W = size(basis,1);

for u = 1:numel(urange)
    
    load(['BHM 2D U ' num2str(U(urange(u))) ' N ' num2str(Nv) ' ' AnsStr ' Logs.mat'],'AnsatzObj');
    
    psi_mod = AnsatzObj.Modifier{1}.PsiGenerate(basis); JastFac = AnsatzObj.Modifier{1}.Js;
    
    % Set up correlator vectors and matrices using Jastrow factors.
    
    C_vis = cell(Nv,1); C_link = cell(Nv,Nh);
    
    for v = 1:Nv
        C_vis{v} = exp(-0.5*JastFac(v,v)*(0:Nmax).'.^2);
        for h = 1:Nh
            if h == v
                C_link{v,h} = eye(Nmax+1);
            else
                C_link{v,h} = exp(-0.5*JastFac(v,h)*((0:Nmax).'*(0:Nmax)));
            end
        end
    end
    
    % Introduce finite S in C_link{v,v}, which start as identity.
    
    for s = 1:numel(S_vals)
        for v = 1:Nv
            C_link{v,v} = exp(-S_vals(s)*(((0:Nmax).' - (0:Nmax)).^2));
        end
        psi = zeros(W,1);
        for b = 1:W
            cfg = basis(b,:);
            amp = 1;
            for v = 1:Nv
                amp = amp * C_vis{v}(cfg(v)+1);
            end
            h_amp = zeros(Nh,1);
            for h = 1:Nh
                hv = ones(Nmax+1,1);
                for v = 1:Nv
                    hv = hv .* (C_link{v,h}(cfg(v)+1,:).');
                end
                h_amp(h) = sum(hv);
            end
            amp = amp * prod(h_amp);
            psi(b) = amp;
        end
        psi_s{s,u} = psi / sqrt(sum(abs(psi).^2));
        S_fid(s,u) = abs(psi_s{s,u}'*psi_mod)^2;
    end
    
    figure(u); semilogy(S_vals,1-S_fid(:,u),'-o');
    xlabel('S'); ylabel('Infidelity'); 
    title(['U = ' num2str(U(urange(u)))]);    
end