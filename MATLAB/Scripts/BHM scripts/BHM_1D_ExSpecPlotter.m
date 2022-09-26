N = 10; U = [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 25/4 20/3 15/2 35/4 10 12 16 20 30 50 100];

Q = ((0:N) - round(N/2))/N; urange = [1 5 8 14 19 23]; 

AnsStr = 'BECR-Jast FT';

load(['BHM 1D N ' num2str(N) ' U Scan ' AnsStr ' excitation observable values.mat']);

Nw = 1000; 

LorentzProf = @(w,w0,dw) dw ./ (pi * ((w - w0).^2 + dw^2));
GaussProf = @(w,w0,dw) (1 / (dw * sqrt(pi))) * exp( -((w-w0).^2) / (2*dw^2));

StructFac = zeros(N,Nw,numel(urange)); EneEx = zeros(N,N,numel(urange));

for u = 1:numel(urange)
    EScale = N/U(urange(u)); w = EScale * linspace(1/Nw,1,Nw); dw = 0.1/U(urange(u));
    ExOL = ExOpEval{u}; HExOL = HExOpEval{u};
    for q = 1:N
        [EigVecQ,EigValQ] = eig(HExOL(:,:,q),ExOL(:,:,q));
        EneExQ = real(diag(EigValQ/N - EneGS(u)));
        EneEx(q,:,u) = reshape(EneExQ,1,N);
        for n = 1:N
            if real(EneExQ(n))>0
                EigVecQN = EigVecQ(:,n); EigVecQN = EigVecQN/sqrt(sum(abs(EigVecQN.^2)));
                OLFac = abs(EigVecQN'*ExOL(:,1,q))^2;
                StructFac(q,:,u) =  StructFac(q,:,u) ...
                + OLFac * LorentzProf(w,EneExQ(n),dw);
            end
        end
    end
    % PInds = [(round(N/2)+1):N 1:(round(N/2)+1)];
    % figure(u); imagesc(Q,w,StructFac(PInds,:,u).'); ax = gca; ax.YDir = 'normal'; colorbar;
    % xlabel('Q/2\pi'); ylabel('\omega');
    % title(['U = ' num2str(U(urange(u)))]);
    for q = 1:round((N+1)/2)
        figure(u); subplot(2,3,q);
        plot(w,StructFac(q,:,u));
        title(['U = ' num2str(U(urange(u))) ' Q/2\pi = ' num2str(Q(q+round(N/2)))]);
    end
end

% save(['BHM 1D N ' num2str(N) ' ' AnsStr ' structure factor values.mat'],'urange','StructFac','EneEx');
