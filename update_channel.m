function [sigma2,h0,H1,h2,h,PL2] = update_channel(v,Ups,sys,chan)

% global K M T N Na sigma2_r sigma2_u
% global u r e0 e1 e2 xi0 ONE
% global sigma2 h0 h1 h2 g0 g1 g2 PL2 H% to compute

K = sys.K; N = sys.N; Na = sys.Na; Nt = sys.Nt; T = sys.T;
u = sys.u; r = sys.r;
g0 = chan.g0; G1 = chan.G1; g2 = chan.g2;
e0 = sys.e0; e1 = sys.e1; e2 = sys.e2; xi0 = sys.xi0; ONE = sys.ONE;
sigma2_r = sys.sigma2_r; sigma2_u = sys.sigma2_u;

%% Channels h0 h1 h2
h0 = zeros(Nt,K,T); H1 = zeros(N,Nt,T); h2 = zeros(N,K); PL2 = zeros(K,1);
for t = 1:T
    for k = 1:K
        d0 = norm(u(:,k) - v(:,t));
        PL0 = xi0*d0^(-e0);
        h0(:,k,t) = sqrt(PL0)*g0(:,k,t);
        
        if N > 0
            d1 = norm(v(:,t) - r);
            PL1 = xi0*d1^(-e1);
            H1(:,:,t) = sqrt(PL1).*G1(:,:,t);
            
            d2 = norm(u(:,k) - r);
            PL2(k) = xi0*d2^(-e2);
            h2(:,k) = sqrt(PL2(k)).*g2(:,k);
        end
    end
end

%% effective channel
h = zeros(Nt,K,T);
for t = 1:T
    for k = 1:K
        h(:,k,t) = h0(:,k,t) + (h2(:,k)'*Ups(:,:,t)*H1(:,:,t))';
    end
end

%% effective noise
sigma2 = zeros(K,T);
for t = 1:T
    for k = 1:K
        if N > 0 && Na > 0
            sigma2(k,t) = sigma2_u + sigma2_r*norm(h2(:,k)'*Ups(:,:,t)*ONE)^2;
        else
            sigma2(k,t) = sigma2_u;
        end
    end
end
end % EOF