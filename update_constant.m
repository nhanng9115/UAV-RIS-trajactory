function [c0,c1,c2,c3] = update_constant(Ups,w,sys,chan)

% global c0 c1 c2 c3% to compute
% global K M T N Na xi0 g0 G1 g2 h1 PL2  AA sigma2 sigma2_r

K = sys.K; N = sys.N; Na = sys.Na; T = sys.T;
g0 = chan.g0; G1 = chan.G1; xi0 = sys.xi0; h2 = chan.h2; H1 = chan.H1;
sigma2_u = sys.sigma2_r;

c0 = zeros(K,T); c1 = zeros(K,T); c2 = zeros(K,T);
for t = 1:T
    for k = 1:K
%         g0(:,k,t)'*w(:,t)
%         g0(:,k,t)'
%         w(:,t)
        c0(k,t) = abs(g0(:,k,t)'*w(:,t))^2*xi0;
        
        if N > 0
            c1(k,t) = abs(h2(:,k)'*Ups(:,:,t)*G1(:,:,t)*w(:,t))^2*xi0;
            c2(k,t) = 2*xi0*real(g0(:,k,t)'*w(:,t) * h2(:,k)'*Ups(:,:,t)*G1(:,:,t)*w(:,t));
        end
    end
end


c3 = zeros(Na,T);
if Na > 0
    for t = 1:T
        for nn = 1:Na
            c3(nn,t) = abs(H1(nn,:,t)*w(:,t))^2;
        end
    end
end

end % EOF