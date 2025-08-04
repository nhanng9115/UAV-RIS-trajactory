function [D,Dtilde,d,dtilde,z,ztilde,Xi,h12tilde,btilde0] = update_Q(w0,Ups0,sys,chan)

% global Q Qtilde q qtilde c4 c5 Xi h12tilde% to compute
% global K M T N Na h0 H1 h2 sigma2_u sigma2_r AA ONE
test_derivation = 0;

K = sys.K; N = sys.N; Na = sys.Na; ONE = sys.ONE; T = sys.T;
sigma2_r = sys.sigma2_r; sigma2_u = sys.sigma2_r;
h0 = chan.h0; H1 = chan.H1; h2 = chan.h2;

%% compute H2tilde and h12tilde
H2tilde = zeros(N,N,K,T); h12tilde = zeros(N,K,T);
D = zeros(N,N,K,T); d = zeros(N,K,T);
Dtilde = zeros(N,N,K,T); dtilde = zeros(N,K,T);
z = zeros(K,T); ztilde = zeros(K,T);
btilde0 = zeros(K,T);

%% First, compute Q,q,c4 for all  k
for t = 1:T
    phi0 = zeros(N,1);
    phi0(Na+1:N) = diag(Ups0(Na+1:N,Na+1:N,t));
    for k = 1:K
        H2tilde(:,:,k) = diag(h2(:,k)');
        
        %% compute Q and q
        h0bar(k,t) = h0(:,k,t)'*w0(:,t);
        h1bar(:,k,t) = H1(:,:,t) * w0(:,t);
        
        h12tilde(:,k,t) = H2tilde(:,:,k) * h1bar(:,k,t);
        h0barbar(k,t) = h0bar(k,t) + phi0.'*h12tilde(:,k,t);
        
        D(:,:,k,t) = conj(h12tilde(:,k,t))*h12tilde(:,k,t).';
        d(:,k,t) = conj(h12tilde(:,k,t))*h0barbar(k,t);
        z(k,t) = abs(h0barbar(k,t))^2;
    end
end

%% Second, compute Qtilde, qtile, c5
for t = 1:T
    alpha0 = diag(Ups0(:,:,t));
    for k = 1:K
        Dtilde(:,:,k,t) = sigma2_r * ONE*conj(H2tilde(:,:,k))*H2tilde(:,:,k).'*ONE;
        dtilde(:,k,t) = 0;
        ztilde(k,t) = sigma2_u;
        
        % compute btilde
        btilde0(k,t) = alpha0'*Dtilde(:,:,k,t)*alpha0 + 2*real(alpha0'*dtilde(:,k,t));
    end
end

%% compute Xi
Xi = zeros(N,N,T);
for t = 1:T
    if Na > 0
        for nn = 1:Na
            xi_n = sigma2_r + norm(H1(nn,:,t))^2 * norm(w0(:,t))^2;
            Xi(nn,nn,t) = xi_n;
        end
    end
end

%% check the derivations in (68) and (69)
if test_derivation == 1
    
    for t = 1:T
        Ups0_t = Ups0(:,:,t);
        Psi0_t = Ups0_t.*ONE;
        psi_t = diag(Psi0_t);
        
        for k = 1:K
            
            h0bar(k,t) = h0(:,k,t)'*w0(:,t);
            h1bar(:,k,t) = H1(:,:,t) * w0(:,t);
            
            h12tilde(:,k,t) = H2tilde(:,:,k) * h1bar(:,k,t);
            h0barbar(k,t) = h0bar(k,t) + phi0.'*h12tilde(:,k,t);
           
            
            % test numerator of SINR ==> ok
            Num = abs( h0bar(k,t) + h2(:,k)'*Ups0_t*h1bar(:,k,t) )^2;
            Num_ana = real(psi_t'*D(:,:,k,t)*psi_t + 2*real(psi_t'*d(:,k,t)) + z(k,t));
            
            Denom = sigma2_r*norm(h2(:,k)'*Psi0_t)^2 + sigma2_u
            Denom_ana = real(ztilde(k,t) + psi_t'*Dtilde(:,:,k,t)*psi_t + 2*real(psi_t'*dtilde(:,k,t)))

        end
    end
end

end % EOF