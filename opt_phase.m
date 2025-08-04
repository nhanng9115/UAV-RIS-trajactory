%% Result Proposed method
function [Phi_new] = opt_phase(b_old,w_old,Psi,sys,cons,chan)

% global K M N T Na h12tilde h0 h2 h1
T = sys.T; N = sys.N; Na = sys.Na;
h12tilde = cons.h12tilde; 
h0 = chan.h0; H1 = chan.H1; h2 = chan.h2;


Phi_new = zeros(N,N,T);

for t = 1:T
    bt = b_old(:,t);
    [~,k0] = max(bt);
    if Na == 0
        theta = angle(h0(:,k0,t)'*w_old(:,t)) - angle(h12tilde(:,k0,t));
    else
        theta = angle(h0(:,k0,t)'*w_old(:,t) + h2(:,k0)'*Psi(:,:,t)*H1(:,:,t)*w_old(:,t)) - angle(h12tilde(:,k0,t));
    end
    alpha = exp(1i.*theta);
    Phi_new(:,:,t) = diag(alpha);
end

%% Get results
% disp('*********** PASIVE RIS **************************');

end % EOF


