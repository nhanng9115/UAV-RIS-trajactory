%% Result Proposed method
function [obj_value,Ups_new,btilde_new] = opt_actRIS(Ups_old,btilde_old,b_old,sys,cons,chan)

% global ops_soc K M N T amax Na pmax_r sigma2_r

Q = cons.Q; Qtilde = cons.Qtilde; q = cons.q; qtilde = cons.qtilde; Xi = cons.Xi; c4 = cons.c4; c5 = cons.c5;

K = sys.K; T = sys.T; Na = sys.Na;
amax = sys.amax; pmax_r = sys.pmax_r; sigma2_r = sys.sigma2_r;
ops_soc = sys.ops_soc;

ZERO = 0;

yalmip('clear') % clearing YALMIPs internal database

%% Variables
alpha = sdpvar(Na,T,'full','complex');
tau = sdpvar(1,1,'full','real');
b = sdpvar(K,T,'full','real');
btilde = sdpvar(K,T,'full','real');

%% objective function
obj = tau;

%% Constraints
F = [];
F = [F, tau >= 0];
F = [F, abs(alpha) <= amax];
F = [F, btilde >= ZERO, b >= ZERO];


%% (72)
if Na > 0
    for t = 1:T
        F = [F, alpha(1:Na,t)'*Xi(1:Na,1:Na,t)*alpha(1:Na,t) <= pmax_r];
    end
end

%% (71)
for k = 1:K
    diff_rate = 0;
    for t = 1:T
        alpha_old = diag(Ups_old(1:Na,1:Na,t));
        
        Qtmp = sqrtm(Q(1:Na,1:Na,k,t) + Qtilde(1:Na,1:Na,k,t));
        Fqua = -get_Fqua(Qtmp*alpha(1:Na,t),zeros(Na,1),Qtmp*alpha_old);
        %Fqua = -get_Fqua(Qtmp*alpha(1:Na,t),zeros(Na,1),Qtmp*alpha_old);
        qtmp = q(1:Na,k,t) + qtilde(1:Na,k,t);
        RHS_1 = Fqua + 2*real(alpha(:,t)'*qtmp);
        %F = [F, b(k,t) <= RHS_1];
        
        rate_1 = log(RHS_1 + c4(k,t) + c5(k,t));
        
        Q1 = sqrtm(Qtilde(1:Na,1:Na,k,t));
        a1 = Q1*alpha(:,t);
        F = [F, sum(a1.*conj(a1)) <= b(k,t)];
        F = [F, btilde(k,t) <= b(k,t) + 2*real(alpha(1:Na,t)'*qtilde(1:Na,k,t))];
        D = real(btilde_old(k,t) + c5(k,t));
        rate_2 = log(D) + real(btilde(k,t) - btilde_old(k,t))/D;
        
        
        diff_rate = diff_rate + b_old(k,t)*(rate_1 - rate_2);
    end
    F = [F, 1/T*diff_rate >= tau];
end

%% Get results
% disp('*********** SUB ACTIVE RIS **************************');
optimize(F,-obj,ops_soc);
obj_value = double(obj);
alpha_new = double(alpha);
% Ups_new = diag(alpha_new);
btilde_new = double(btilde);
Ups_new = zeros(Na,Na,T);
for t = 1:T
    Ups_new(:,:,t) = diag(alpha_new(:,t));
end
end % EOF
