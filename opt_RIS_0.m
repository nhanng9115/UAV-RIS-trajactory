%% Result Proposed method
function [obj_value,gammatilde_new,Ntilde_new,Ups_new] = opt_RIS_0(b_old,Ups_old,gammatilde_old,Ntilde_old)

global ops_soc K M T N amax Na pmax_r AA sigma2_r
global Q Qtilde q qtilde Xi c4 c5

% if M > 1
%     f = 1e5;
%     Qtilde = Qtilde*f; Q = Q*f;  q = q*f; qtilde = qtilde*f; c4 = c4*f; c5 = c5*f;
% end

ZERO = 0;

yalmip('clear') % clearing YALMIPs internal database

%% Variables
alpha = sdpvar(N,T,'full','complex');
tau = sdpvar(1,1,'full','real');
gammatilde = sdpvar(M,T,K,'full','real');
Ntilde = sdpvar(M,T,K,'full','real');

%% objective function
obj = tau;

%% Constraints
F = [];

F = [F, tau >= ZERO, gammatilde >= ZERO, Ntilde >= ZERO];

F = [F, abs(alpha) <= amax];

%% (72)
if Na > 0
    for t = 1:T
        F = [F, alpha(:,t)'*Xi(:,:,t)*alpha(:,t) <= pmax_r*sigma2_r];
    end
end

%% (70a)
for k = 1:K
    bk = b_old(:,k,:); bk = bk(:);
    gammatildek_mat = gammatilde(:,:,k); gammatildek_vec = gammatildek_mat(:);
    F = [F, 1/T*bk'*log(1 + gammatildek_vec) >= tau]; % (50a)
end

%% (71)
Fqol = get_Fqol(Ntilde,gammatilde,Ntilde_old,gammatilde_old);
for t = 1:T
    alpha_old = diag(Ups_old(:,:,t));
    for k = 1:K
        for m = 1:M
            
            % (72a)
            Dmkt = alpha(:,t)'*(Qtilde(:,:,m,k,t) + 100)*alpha(:,t) + 2*real(alpha(:,t)'*qtilde(:,m,k,t)) + c5(m,k,t);
            F = [F, Dmkt + Fqol(m,t,k) <= ZERO];
            
            % (72b)
            Qtmp = sqrtm(Q(:,:,m,k,t));
            Fqua = get_Fqua(Qtmp*alpha(:,t),zeros(N,1),Qtmp*alpha_old);
            %F = [F, Fqua <= ZERO];
            F = [F, Ntilde(m,t,k)^2 + Fqua - 2*real(alpha(:,t)'*q(:,m,k,t)) - c4(m,k,t) <= ZERO];
        end
    end
end

%% Get results
disp('*********** RIS **************************');
optimize(F,-obj,ops_soc)
obj_value = double(tau)
alpha_new = double(alpha);
gammatilde_new = double(gammatilde);
Ntilde_new = double(Ntilde);
Ups_new = zeros(N,N,T);
for t = 1:T
    Ups_new(:,:,t) = diag(alpha_new(:,t));
    for n = 1:N
       if ~ismember(n,AA)
           Ups_new(n,n,t) = Ups_new(n,n,t)/abs(Ups_new(n,n,t));
           alpha_new(n,t) = alpha_new(n,t)/abs(alpha_new(n,t));
       end
    end
end

end % EOF

% tic
%     Xi_cells = num2cell(Xi,[1,2]);
%     alpha_cells = num2cell(alpha,1);
%     alpha_block = blkdiag(alpha_cells{1:T});
%
%     XXi = blkdiag(Xi_cells{1:T});
%     alpha_all = alpha(:);
%     LHS = alpha_all'*XXi*alpha_block;
%     F = [F, LHS <= pmax_r];
%     toc


