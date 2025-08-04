%% Result Proposed method
function [obj_value,w_new,gamma_new] = opt_TXBF(w_old,gamma_old,Ups_old,b_old,sys,chan)

% global K M T pmax_r sigma2_r ops_soc sigma2 H h1 Na
K = sys.K; Nt = sys.Nt; Na = sys.Na; T = sys.T;
pmax = sys.pmax; pmax_r = sys.pmax_r; 
sigma2_r = sys.sigma2_r; sigma2 = sys.sigma2; sigma2_u = sys.sigma2_u;
H1 = chan.H1; h = chan.h; h0 = chan.h0;
ops_soc = sys.ops_soc;

yalmip('clear')

w = sdpvar(Nt,T,'full','complex');
gamma = sdpvar(K,T,'full','real');
tau = sdpvar(1,1,'full','real');

obj = tau;
F = [];
F = [F, tau >= 0, gamma >= 0];

for t = 1:T
    %for k = 1:K
        %w_kt = w(Nt*(k-1)+1:Nt*k,t);
        F = [F, norm(w(:,t))^2 <= pmax];
    %end
end

%% (63)
% if Na > 0
%     for t = 1:T
%         sum_pris = 0;
%         for nn = 1:Na
%             xi_n = sigma2_r + norm(H1(nn,:,t))^2 * norm(w(:,t))^2;
%             sum_pris = sum_pris + abs(Ups_old(nn,nn,t))^2*xi_n;
%         end
%         F = [F, sum_pris <= pmax_r];
%     end
% end

%% (59), another way, based on Q. Wu paper
for k = 1:K
    
    sum_rate = 0;
    for t = 1:T
        
        %w_old_2D = w_old_3D(:,:,t);
        w_old_t = w_old(:,t); % (13)
        
        if sys.N == 0
            hk = h0(:,k,t);
        else
            hk = h(:,k,t);
        end
    
        %Hhat = zeros(K*Nt,K*Nt);
        %Hhat(Nt*(k-1)+1:Nt*k,Nt*(k-1)+1:Nt*k) = hk*hk'; % (14)
        Hhat = hk*hk';
        
        F_qol = w_old_t'*Hhat*w_old_t/gamma_old(k,t)^2 * gamma(k,t) - 2*real(w_old_t'*Hhat*w)/gamma_old(k,t); % (17)
        F = [F, F_qol <= 0];
        F = [F, sigma2(k,t) + F_qol <= 0]; % (18) ??????????????????????????????????
        
        %rate_kt = log( 1 + w(:,t)'*Hhat*w(:,t) / sigma2(k,t));
        
        rate_kt = log(1 + gamma(k,t));
        sum_rate = sum_rate + b_old(k,t)*rate_kt;
    end
    F = [F, 1/T*sum_rate >= tau];
end

%% Get results
disp('*********** TX BF **************************');

optimize(F,-obj,ops_soc)
obj_value = double(obj)
gamma_new = double(gamma_old);
w_new = double(w);
% w_new = zeros(Nt,K,T);
% for t = 1:T
%     w_new_mat = double(w(:,t));
%     w_new(:,:,t) = reshape(w_new_mat, Nt, []);
% end


end % EOF




