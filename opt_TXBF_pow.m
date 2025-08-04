%% Result Proposed method
function [obj_value,w_new] = opt_TXBF_pow(Ups_old,b_old,sys,chan)

% global K M T pmax_r sigma2_r ops_soc sigma2 H h1 Na
K = sys.K; Nt = sys.Nt; Na = sys.Na; T = sys.T;
pmax = sys.pmax; pmax_r = sys.pmax_r;
sigma2_r = sys.sigma2_r; sigma2 = sys.sigma2; sigma2_u = sys.sigma2_u;
H1 = chan.H1; h = chan.h; h0 = chan.h0;
ops_soc = sys.ops_soc;

yalmip('clear')

p = sdpvar(T,1,'full','real');
tau = sdpvar(1,1,'full','real');

obj = tau;
F = [];
F = [F, tau >= 0, p >= 0, p <= pmax];


%% (59), another way, based on Q. Wu paper
% for k = 1:K
    sum_rate = 0;
    for t = 1:T
        bt = b_old(:,t);
        [~,k0] = max(bt);
        if sys.N == 0
            hk = h0(:,k0,t);
        else
            hk = h(:,k0,t);
        end
        rate_kt = log(1 + p(t)*norm(hk)^2/sigma2(k0,t));
        sum_rate = sum_rate + rate_kt;
    end
    F = [F, 1/T*sum_rate >= tau];
% end


%% (63)
if Na > 0
    for t = 1:T
        sum_pris = 0;
        for nn = 1:Na
            xi_n = sigma2_r + norm(H1(nn,:,t))^2 * p(t);
            sum_pris = sum_pris + abs(Ups_old(nn,nn,t))^2*xi_n;
        end
        F = [F, sum_pris <= pmax_r];
    end
end

%% Get results
disp('*********** TX BF **************************');

optimize(F,-obj,ops_soc)
obj_value = double(obj)
p_new = double(p);
w_new = zeros(Nt,T);
for t = 1:T
    bt = b_old(:,t);
    [~,k0] = max(bt);
    if sys.N == 0
        hk = h0(:,k0,t);
    else
        hk = h(:,k0,t);
    end
    w_new(:,t) = sqrt(pmax) * hk/norm(hk);
end



end % EOF




