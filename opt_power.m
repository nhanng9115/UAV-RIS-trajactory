%% Result Proposed method
function [obj_value,w_new] = opt_power(p_old,Ups_old,b_old,sys,chan)

K = sys.K; N = sys.N; Na = sys.Na; T = sys.T;
pmax = sys.pmax; pmax_r = sys.pmax_r; sigma2_r = sys.sigma2_r; sigma2 = sys.sigma2;
H1 = chan.H1; h = chan.h;
% ops_soc = sys.ops_soc;
% 
% % pmax
% yalmip('clear')
% p = sdpvar(1,T,'full','real');
% tau = sdpvar(1,1,'full','real');
% 
% obj = tau;
% F = [];
% 
% F = [F, tau >= 0];
% F = [F, p >= 0, p <= pmax];
% 
% %% (63)
% if Na > 0
%     for t = 1:T
%         sum_active_power = 0;
%         for nn = 1:Na
%             w_bar = abs(H1(nn,:,t).').^2;
%             xi_n = sigma2_r + p(t).'*w_bar;
%             sum_active_power = sum_active_power + abs(Ups_old(nn,nn,t))^2*xi_n;
%         end
%         F = [F, sum_active_power <= pmax_r];
%     end
% end
% 
% %% (59), another way, based on Q. Wu paper
% for k = 1:K
%     diff_rate = 0;
%     for t = 1:T
%         absH = abs(h(:,k,t)).^2;
%         % compute denumerator of D_kj
%         
%         % compute rate 2
%         
%         rate_2 = log(sigma2(k,t));
%         
%         rate_1 = log(sigma2(k,t) + sum(p(t).*absH));
%         diff_rate = diff_rate + b_old(k,t)*(rate_1 - rate_2);
%     end
%     F = [F, 1/T*diff_rate >= tau];
% end
% 
% %% Get results
% % disp('*********** POWER **************************');
% 
% optimize(F,-obj,ops_soc);
% obj_value = double(obj);

% p_new = double(p);
% gamma_new = double(gamma);

obj_value = 0;

Nt = sys.Nt;
if Nt  == 1
    w_new = p_new;
else
    w_new = zeros(Nt,T);
    for t = 1:T
        bt = b_old(:,t);
        [~,k0] = max(bt);
        hk = h(:,k0,t);
        
        % compute p
        % if Na > 0
        sum_pris = 0; sum_normh1 = 0;
        for nn = 1:Na
            sum_pris = sum_pris + abs(Ups_old(nn,nn,t))^2*sigma2_r;
            sum_normh1 = sum_normh1 + norm(H1(nn,:,t))^2;
        end
        ptilde = (pmax_r - sum_pris)/sum_normh1;
        p = min(pmax,ptilde);
        w_new(:,t) = sqrt(p) * hk/norm(hk);
    end
end

end % EOF




