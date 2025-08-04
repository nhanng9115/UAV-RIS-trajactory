function [rate_min,gamma] = compute_rate(Ups,w,b,sys,chan)

% global M K T H sigma2
K = sys.K; T = sys.T; sigma2 = sys.sigma2; h = chan.h;

rate = zeros(K,1); gamma = zeros(K,T);
for k = 1:K
    for t = 1:T
        intf_power = sigma2(k,t);
        gamma(k,t) = abs(h(:,k,t)'*w(:,t))^2 / intf_power;
        rate(k) = rate(k) + 1/T*b(k,t)*log(1 + gamma(k,t));
    end
end

rate_min = min(rate);

end % EOF