function F_qua_2var = get_Fqua_2var(x,x0)
global M T K
F_qua_2var = sdpvar(M,M,T,K,'full','real');
% F_qua_2var = zeros(M,M,T,K);
for t = 1:T
    for k = 1:K
        for m = 1:M
            for l = 1:M
                if l ~= m
                    F_qua_2var(m,l,t,k) = norm(x0(:,m,t)-x0(:,l,t))^2 - 2*(x0(:,m,t)-x0(:,l,t)).'*(x(:,m,t)-x(:,l,t));
                else
                    F_qua_2var(m,l,t,k) = -Inf;
                end
            end
        end
    end
end
end % EOF