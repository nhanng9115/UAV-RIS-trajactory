function F_qua_2var = get_Fqua_2var_fast(x,x0)
global M T K
F_qua_2var = sdpvar(M,M,K,'full','real');
% F_qua_2var = zeros(M,M,K);
for t = 1:T
    for k = 1:K
        for m = 1:M
            for l = 1:M
                if l ~= m
                    F_qua_2var(m,l,k) = -norm(x0(:,m)-x0(:,l))^2 + 2*(x0(:,m)-x0(:,l)).'*(x(:,m)-x(:,l));
                else
                    F_qua_2var(m,l,k) = -Inf;
                end
            end
        end
    end
end
end % EOF