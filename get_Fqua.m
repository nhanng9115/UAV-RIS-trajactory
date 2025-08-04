function F_qua = get_Fqua(x,c,x0)
global M T K

if length(size(x)) > 2
    F_qua = sdpvar(M,T,K,'full','real');
    for m = 1:M
        for t = 1:T
            for k = 1:K
                F_qua(m,t,k) = 2*(c(:,k) - x0(:,m,t))'*(x(:,m,t) - x0(:,m,t)) - norm(x0(:,m,t) - c(:,k))^2;
            end
        end
    end
else
    F_qua = 2*(c - x0)'*(x - x0) - norm(x0 - c)^2;
end
end % EOF