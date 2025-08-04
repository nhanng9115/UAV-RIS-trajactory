function F_qua = get_Fqua_all(xall,call,x0all)
global M T K
F_qua = sdpvar(M,K,'full','real');

for m = 1:M
    for k = 1:K
        x = xall(:,m); x0 = x0all(:,m); c = call(:,k);
        F_qua(m,k) = 2*(c - x0)'*(x - x0) - norm(x0 - c)^2;
    end
end
end % EOF