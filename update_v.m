function [v0hat0,v1hat0,v0check0,v1check0,v1tilde0] = update_v(v0,sys,cons)

% global M K T u r e0 e1 c0

K = sys.K; N = sys.N; T = sys.T;
u = sys.u; r = sys.r;
e0 = sys.e0; e1 = sys.e1;

v0hat0 = zeros(K,T); v1hat0 = zeros(T,1); v0check0 = zeros(K,T); v1check0 = zeros(T,1); v1tilde0 = zeros(T,1);

for t = 1:T
    for k = 1:K
        v0hat0(k,t) = norm(v0(:,t) - u(:,k))^(-e0/2);
        v0check0(k,t) = norm(v0(:,t) - u(:,k))^(-e0/2);
        %acheck0(k,t) = c0(k,t)*v0check0(k,t)^2;
    end
    v1hat0(t) = norm(v0(:,t) - r)^(-e1/2);
    v1check0(t) = norm(v0(:,t) - r)^(-e1/2);
    v1tilde0(t) = norm(v0(:,t) - r)^(-e1/2);
end


end %EOF