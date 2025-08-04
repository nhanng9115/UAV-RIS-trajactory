function F_bil = get_Fbil(xall,yall,call,x0all,y0all)

global M K T
if length(size(xall)) == 2
    if call > 0
        F_bil = 0.25*abs(call)*(xall-yall)^2 + 0.25*abs(call)*(x0all+y0all)^2 - 0.5*(x0all+y0all)*(xall+yall);
    else
        F_bil = 0.25*abs(call)*(xall-yall)^2 + 0.25*abs(call)*(x0all+y0all)^2 - 0.5*(x0all+y0all)*(xall+yall);
    end
else
    F_bil = sdpvar(M,T,K,'full','real');
    for m = 1:M
        for t = 1:T
            for k = 1:K
                x = xall(m,t,k); x0 = x0all(m,t,k); y = yall(m,t,k); y0 = y0all(m,t,k); c = call(m,t,k);
                if c > 0
                    t0 = y0/x0;
                    F_bil(m,t,k) = 0.5*abs(c)*(t0*x^2 + 1/t0*y^2);
                else
                    F_bil(m,t,k) = 0.25*abs(c)*(x-y)^2 + 0.25*abs(c)*(x0+y0)^2 - 0.5*(x0+y0)*(x+y);
                end
            end
        end
    end
end
end % EOF