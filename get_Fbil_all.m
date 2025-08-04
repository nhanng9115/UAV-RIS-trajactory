function F_bil = get_Fbil_all(xall,yall,call,x0all,y0all)

% t0 = y0./x0;
% F_bil_pos = 0.5*abs(c).*(t0.*x(:,:,:).^2 + 1./t0.*y.^2);
% F_bil_neg = 0.25*abs(c).*(x-y).^2 + 0.25*abs(c).*(x0+y0).^2 - 0.5*(x0+y0).*(x+y);
%
% c_pos = max(c,0); sign_pos = min(c_pos,1);
% c_neg = min(c,0); sign_neg = min(c_neg,0);
%
% F_bil = F_bil_pos.*sign_pos + F_bil_neg.*sign_neg;
global M K T

F_bil = sdpvar(M,K,'full','real');
for m = 1:M
    for k = 1:K
        x = xall(m,k); x0 = x0all(m,k); y = yall(m,k); y0 = y0all(m,k); c = call(m,k);
        if c > 0
            t0 = y0/x0;
            F_bil(m,k) = 0.5*abs(c)*(t0*x^2 + 1/t0*y^2);
        else
            F_bil(m,k) = 0.25*abs(c)*(x-y)^2 + 0.25*abs(c)*(x0+y0)^2 - 0.5*(x0+y0)*(x+y);
        end
    end
end

end % EOF