%% Result Proposed method
function [obj_value,v_new,v0hat_new,v1hat_new,v0check_new,v1check_new,v1tilde_new] = ...
    opt_traj(v_old,v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old,Ups_old,b_old,sys,cons)

ZERO = 0;

% global ops_soc K M T Na pmax_r e0 e1 u r AA dmin dmax c0 c1 c2 c3 Dmin hmin hmax sigma2 Smax sigma2_r
K = sys.K; T = sys.T; Na = sys.Na;
e0 = sys.e0; e1 = sys.e1; e2 = sys.e2;
u = sys.u; r = sys.r; dmin = sys.dmin; dmax = sys.dmax; hmin = sys.hmin; Smax = sys.Smax;
c0 = cons.c0; c1 = cons.c1; c2 = cons.c2; c3 = cons.c3;
sigma2 = sys.sigma2; sigma2_r = sys.sigma2_r; pmax_r = sys.pmax_r;
ops_soc = sys.ops_soc;

yalmip('clear') % clearing YALMIPs internal database

%% Variables
v = sdpvar(3,T,'full','real');
tau = sdpvar(1,1,'full','real');
v0hat = sdpvar(K,T,'full','real');
v1hat = sdpvar(T,1,'full','real');
x1check = sdpvar(T,1,'full','real');
ahat = sdpvar(K,T,'full','real');
v1check = sdpvar(T,1,'full','real');
v0check = sdpvar(K,T,'full','real');
x0 = sdpvar(K,T,'full','real');
v1tilde = sdpvar(T,1,'full','real');
x1tilde = sdpvar(T,1,'full','real');

%% objective function
obj = tau;

%% Constraints
F = [];
F = [F, tau >= 0];
F = [F, v0hat >= ZERO, v1hat >= ZERO];
F = [F, v0check >= ZERO, v1check >= ZERO, v1tilde >= ZERO];
F = [F, x0 >= ZERO, x1check >= ZERO, x1tilde >= ZERO, ahat >= ZERO];

F = [F, v([1,2],:,:) <= dmax, v([1,2],:,:) >=  dmin];
F = [F, v(3,:) == hmin]; % 2D
F = [F, v(:,1) == v(:,T)]; % (24b)

%% (19d)
for t = 1:T-1
    F = [F, norm(v(:,t+1) - v(:,t)) <= Smax];
end

%% compute Fpow for constraints (43), (44), (46) and (47)
F00_pow = get_Fpow(v0hat,-2/e0,v0hat_old); % for (43)
F10_pow = get_Fpow(v1hat,-2/e1,v1hat_old); % for (44)

%% (43b)
for t = 1:T
    F = [F, cone(v(:,t) - r, -F10_pow(t))];  % (41)
end

F01_pow = get_Fpow(v0check,4/e0,v0check_old); % for (46)
F11_pow = get_Fpow(v1check,4/e1,v1check_old); % for (47)
F = [F, F00_pow <= 0, F01_pow <= 0, F10_pow <= 0, F11_pow <= 0];



%% constraint (33)-(44)
for k = 1:K
    diff_rate = 0;
    for t = 1:T
        
        %% Constraint %% (43a)
        F = [F, cone(v(:,t) - u(:,k), -F00_pow(k,t))]; % (40)
        
        %% Constraint (43c), (43d)
        if c2(k,t) < 0
            F0_qua = get_Fqua(v(:,t),u(:,k),v_old(:,t));
            F = [F, -F0_qua >= x0(k,t)];
            F = [F, cone([1,0.5*(-F01_pow(k,t) - x0(k,t))], 0.5*(-F01_pow(k,t) + x0(k,t)))];  % (43)
            
            F1_qua = get_Fqua(v(:,t),r,v_old(:,t));
            F = [F, -F1_qua >= x1check(t)];
            F = [F, cone([1,0.5*(-F11_pow(t) - x1check(t))], 0.5*(-F11_pow(t) + x1check(t)))];  % (43)
        else
            F = [F, v0check(k,t) == v0hat(k,t), v1check(t) == v1hat(t)]; % (40)
        end
        
        %% (31a)
        Fquav0hat = v0hat_old(k,t)^2 + 2*v0hat_old(k,t)*(v0hat(k,t) - v0hat_old(k,t));
        Fquav1hat = v1hat_old(t)^2 + 2*v1hat_old(t)*(v1hat(t) - v1hat_old(t));
        %Fbilv0v1hat = get_FbilLB(v0check(k,t),v1check(t),c2(k,t),v0check_old(k,t),v1check_old(t));
        %RHS = c0(k,t)*Fquav0hat + c1(k,t)*Fquav1hat + abs(c2(k,t))*Fbilv0v1hat;
        %F = [F, ahat(k,t) <= RHS]; % (31a)
        
        if c2(k,t) > 0
            Fbilv0v1hat = get_FbilLB(v0check(k,t),v1check(t),c2(k,t),v0check_old(k,t),v1check_old(t));
            RHS = c0(k,t)*Fquav0hat + c1(k,t)*Fquav1hat + abs(c2(k,t))*Fbilv0v1hat;
            F = [F, ahat(k,t) <= RHS]; % (31a)
        else
            RHS = c0(k,t)*Fquav0hat + c1(k,t)*Fquav1hat;
            F_bil = get_Fbil(v0check(k,t),v1check(t),1,v0check_old(k,t),v1check_old(t));
            F = [F, ahat(k,t) + abs(c2(k,t))*F_bil <= RHS]; % (31a)
        end
        
        
        %% Constraint (33)
        % compute 1st rate term
        rate_1 = log(sigma2(k,t) + ahat(k,t));
        
        % compute second rate term
        rate_2 = log(sigma2(k,t));
        
        diff_rate = diff_rate + b_old(k,t)*(rate_1 - rate_2);
    end
    F = [F, 1/T*diff_rate >= tau];
end
%% ---------------------------------------------------------------------------------

%% (56)
F22_pow = get_Fpow(v1tilde,4/e1,v1tilde_old); % for (47)
if Na > 0
    for t = 1:T
        active_power = 0;
        for nn = 1:Na
            F = [F, cone([1,0.5*(-F22_pow(t) - x1tilde(t))], 0.5*(-F22_pow(t) + x1tilde(t)))];   % (46)
            xi_n = sigma2_r + c3(nn,t)*(v1tilde(t))^2;
            %xi_n = xi_n + c3(n,t)*x1(t)^2;
            active_power = active_power + abs(Ups_old(nn,nn,t))^2*xi_n;
        end
        F = [F, active_power <= pmax_r];
    end
end


%% Get results
% disp('*********** TRAJECTORY **************************');
optimize(F,-obj,ops_soc);
obj_value = double(obj);
v_new = double(v);
v0hat_new = double(v0hat);
v1hat_new = double(v1hat);
v0check_new = double(v0check);
v1check_new = double(v1check);
v1tilde_new = double(v1tilde);

end % EOF




