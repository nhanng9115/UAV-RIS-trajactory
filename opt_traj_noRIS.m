%% Result Proposed method
function [obj_value,v_new,v0hat_new,v0check_new] = opt_traj_noRIS(v_old,v0hat_old,v0check_old,b_old,sys,cons)

ZERO = 0;

% global ops_soc K M T e0 u dmin dmax c0 Dmin hmin hmax Smax sigma2
K = sys.K; T = sys.T;
e0 = sys.e0; u = sys.u; dmin = sys.dmin; dmax = sys.dmax; hmin = sys.hmin; Smax = sys.Smax;
c0 = cons.c0;
sigma2 = sys.sigma2;
ops_soc = sys.ops_soc;

yalmip('clear') % clearing YALMIPs internal database

%% Variables
v = sdpvar(3,T,'full','real');
tau = sdpvar(1,1,'full','real');
v0hat = sdpvar(K,T,'full','real');
x0 = sdpvar(K,T,'full','real');


%% objective function
obj = tau;

%% Constraints
F = [];
F = [F, tau >= 0];
F = [F, v0hat >= ZERO, x0 >= ZERO, v0hat >= ZERO];


F = [F, v([1,2],:,:) <= dmax, v([1,2],:,:) >=  dmin];
% F = [F, hmin <= v(3,:,:) <= hmax]; % 3D trajectory
F = [F, v(3,:) == hmin]; % 2D

%% (19c)
F = [F, v(:,1) == v(:,T)];

%% (19d)
for t = 1:T-1
    F = [F, cone(v(:,t+1) - v(:,t),Smax)];
end

%% compute Fpow for constraints (43), (44), (46) and (47)
F00_pow = get_Fpow(v0hat,-2/e0,v0hat_old); % for (40)

% F = [F, F00_pow <= ZERO, F01_pow <= ZERO];


%% (33)
for k = 1:K
    diff_rate = 0;
    for t = 1:T
        
        %% (43a)
        F = [F, cone(v(:,t) - u(:,k), -F00_pow(k,t))]; % (40)
        
        
        %% (31a,b)
        Fquav0hat = v0hat_old(k,t).^2 + 2*v0hat_old(k,t).*(v0hat(k,t) - v0hat_old(k,t));
        rate_1 = log(sigma2(k,t) + sum(c0(k,t).*Fquav0hat));
        %% RATE
        rate_2 = log(sigma2(k,t));
        
        diff_rate = diff_rate + b_old(k,t)*(rate_1 - rate_2);
        
    end
    F = [F, 1/T*diff_rate >= tau];
end
%% ---------------------------------------------------------------------------------


% %% Get results
% disp('*********** TRAJECTORY **************************');
optimize(F,-obj,ops_soc);
obj_value = double(obj);
v_new = double(v);
v0hat_new = double(v0hat);

v0check_new = 0;
end % EOF




