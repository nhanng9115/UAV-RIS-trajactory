function [obj_value,b_new] = opt_association(b_old,gamma_old,sys)

% global K M T ops_soc

K = sys.K; T = sys.T;
ops_soc = sys.ops_soc; 

% mu = 50;

yalmip('clear')
b = sdpvar(K,T,'full','real');
tau = sdpvar(1,1,'full','real');

obj = tau; 

%% constraints
F = [];

F = [F, tau >= 0];
F = [F, b >= 0, b <= 1];
F = [F, sum(b,1) <= 1];

for k = 1:K
    rate_k = log(1+gamma_old(k,:));
    F = [F, 1/T*sum(b(k,:).*rate_k) >= tau]; % (57a)
end



%% Get results
% disp('*********** SCHEDULING **************************');

optimize(F,-obj,ops_soc);
obj_value = double(obj);
b_new = abs(double(b));


end %EOF