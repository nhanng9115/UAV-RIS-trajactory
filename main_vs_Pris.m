
% run addpathyalmip.m
% close all
clear all
clc

% load simulation parameters
global g0 g1 g2 pmax_r linestyle sigma2_r ops_soc M K T N0

config;
n_chan = 5;
pmax = db2pow(20); M = 1; T = 50; N = N0; G = db2pow(3); 

pr_max_vec = -6:2:6; % number of iterations
ops_soc = sdpsettings('solver','mosek','verbose',0,'debug',1) ;
scheduling = 1; traj = 1; power = 0; RIS_opt = 1; run_RIS = 1;

for rr = 3
    config_RIS(rr)
    minrate = zeros(length(pr_max_vec),n_chan);
    
    for cc = 1:n_chan
        cc
        
        %% load small scale-fading channels
        data_location = strcat('./channels'); %mkdir(data_location);
        file_name = strcat(data_location,'/small_scale',num2str(cc));
        data = load(file_name);
        g0all = data.g0; g1all = data.g1; g2all = data.g2;
        g0 = G^2*g0all(1:M,:,1:T); g1 = G*g1all(1:N,1:M,1:T); g2 = G*g2all(1:N,:,1:T);
        
        for pp = 1:length(pr_max_vec)
            %pmax = db2pow(pmax_vec(pp))
            pmax_r = db2pow(pr_max_vec(pp))/sigma2_r;
            %% Start inner iterations
            [minrate(pp,cc)] = optimize_alg(rr,power,traj,scheduling,run_RIS,pmax);
            if rr < 2 % no RIS and passive RIS
               minrate(2:length(pr_max_vec),cc) = minrate(1,cc);
               break;
            end
            %alpha_old
        end

    end
    
    minrate_mean = mean(minrate,2);
    plot(pr_max_vec,minrate_mean, linestyle, 'LineWidth', 2); hold on
end
xlabel('Maximum RIS transmit power $\rho_{\mathrm{r,max}}$ [dBm]','Interpreter','latex')
ylabel('Minimum rate [bits/s/Hz]')
legend('No RIS','Passive RIS','Active RIS','Hybrid RIS')