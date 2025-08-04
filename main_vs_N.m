
% run addpathyalmip.m
% close all
clear% all
clc

% load simulation parameters
global g0 g1 g2 linestyle N
scheduling = 1; traj = 1; power = 1; RIS_opt = 1;
config;
n_chan = 1;
N_vec = 20:20:60; % number of iterations
% figure
for rr = 4
    
    minrate = zeros(length(N_vec),n_chan);
    
    for cc = 1:n_chan
        cc
        
        %% load small scale-fading channels
        data_location = strcat('./channels'); %mkdir(data_location);
        file_name = strcat(data_location,'/small_scale',num2str(cc));
        load(file_name,'g0','g1','g2');
        g1_all = g1; g2_all = g2;
        
        for nn = 1:length(N_vec)
            N = N_vec(nn)
            g1 = g1_all(1:N,:); g2 = g2_all(1:N,:); 
            config_RIS(rr);
            minrate(nn,cc) = optimize_alg(rr,scheduling,RIS_opt,traj,power)
        end
        
    end
    
    minrate_mean = mean(minrate,2);
    plot(N_vec,minrate_mean, linestyle, 'LineWidth', 2); hold on
end
xlabel('Maximum UAV transmit power $\rho_{\mathrm{r,max}}$ [dBm]','Interpreter','latex')
ylabel('Minimum rate [bits/s/Hz]')
legend('No RIS','Passive RIS','Active RIS','Hybrid RIS')