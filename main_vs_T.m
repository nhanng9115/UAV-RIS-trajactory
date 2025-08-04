
% run addpathyalmip.m
% close all
clear %all
clc

% load simulation parameters
sys = config(1,1);
sys.ops_soc = sdpsettings('solver','mosek','verbose',0,'debug',1) ;
% sys.T = 50;

n_chan = 8;
T_vec = [5,10,20:20:100]; % number of iterations
linestyle = {':k','--g','-cs','-rd','-bo',':rs','-r*'};
opt_BF = 1; opt_tra = 1; opt_scheduling = 1;

RIS_vec = [0:4];
% figure
for rr = 1:length(RIS_vec)
    RIS = RIS_vec(rr)
    
    % configurations of RIS
    [sys.N,sys.Na,sys.amax,sys.AA,sys.ONE] = config_RIS(RIS,sys.N0);
    
    minrate = zeros(length(T_vec),n_chan);
    
    
    for pp = 1:length(T_vec)
        %pp
        sys.T = T_vec(pp);
        parfor cc = 1:n_chan
            %cc
            [minrate(pp,cc)] = optimize_alg(RIS,opt_scheduling,opt_BF,opt_tra,sys,cc);
            %plot([1:length(minrate_vec)],minrate_vec, linestyle{RIS+1}, 'LineWidth', 1); hold on
        end
        
    end
    minrate_mean = mean(minrate,2)
    plot(T_vec,minrate_mean, linestyle{RIS+1}, 'LineWidth', 2); hold on
end

xlabel('Maximum transmit power of the UAV $(p_{\mathrm{max}}^{\mathrm{UAV}})$ [dBm]','Interpreter','latex','FontSize',12)
ylabel('Minimum rate [nats/s/Hz]','Interpreter','latex','FontSize',12)
legend('No RIS','Passive RIS','Hybrid RIS, 4 active elements')