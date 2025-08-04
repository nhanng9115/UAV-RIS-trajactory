
% run addpathyalmip.m
close all
clear %all
warning off
clc

% load simulation parameters
UE_fix = 1;

%% Hybrid RIS configuration
sys = config(UE_fix,1);
sys.hmin = 50; % min height of UAV
sys.hv = sys.hmin;% initial height of UAV
linestyle = {':k','--g','cs','md','bh','-cs','--m'};
G = db2pow(0); % antenna gain, = 0 if not considered

T = 50;

n_chan = 5;
n_iter = 5; % number of iterations
V = zeros(3,T);

sys.T = T;
sys.ops_soc = sdpsettings('solver','mosek','verbose',0,'debug',1) ;

trajectory = 1; TX_BF = 1; association = 1;

% figure
for RIS = 2

    obj_UE = zeros(n_iter,n_chan);
    obj_TXBF = zeros(n_iter,n_chan);
    obj_RIS = zeros(n_iter,n_chan);
    obj_traj = zeros(n_iter,n_chan);
    minrate = zeros(n_iter,n_chan);

    % configurations of RIS
    [sys.N,sys.Na,sys.amax,sys.AA,sys.ONE] = config_RIS(RIS,sys.N0);
    N = sys.N;

    for cc = 1:n_chan
        cc

        % load small scale-fading channels
        data_location = strcat('./channels/',num2str(32)); %mkdir(data_location);
        file_name = strcat(data_location,'/small_scale',num2str(cc));
        data = load(file_name);
        g0all = data.g0; G1all = data.G1; g2all = data.g2;
        chan.g0 = G^2*g0all(:,:,1:T); chan.G1 = G*G1all(1:N,:,1:T); chan.g2 = G*g2all(1:N,:,1:T);

        %% Initialization
        [b_old,gamma,w_old,Ups_old,v_old,v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old,btilde_old,sys,chan,cons] = initialize(sys,chan);
        gamma_old = gamma;
        % if cc == 1
        plot_UE_location(v_old,sys,1,linestyle{RIS+1})
        % end

        %% Start inner iterations
        for ii = 1:n_iter
            ii
            %% UE scheduling -------------------------------------------------------------------
            if association == 1
                [obj_UE(ii,cc),b_new] = opt_association(b_old,gamma,sys);
                b_old = b_new;
            end

            %% Power opt -------------------------------------------------------------------
            if TX_BF == 1
                [obj_TXBF(ii,cc),w_new] = opt_power(w_old,Ups_old,b_old,sys,chan);
                w_old = w_new;
                if N > 0
                    [cons.Q,cons.Qtilde,cons.q,cons.qtilde,cons.c4,cons.c5,cons.Xi,cons.h12tilde,~] = update_Q(w_old,Ups_old,sys,chan);
                end
                [cons.c0,cons.c1,cons.c2,cons.c3] = update_constant(Ups_old,w_old,sys,chan);
            end

            % traj
            if trajectory == 1
                if RIS == 0
                    [obj_traj(ii,cc),v_new,v0hat_new,v0check_new] = ...
                        opt_traj_noRIS(v_old,v0hat_old,v0check_old,b_old,sys,cons);
                else
                    [obj_traj(ii,cc),v_new,v0hat_new,v1hat_new,v0check_new,v1check_new,v1tilde_new] = ...
                        opt_traj(v_old,v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old,Ups_old,b_old,sys,cons);
                end
                % update variables and constants and channels
                v_old = v_new;
                [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
                [cons.c0,cons.c1,cons.c2,cons.c3] = update_constant(Ups_old,w_old,sys,chan);
                [v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old] = update_v(v_old,sys,cons);
            end

            % RIS
            if RIS > 0
                [obj_RIS(ii,cc),Ups_new,btilde_new] = RIS_schemes(RIS,Ups_old,w_old,v_old,btilde_old,b_old,sys,chan,cons);
                % update variables and constants and channels
                Ups_old = Ups_new; btilde_old = btilde_new;
                [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
                [cons.Q,cons.Qtilde,cons.q,cons.qtilde,cons.c4,cons.c5,cons.Xi,cons.h12tilde,~] = update_Q(w_old,Ups_old,sys,chan);
            end


            %% Compute min rate
            [minrate(ii,cc),gamma] = compute_rate(Ups_old,w_old,b_old,sys,chan);
        end
        plot_UE_location(v_old,sys,1,linestyle{RIS+1})
    end
end