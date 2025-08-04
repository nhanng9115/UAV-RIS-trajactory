function [minrate1] = optimize_alg(RIS,opt_scheduling,opt_BF,opt_tra,sys,cc)

% global N Na
% sys.pmax = db2pow(pmax_dBm);

%% load small scale-fading channels
if RIS == 0
    file_name = strcat('./channels/',num2str(32),'/small_scale',num2str(cc));
else
    file_name = strcat('./channels/',num2str(sys.N),'/small_scale',num2str(cc));
end
data = load(file_name);
chan.g0 = data.g0; chan.G1 = data.G1; chan.g2 = data.g2;

%% Initialization
[b_old,gamma,w_old,Ups_old,v_old,v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old,btilde_old,sys,chan,cons] = initialize(sys,chan);

run_RIS = 1; run_TX = 1; run_loc = 1;
minrate = 0; ii = 0;
minrate0 = -Inf;

while (abs(minrate - minrate0) >= 1e-3 && ii <= 20) || (ii <= 5)
    
    minrate0 = minrate;  ii = ii + 1;
    
    %% UE scheduling -------------------------------------------------------------------
    if opt_scheduling == 1
        [obj_UE,b_new] = opt_association(b_old,gamma,sys);
        b_old = b_new;
    end
    
    %% Power opt -------------------------------------------------------------------
    if opt_BF == 1
        [obj_TXBF,w_new] = opt_power(w_old,Ups_old,b_old,sys,chan);
        %[obj_TXBF,w_new,gamma_new] = opt_TXBF(w_old,gamma_old,Ups_old,b_old,sys,chan);
        %[obj_TXBF,w_new] = opt_TXBF_pow(Ups_old,b_old,sys,chan);
        w_old = w_new;% gamma_old = gamma_new;
        if sys.N > 0
            [cons.Q,cons.Qtilde,cons.q,cons.qtilde,cons.c4,cons.c5,cons.Xi,cons.h12tilde,~] = update_Q(w_old,Ups_old,sys,chan);
        end
        [cons.c0,cons.c1,cons.c2,cons.c3] = update_constant(Ups_old,w_old,sys,chan);
    end
    
    % traj
    if opt_tra == 1
        if RIS == 0
            [obj_traj,v_new,v0hat_new,v0check_new] = ...
                opt_traj_noRIS(v_old,v0hat_old,v0check_old,b_old,sys,cons);
        else
            [obj_traj,v_new,v0hat_new,v1hat_new,v0check_new,v1check_new,v1tilde_new] = ...
                opt_traj(v_old,v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old,Ups_old,b_old,sys,cons);
        end
        % update variables and constants and channels
        v_old = v_new;
        if v_new(3) < 50
            disp('************INFEASIBILITY*********');
            break
        end
        [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
        [cons.c0,cons.c1,cons.c2,cons.c3] = update_constant(Ups_old,w_old,sys,chan);
        [v0hat_old,v1hat_old,v0check_old,v1check_old,v1tilde_old] = update_v(v_old,sys,cons);
    end
    
    % RIS
    if RIS > 0
        [obj_RIS,Ups_new,btilde_new] = RIS_schemes(RIS,Ups_old,w_old,v_old,btilde_old,b_old,sys,chan,cons);
        % update variables and constants and channels
        Ups_old = Ups_new; btilde_old = btilde_new;
        [sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v_old,Ups_old,sys,chan);
        [cons.Q,cons.Qtilde,cons.q,cons.qtilde,cons.c4,cons.c5,cons.Xi,cons.h12tilde,~] = update_Q(w_old,Ups_old,sys,chan);
        %abs(diag(Ups_new(:,:,1)))
    end
    
    
    %% Compute min rate
    [minrate,gamma] = compute_rate(Ups_old,w_old,b_old,sys,chan);
end


minrate1 = minrate;

end % EOF