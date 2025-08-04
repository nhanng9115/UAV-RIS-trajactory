function [b0,gamma0,w0,Ups0,v0,v0hat0,v1hat0,v0check0,v1check0,v1tilde0,btilde0,sys,chan,cons] = initialize(sys,chan)

% global K M N T Na AA dmin dmax u hv UAV_fix sigma2 H r vmax

K = sys.K; N = sys.N; T = sys.T; Na = sys.Na; Nt = sys.Nt; pmax = sys.pmax;
vmax = sys.vmax; dmax = sys.dmax; u = sys.u; hv = sys.hv; r = sys.r; D = sys.D;
g0 = chan.g0;% G1 g2

traj_init_mode = 1; % circular initalization

%% initialize b
% b0 = ones(K,T);
b0 = [ones(1,T);zeros(K-1,T)];

%% Initialize p
p0 = pmax*ones(T,1);
w0 = zeros(Nt,T);
for t = 1:T
    %for k = 1:K
        w0(:,t) = sqrt(pmax)*g0(:,1,t)/norm(g0(:,1,t));
        %norm(w0(:,t))^2;
    %end
end
% w0 = zeros(Nt,K,T);
% for t = 1:T
%     Wrand = 500; % number of random precoding vector
%     W0 = randn(K*Nt,Wrand)+ 1j*randn(K*Nt,Wrand);
%     [~,minpos] = min(abs(diag(W0'*W0) - pmax));
%     w0tmp = W0(:,minpos);
%     %w0 = zeros(Nt,K);
%     for k = 1:K
%         w0(:,k,t) = w0tmp(Nt*(k-1)+1:Nt*k);
%     end
% end

%% UAV trajectory
v0 = zeros(3,T); dd = 0;
if traj_init_mode == 1 % circular initalization
    ru_vec = zeros(K,1);
    cg = sum(u,2)/K;
    for k = 1:K
        ru_vec(k) = norm(u(:,k) - cg);
    end
    ru = max(ru_vec);
    if T <= 50
        rcp = ru;
    else
        rcp = 1.5*ru;
    end
    rmax = vmax*T/(2*pi);
    rtrj = min(rmax,rcp/2);
    for t = 1:T
        theta = 2*pi*(t-1)/(T-1);
        xtrj_m = cg(1);
        ytrj_m = cg(2);
        v0x = xtrj_m + rtrj*cos(theta) - dd*dmax;
        v0y = ytrj_m + rtrj*sin(theta) - dd*dmax;
        v0(:,t) = [v0x;v0y;hv];
    end
else
    v0(:,1,:) = dmax*rand(3,1,T);
end

%% Ups0, Psi0
Ups0 = zeros(N,N,T);
if N > 0
    for t = 1:T
        Ups_rand = 0*diag(exp(1i*2*pi.*ones(N,1))); Psi0 = zeros(N,N,T);
        Ups0(:,:,t) = Ups_rand;
        for nn = 1:Na
            alpha_active = 1;
            Ups0(nn,nn,t) = sqrt(alpha_active)*Ups_rand(nn,nn);
            Psi0(nn,nn,t) = Ups_rand(nn,nn);
        end
    end
end

%% Initialize channels h0 h1 h2, effective channel h, and effective noise
[sys.sigma2,chan.h0,chan.H1,chan.h2,chan.h,chan.PL2] = update_channel(v0,Ups0,sys,chan);

%% Initialize gammatilde and Ntilde
btilde0 = 0;
if N > 0
    [cons.Q,cons.Qtilde,cons.q,cons.qtilde,cons.c4,cons.c5,cons.Xi,cons.h12tilde,btilde0] = update_Q(w0,Ups0,sys,chan);
end

%% Initialize constants c
% update_constant(Ups0,p0);
[cons.c0,cons.c1,cons.c2,cons.c3] = update_constant(Ups0,w0,sys,chan);

[v0hat0,v1hat0,v0check0,v1check0,v1tilde0] = update_v(v0,sys,cons);

%% Initialize gamma
gamma0 = zeros(K,T);
for t = 1:T
    for k = 1:K
        gamma0(k,t) = abs(chan.h(:,k,t)'*w0(:,t))^2 / sys.sigma2(k,t);
    end
end
end % EOF