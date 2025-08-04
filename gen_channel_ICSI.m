
% function gen_channel()% % clear all
clc
sys = config(1);

K = sys.K; N0 = sys.N0; T0 = sys.T0; kappa0 = sys.kappa0; kappa1 = sys.kappa1; kappa2 = sys.kappa2;
Nt = sys.Nt; T = T0; N = N0;

data_location = strcat('./channels_ICSI/',num2str(sys.N0));
if ~exist(data_location, 'dir')
    mkdir(data_location)
end

n_chan = 10;

for ii = 1:n_chan
    
    % UAV-UE
    channel0 = 0; % Rician
    if channel0 == 0 % Rayleigh
        g0 = 1/sqrt(2)*(randn(Nt,K,T) + 1i*randn(Nt,K,T)); % Rayleigh fading for NLOS
    else
        g0_LoS = sqrt(kappa0/(kappa0+1))*ones(Nt,K,T);
        g0_NLoS = sqrt(1/(kappa0+1))*(1/sqrt(2)*(randn(Nt,K,T) + 1i*randn(Nt,K,T)));
        g0 = g0_LoS + g0_NLoS;% Rician fading for NLOS
    end
    % add error
    g0 = g0 - eps*1/sqrt(2)*(randn(Nt,K,T) + 1i*randn(Nt,K,T));
    
    
    % UAV-RIS
    G1 = zeros(N,Nt,T);
    for t = 1:T
        kappa1 = 1e10;
        g1_NLoS = sqrt(1/(kappa1+1))*1/sqrt(2)*(randn(N,Nt) + 1i*randn(N,Nt));
        g1_LoS = sqrt(kappa1/(kappa1+1))*LoS_channel(Nt, N, 'UAV-RIS');
        G1(:,:,t) = g1_LoS + g1_NLoS;% Rician fading for NLOS
    end
    
    % RIS-UE
    g2 = zeros(N,K,T);
    for t = 1:T
        for k = 1:K
            g2_NLoS = sqrt(1/(kappa2+1))*1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
            g2_LoS = sqrt(kappa2/(kappa2+1))*LoS_channel(N, 1, 'RIS-UE');
            g2(:,k,t) = (g2_LoS + g2_NLoS).'; % RIS-AP
        end
    end
    
    % save large scale params
    file_name = strcat(data_location,'/small_scale',num2str(ii));
    save(file_name,'g0','G1','g2');
end

% end % EOF