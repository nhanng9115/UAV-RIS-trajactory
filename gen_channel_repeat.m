
% function gen_channel()% % clear all
clc
config

global K M0 N0 T0 kappa
M = M0; T = T0; N = N0;
kappa0 = 5;

data_location = strcat('./channels');
if ~exist(data_location, 'dir')
    mkdir(data_location)
end

n_chan = 20;

for ii = 1:n_chan
    
    % UAV-UE
    channel0 = 0; % Rician
    if channel0 == 0 % Rayleigh
        g0 = 1/sqrt(2)*(randn(M,K) + 1i*randn(M,K)); % Rayleigh fading for NLOS
    else
        g0_LoS = sqrt(kappa0/(kappa0+1))*ones(M,K,T);
        g0_NLoS = sqrt(1/(kappa0+1))*(1/sqrt(2)*(randn(M,K,T) + 1i*randn(M,K,T)));
        g0 = g0_LoS + g0_NLoS;% Rician fading for NLOS
    end
    g0 = repmat(g0,1,1,T);
    
    % UAV-RIS
    g1 = zeros(N,M);
    for m = 1:M
        kappa_LoS = 1e10;
        g1_NLoS = sqrt(1/(kappa_LoS+1))*1/sqrt(2)*(randn(N,1) + 1i*randn(N,1));
        g1_LoS = sqrt(kappa_LoS/(kappa_LoS+1))*LoS_channel(1, N, 'UAV-RIS');
        g1(:,m) = g1_LoS + g1_NLoS;% Rician fading for NLOS
    end
    g1 = repmat(g1,1,1,T);
    
    % RIS-UE
    g2 = zeros(N,K);
    for k = 1:K
        g2_NLoS = sqrt(1/(kappa+1))*1/sqrt(2)*(randn(1,N) + 1i*randn(1,N));
        g2_LoS = sqrt(kappa/(kappa+1))*LoS_channel(N, 1, 'RIS-UE');
        g2(:,k) = (g2_LoS + g2_NLoS).'; % RIS-AP
    end
    g2 = repmat(g2,1,1,T);
    
    % save large scale params
    file_name = strcat(data_location,'/small_scale',num2str(ii));
    save(file_name,'g0','g1','g2');
end

% end % EOF