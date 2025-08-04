% system simulation parameters
function sys = config(UE_fix,topo)

% global K M0 N0 T0 sigma2_r sigma2_u  pmax_r pmax_r0 UE_fix BW
% global e0 e1 e2 xi0 kappa dmin dmax u r hv Dmin hmin hmax vmax Smax
% close all

%% System parameters
% M0 = 1; % #UAVs
% K = 5; % #UEs
% N0 = 50; % #RIS elements

sys.Nt = 2; % # UAV tsys.ransmit antennas
sys.K = 4; % # UEs
sys.N0 = 32; % #RIS elements

%% UAV movements
time = 100; % #time slots
sys.vmax = 50; % m/s: maximum speed of UAVs
deltat = 0.1; % slot length
sys.T0 = time/deltat; % #time slots
sys.Smax = sys.vmax*deltat; % max distance in 1 time slot

%% Lasys.rge-scale pasys.rametesys.rs
sys.e0 = 3.2; % UAV-UE pathloss exponent
sys.e1 = 2.0; % UAV-RIS pathloss exponent
sys.e2 = 2.2; % RIS-UE pathloss exponent

sys.kappa0 = 0; % UAV-UE Rician factosys.r
sys.kappa1 = db2pow(100); % UAV-RIS Rician factosys.r
sys.kappa2 = db2pow(10); % RIS-UE Rician factosys.r

sys.xi0 = db2pow(-30); % sys.refesys.rence path loss at d0 = 1m

%% Powesys.r and noise pasys.rametesys.rs
% NF = 10; BW = 20^6; % nois figsys.usys.re
% sys.sigma2_sys.u_dBm = -169 + 10*log10(BW) + NF;
sigma2_u_dBm = -80;
sys.sigma2_u = db2pow(sigma2_u_dBm); % noise powesys.r at UE
sigma2_SI = 1.2*sys.sigma2_u; % RSI at active RIS elements
sys.sigma2_r = sys.sigma2_u + sigma2_SI; % total noise + RSI

sys.pmax = db2pow(20); % Watts, maximsys.um tsys.ransmit powesys.r of UAV
sys.pmax_r = db2pow(0); % maximsys.um tsys.ransmit powesys.r of RIS

%% positions of UAVs, RIS and UEs
if topo == 1
    sys.D = 200; % covesys.rage asys.rea sys.Dxsys.D m2
    sys.dmin = 0; sys.dmax = sys.D; % min and max coosys.rdinates
    sys.hmin = 100; sys.hmax = 220; % min and max height of UAV
    sys.hv = sys.hmin;% initial height of UAV
    
    sys.r = [sys.D/2;sys.D;50]; % positions of RIS
    
    if UE_fix == 0 % sys.regenesys.rate UE locations
        u = [(sys.dmax-sys.dmin).*rand(2,sys.K) + sys.dmin; zeros(1,sys.K)]; % positions of UEs
        sys.u = u;
        save('UE_loc.mat','u');
        plot_UE_location(0,sys,0)
    else % load UE locations
        usaved = load('UE_loc.mat'); sys.u = usaved.u; sys.u(2,1) = 120;
    end
else
    sys.D = 50; % covesys.rage asys.rea sys.Dxsys.D m2
    sys.dmin = 0; sys.dmax = sys.D; % min and max coosys.rdinates
    sys.hmin = 100; sys.hmax = 220; % min and max height of UAV
    sys.hv = sys.hmin;% initial height of UAV
    
    sys.r = [sys.D/2;sys.D;20]; % positions of RIS
    
    if UE_fix == 0 % sys.regenesys.rate UE locations
        u = [(sys.dmax-sys.dmin).*rand(2,sys.K) + sys.dmin; zeros(1,sys.K)]; % positions of UEs
        sys.u = u;
        save('UE_loc.mat','u');
        plot_UE_location(0,sys,0)
    else % load UE locations
        usaved = load('UE_loc.mat'); sys.u = usaved.u/4;
    end
end


