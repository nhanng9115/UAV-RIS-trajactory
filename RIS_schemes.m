function [obj_RIS,Ups_new,btilde_new] = RIS_schemes(RIS,Ups_old,w_old,v_old,btilde_old,b_old,sys,chan,cons)

N = sys.N; Na = sys.Na; T = sys.T;

if RIS == 1 % Passive RIS
    [Phi_new] = opt_phase(b_old,w_old,Ups_old,sys,cons,chan);
    obj_RIS = 0; Ups_new = Phi_new; btilde_new = 0;
    
else%if RIS == 2 || RIS == 3 || RIS == 4 % Hybrid RIS
    
    [obj_RIS,Psi_new,btilde_new] = opt_actRIS(Ups_old,btilde_old,b_old,sys,cons,chan);
    Psi = zeros(N,N,T); Psi(1:Na,1:Na,:) = Psi_new(1:Na,1:Na,:);
    
    [Phi_new] = opt_phase(b_old,w_old,Psi,sys,cons,chan);
    Phi = zeros(N,N,T); Phi(Na+1:N,Na+1:N,:) = Phi_new(Na+1:N,Na+1:N,:);

    Ups_new = Psi + Phi;
end

end % EOF