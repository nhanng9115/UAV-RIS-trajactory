function [v] = Rankone(V,p,b)
%INPUT:
% V  =(L*N1) x (L*N1) x (K) matrix

%OUTPUT
%W   =L*N1 x K matrix

%Using Gaussian Randomization Technology

N1=size(V,1); N = N1-1;
% v1 = zeros(N1,1);

%% Gaussian Randomization
%check if V is a rank one matrix
[U, S] = eig(V);  % U: eigvector;   S: eigvalue
us = U*sqrtm(S);
% Index = find(diag(S)>10^(-4));
% if length(Index)==1  %rank one
%     v1=sqrt(S(Index,Index))*U(:,Index);
% else
%     r = 1/sqrt(2)*(randn(N1,1) + 1i*randn(N1,1));
%     v1=U*sqrtm(S)*r;  %Gaussian Randomization
% end
% 
% % v = v1(1:N)./abs(v1(1:N));
% 
% v = exp(1i*angle(v1(1:N)./v1(N+1)));

minrate0 = -Inf;
for dd = 1:50
    r = 1/sqrt(2)*(randn(N1,1) + 1i*randn(N1,1));
    v1 = us*r;  %Gaussian Randomization
    v_tmp = exp(1i*angle(v1(1:N)./v1(N+1)));
    Ups = diag(v_tmp');
    %update_channel(v_old,Ups);
    minrate = compute_rate(Ups,p,b);
    if minrate > minrate0
        v = v_tmp;
        minrate0 = minrate;
%         disp('IMPOVED>>>>>>>>>>>>>>')
    end
end

%% Maximum eigvalue
%%%%%%check if Q(:,:,k) is a rank one matrix%%%%%%
% [U, S]=eig(V(:,:));  % S: eigvector;   D %eigvalue
% [value,Index]=sort(diag(S));
% v1=sqrt(S(Index(length(diag(S))),Index(length(diag(S)))))*U(:,Index(length(diag(S))));
% v = v1(1:N);%./abs(v1(1:N));
% 
% 
% [U,S,V] = svd(V);
% u = U(:,1);
% v1 = -sqrt(S(1,1))*u;
% v = v1(1:N)./abs(v1(1:N));