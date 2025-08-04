function plot_UE_location(V,sys,plot_v,linestyle)

u = sys.u; r = sys.r; K = sys.K; dmin = sys.dmin; dmax = sys.dmax; T = sys.T;

figure
% for k = 1:K
    plot(u(1,:),u(2,:),'rp','MarkerFaceColor','r','MarkerSize',12,'LineWidth',1); hold on;
% end
plot(r(1,1),r(2,1),'gs','MarkerFaceColor','c','MarkerSize',40); hold on;


if plot_v == 1
    T = sys.T;
    color = {'ko','gd','rp','b^'};
    V1 = permute(V, [2 1]);
%     for tt = 1:T
        plot(V1([1:3:T,T],1),V1([1:3:T,T],2),linestyle,'LineWidth',2); hold on;
%     end
end
xlim([dmin dmax])
ylim([dmin dmax])
xlabel('x [m]')
ylabel('y [m]')

end % EOF
% function plot_UE_location(V)
%
% global u r K dmin dmax M T
%
% figure
% for k = 1:K
%     plot(u(1,k),u(2,k),'r*','MarkerFaceColor','r','MarkerSize',8,'LineWidth',1); hold on;
% end
% plot(r(1,1),r(2,1),'gs','MarkerFaceColor','c','MarkerSize',40); hold on;
%
% if nargin >= 1
%     color = {'ko','gd','rp','b^'};
%     for tt = 1:T
%         for m = 1:M
%             if m == 1
%                 plot(V(1,m,tt),V(2,m,tt),'bs','LineWidth',2); hold on;
%             else
%                 plot(V(1,m,tt),V(2,m,tt),'r^','LineWidth',2); hold on;
%             end
%         end
%     end
% end
% xlim([dmin dmax])
% ylim([dmin dmax])
% xlabel('x [m]')
% ylabel('y [m]')
%
% end % EOF