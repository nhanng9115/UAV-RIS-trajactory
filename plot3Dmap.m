%% plot deployment map
close all
global r u
config

%% Load data
d1 = load('power','Pall'); Pall = d1.Pall;
d2 = load('loc','Vall'); Vall = d2.Vall; 


%% Plot UAV locations
marker = {'o','d','p'};
figure
n_chan = 10;

%% x-y coordinate
plot(r(1),r(2),'s','MarkerSize',30,'MarkerFaceColor','c'); hold on
plot(u(1,1),u(2,1),'^','MarkerSize',15,'MarkerFaceColor','r'); hold on;

for cc = 1:n_chan
    for rr = 0:2
        vv = Vall(:,cc,rr+1);
        h = plot(vv(1),vv(2),marker{rr+1},'MarkerSize',8,'LineWidth',2); hold on;
    end
end
for k = 2:K
    plot(u(1,k),u(2,k),'^','MarkerSize',15,'MarkerFaceColor','r'); hold on;
end

xlim([0 200]); ylim([0 200]);
xlabel('x [m]'); ylabel('y [m]');

%% y-z coordinate
figure
plot(r(1),r(3),'s','MarkerSize',30,'MarkerFaceColor','c'); hold on
plot(u(1,1),u(3,1),'^','MarkerSize',15,'MarkerFaceColor','r'); hold on;


for cc = 1:n_chan
    for rr = 0:2
        vv = Vall(:,cc,rr+1);
        h = plot(vv(1),vv(3),marker{rr+1},'MarkerSize',8,'LineWidth',2); hold on;
    end
end
for k = 2:K
    plot(u(1,k),u(3,k),'^','MarkerSize',15,'MarkerFaceColor','r'); hold on;
end

%% Plot power allocation
figure
PP = mean(Pall,2); P = permute(PP,[1 3 2]); x = [12 100 174 53];
bar(x,P,'group','LineWidth',2); hold on;
xlabel('UEs'); ylabel('Averaged allocated power [mW]');
xticklabels({'UE 1','UE 2','UE 3','UE 4'})
legend('W/o RIS','Passive RIS','Hybrid RIS');
% marker = {'o','d','p'};
% scatter3(r(1),r(2),r(3),5000,'sc','filled'); hold on
% 
% for k = 1:K
%     scatter3(u(1,:),u(2,:),u(3,:),200,'rs','filled'); hold on;
% end
% 
% for cc = 1:n_chan
%     for rr = 0:2
%         vv = Vall(:,cc,rr+1);
%         scatter3(vv(1,:),vv(2,:),vv(3,:),120,marker{rr+1},'filled'); hold on;
%     end
% end
% 
% PP = mean(Pall,2); P = permute(PP,[1 3 2]);
% bar3(P,'grouped')
% title('Grouped Style')