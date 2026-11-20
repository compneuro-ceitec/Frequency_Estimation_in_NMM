%% SETTING (GRAPHIC STUFF)
Path = 'C:\prace\MatCont7p6\Systems';

msty = 'k.';
ms = 12; % 4;
ggreen = [0 0.8 0];
ggray = [0.6 0.6 0.6];
lsty_unst_LC = 'r--';
lsty_st_LC = 'b';
lsty_LP = 'b';
clr_H = ggreen;
lsty_Hsub = '--';
lsty_Hsup = '-';
lsty_LPC = 'r-.';
lsty_Hom = 'k';
lwd = 1;
lwd_orb = 0.5;

lsty_line = 'k:';
lwd_line = 0.25;

lsty_LPC_IPS1 = 'b';
lsty_LPC_IPS2 = 'b-.';
lsty_LPC_APS = 'k';
lsty_LPC_outer = 'r--';
lsty_LPC_STO = '-.';
clr_LPC_STO = [0.6 0.6 0.6];
clr_NS = clr_H;
lsty_NS_sub = '--';

set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
h = figure; pos = get(h,'Position'); close(h) % default plot position


%% FIG. 1: ORBIT (BURSTING)
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\P_bursting_Md_0.668.mat'])

idx = t < 1400; t = t(idx);

r = y(idx,1);
V = y(idx,2);
u = y(idx,3);
Sa = y(idx,4);

rlims = [-0.01 0.14]; rticks = 0:0.07:0.14;
Vlims = [-62 -43]; Vticks = -60:8:-44;
ulims = [-90 110]; uticks = -70:90:110;
Salims = [-0.5 11]; Saticks = 0:5:10;
tlims = [0 1400]; tticks = 0:200:1400;

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.5])
tlo = tiledlayout(12,5,'TileSpacing','tight','Padding','tight');

%
nexttile([3,3])
plot(t,r,'b')
ylabel('$r$'), ylim(rlims), yticks(rticks)
xlim(tlims), xticks(tticks), xticklabels({})
box off

nexttile([3,3])
plot(t,V,'b')
ylabel('$v$'), ylim(Vlims), yticks(Vticks)
xlim(tlims), xticks(tticks), xticklabels({})
box off

nexttile([3,3])
plot(t,u,'r')
ylabel('$u$'), ylim(ulims), yticks(uticks)
xlim(tlims), xticks(tticks), xticklabels({})
box off

nexttile([3,3])
plot(t,Sa,'r')
ylabel('$s_a$'), ylim(Salims), yticks(Saticks)
xlabel('$t$ [ms]'), xlim(tlims), xticks(tticks)
box off

%
rlims = [-0.02 0.14]; 
Vlims = [-62 -23];
ulims = [-110 110];
N = length(r);

nexttile([8,2])
plot3(u,V,r,'b','LineWidth',0.5), hold on

plot3(repelem(ulims(1),N),V,r,'Color',ggray)
plot3(u,repelem(Vlims(2),N),r,'Color',ggray)
plot3(u,V,repelem(rlims(1),N),'Color',ggray)

xlabel('$u$'), xlim(ulims), xticks(uticks), xtickangle(0)
ylabel('$v$'), ylim(Vlims), yticks(Vticks), ytickangle(0)
zlabel('$r$'), zlim(rlims), zticks(rticks), ztickangle(0)

view(51,12)

%
nexttile([4,1])
plot(V,r,'b')
xlabel('$v$'), xlim([-62 -43]), xticks(Vticks), xtickangle(0)
ylabel('$r$'), ylim([-0.01 0.14]), yticks(rticks), ytickangle(0)
box off

%
nexttile([4,1])
plot(Sa,u,'r')
xlabel('$s_a$'), xlim(Salims), xticks(Saticks), xtickangle(0)
ylabel('$u$'), ylim([-90 110]), yticks(uticks), ytickangle(0)
box off





%% FIG. 2: ORBIT (SPIKING)
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\P_spiking_Md_5.mat'])

idx = t < 100; t = t(idx);

r = y(idx,1);
V = y(idx,2);
u = y(idx,3);
Sa = y(idx,4);

rlims = [0 1.05]; rticks = 0:0.5:1.1;
Vlims = [-65 80]; Vticks = -65:65:65;
ulims = [660 680]; uticks = 660:10:680;
Salims = [24 30.5]; Saticks = 25:2.5:30;
tlims = [0 70]; tticks = 0:10:70;

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.5])
tlo = tiledlayout(12,5,'TileSpacing','tight','Padding','tight');

nexttile([3,3])
plot(t,r,'b')
ylabel('$r$'), ylim(rlims), yticks(rticks)
xlim(tlims), xticks(tticks), xticklabels({})
box off

nexttile([3,3])
plot(t,V,'b')
ylabel('$v$'), ylim(Vlims), yticks(Vticks)
xlim(tlims), xticks(tticks), xticklabels({})
box off

nexttile([3,3])
plot(t,u,'r')
ylabel('$u$'), ylim(ulims), yticks(uticks)
xlim(tlims), xticks(tticks), xticklabels({})
box off

nexttile([3,3])
plot(t,Sa,'r')
ylabel('$s_a$'), ylim(Salims), yticks(Saticks)
xlabel('$t$ [ms]'), xlim(tlims), xticks(tticks)
box off

%
rlims = [-0.1 1.05];
Vlims = [-65 100];
ulims = [640 680];
N = length(r);

nexttile([8,2])
plot3(u,V,r,'b','LineWidth',0.5), hold on

plot3(repelem(ulims(1),N),V,r,'Color',ggray)
plot3(u,repelem(Vlims(2),N),r,'Color',ggray)
plot3(u,V,repelem(rlims(1),N),'Color',ggray)

xlabel('$u$'), xlim(ulims), xticks(uticks), xtickangle(0)
ylabel('$v$'), ylim(Vlims), yticks(Vticks), ytickangle(0)
zlabel('$r$'), zlim(rlims), zticks(rticks), ztickangle(0)

view(44,12)

%
nexttile([4,1])
plot(V,r,'b')
xlabel('$v$'), xlim([-65 80]), xticks(Vticks), xtickangle(0)
ylabel('$r$'), ylim([0,1.05]), yticks(rticks), ytickangle(0)
box off

%
nexttile([4,1])
plot(Sa,u,'r')
xlabel('$s_a$'), xlim(Salims), xticks(Saticks), xtickangle(0)
ylabel('$u$'), ylim([660,680]), yticks(uticks), ytickangle(0)
box off



%% FIG. 8a: Delta x Md [left]
fig = figure; hold on
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4.96,4.8])

% Hopf
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1a.mat'])
Hx1 = x(end-2,:);
Hy1 = x(end-1,:);
GH = 118;

load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
Hx2 = x(end-2,:);
Hy2 = x(end-1,:);

plot([flip(Hx1(1:GH)),Hx2],[flip(Hy1(1:GH)),Hy2],'LineStyle',lsty_Hsub,'Color',clr_H,'LineWidth',lwd);
plot(Hx1(GH:end),Hy1(GH:end),'LineStyle',lsty_Hsup,'Color',clr_H,'LineWidth',lwd);

plot(Hx1(GH),Hy1(GH),msty,'MarkerSize',ms)

% LPC
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1a.mat'])
LPCx = x(end-1,:);
LPCy = x(end,:);

load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1b.mat'])
LPCx = [LPCx, x(end-1,:)];
LPCy = [LPCy, x(end,:)];

plot(LPCx,LPCy,lsty_LPC,'LineWidth',lwd);

% Default values
plot(0.668,0.5,'ko','MarkerSize',4), plot(0.668,0.5,'k.','MarkerSize',6)

% Axes
xlabel('$M_{D1}$'); ylabel('$\Delta_{v_{\theta}}$');
xlims = [0,0.8]; xticks(0:0.1:0.8);
ylims = [0,5.5]; yticks([0,0.5,1:6]);
set(gca,'DataAspectRatio',[range(xlims),range(ylims),1],'box','on')
xlim(xlims),ylim(ylims),xtickangle(0)
yline(0.5,lsty_line,'LineWidth',lwd_line)





%% FIG. 8b: Delta x Md [right]
fig = figure; hold on
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[4.96,4.8])

% Hopf
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
Hx1 = x(end-2,:);
Hy1 = x(end-1,:);
GH = 1503;

plot(Hx1(1:GH),Hy1(1:GH),'LineStyle',lsty_Hsub,'Color',clr_H,'LineWidth',lwd);
plot(Hx1(GH:end),Hy1(GH:end),'LineStyle',lsty_Hsup,'Color',clr_H,'LineWidth',lwd);

plot(Hx1(GH),Hy1(GH),msty,'MarkerSize',ms)

% LPC
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2a.mat'])
LPCx = x(end-1,:);
LPCy = x(end,:);

load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2b.mat'])
LPCx = [flip(LPCx), x(end-1,:)];
LPCy = [flip(LPCy), x(end,:)];

plot(LPCx,LPCy,lsty_LPC,'LineWidth',lwd);

% Default values
plot(5,0.5,'ko','MarkerSize',4), plot(5,0.5,'k.','MarkerSize',6)

% Axes
xlabel('$M_{D1}$'); ylabel('$\Delta_{v_{\theta}}$');
xlims = [2.5,20]; xticks(2.5:2.5:20);
ylims = [0,5.5]; yticks([0,0.5,1:6]);
set(gca,'DataAspectRatio',[range(xlims),range(ylims),1],'box','on')
xlim(xlims),ylim(ylims),xtickangle(0)
yline(0.5,lsty_line,'LineWidth',lwd_line)


% =========================================================
% INSET AXES (top-left)
% =========================================================

ax_main = gca;

ax_inset = axes( ...
    'Position',[0.26 0.61 0.23 0.23]); % [left bottom width height]

uistack(ax_inset,'top')
hold(ax_inset,'on')

% Hopf
plot(ax_inset,Hx1(1:GH),Hy1(1:GH), ...
    'LineStyle',lsty_Hsub, 'Color',clr_H, 'LineWidth',lwd);

plot(ax_inset,Hx1(GH:end),Hy1(GH:end), ...
    'LineStyle',lsty_Hsup, 'Color',clr_H, 'LineWidth',lwd);

plot(ax_inset,Hx1(GH),Hy1(GH), ...
    msty,'MarkerSize',ms)

% LPC
plot(ax_inset,LPCx,LPCy, lsty_LPC,'LineWidth',lwd);

% Default values
% plot(5,0.5,'ko','MarkerSize',4), plot(5,0.5,'k.','MarkerSize',6)

% Limits for zoomed region
xlim(ax_inset,[4 5])
ylim(ax_inset,[0 0.5])

box(ax_inset,'on')

% Cleaner inset appearance
set(ax_inset, ...
    'FontSize',8, ...
    'LineWidth',0.8)

% Optional:
xticks(ax_inset,[4:0.5:5])
yticks(ax_inset,[0:0.25:0.5])



%% FIG. 9: LC manifold w.r.t. Md [left]
% LC
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LC_Md_Delta_0.5_left.mat'])
LPC = 78;

V = x(2:4:end-4,1:125);
r = x(1:4:end-5,1:125);
Md = x(end,1:125);
freq = 1000 ./ x(end-1,1:125);

% EP
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\EP_Md_Delta_0.5_b.mat'])
H1 = 295; H2 = 572;

V2 = x(2,1:800);
r2 = x(1,1:800);
Md2 = x(end,1:800);

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.5])
tiledlayout(2,4,'TileSpacing','compact','Padding','compact')


% --- [Md,v,r] ------------------------------------------------------------
nexttile([2,3]), hold on

% =================================
% CREATE SURFACE FROM LIMIT CYCLES
% =================================

% Each column of:
%   V(:,i), r(:,i)
% corresponds to one limit cycle at parameter Md(i).
%
% We construct a surface:
%   rows    -> phase along LC
%   columns -> continuation step

% LC (dense continuation)
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LC_Md_Delta_0.5_left_dense.mat'])
LPC3 = 623;

V3 = x(2:4:end-4,:);
r3 = x(1:4:end-5,:);
Md3 = x(end,:);
freq3 = 1000 ./ x(end-1,:);

Xsurf = repmat(Md3(:)',size(V3,1),1);
Ysurf = V3;
Zsurf = r3;

% --- Stable LC surface
surf(Xsurf(:,LPC3:end), Ysurf(:,LPC3:end), Zsurf(:,LPC3:end), ...
    'FaceColor',[0 0 1], 'EdgeColor','none', 'FaceAlpha',0.03);

% --- Unstable LC surface
surf(Xsurf(:,1:LPC3-1), Ysurf(:,1:LPC3-1), Zsurf(:,1:LPC3-1), ...
    'FaceColor',[1 0 0], 'EdgeColor','none', 'FaceAlpha',0.08);

% ================
% ORIGINAL CURVES
% ================

% Stable LCs
for i = LPC:3:length(Md)
    plot3(Md(i)*ones(1,size(r,1)),V(:,i),r(:,i),'b-','LineWidth',0.5)
end

% Unstable LCs
for i = 1:3:LPC-12
    plot3(Md(i)*ones(1,size(r,1)),V(:,i),r(:,i),'r--','LineWidth',0.5)
end

% LPC
plot3(Md(LPC)*ones(1,size(r,1)),V(:,LPC),r(:,LPC),'r-','LineWidth',1.5)

% Stable EPs
plot3(Md2(1:H1),V2(1:H1),r2(1:H1),'-','Color',ggreen,'LineWidth',1)

plot3(Md2(H2:end),V2(H2:end),r2(H2:end),'-','Color',ggreen,'LineWidth',1)

% Unstable EPs
plot3(Md2(H1:H2),V2(H1:H2),r2(H1:H2),'--','Color',ggreen,'LineWidth',1)

% Hopf points
plot3(Md2([H1,H2]),V2([H1,H2]),r2([H1,H2]),msty,'MarkerSize',ms)

% Axes
xlabel('$M_{D1}$'); ylabel('$v$'); zlabel('$r$')

xlims = [0,1]; xticks(0:0.2:1);
ylims = [-60.5,-43]; yticks(-60:5:-40);
zlims = [0,0.16]; zticks(0:0.04:0.16);
xlim(xlims), ylim(ylims), zlim(zlims),ytickangle(0)

view([13,7])

% Better transparency rendering
set(gcf,'Renderer','opengl')

camlight headlight
lighting gouraud


% --- [Md,r] --------------------------------------------------------------
nexttile

% LC (dense continuation)
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LC_Md_Delta_0.5_left_dense.mat'])
LPC3 = 623;

V3 = x(2:4:end-4,:);
r3 = x(1:4:end-5,:);
Md3 = x(end,:);
freq3 = 1000 ./ x(end-1,:);

Xsurf = repmat(Md3(:)',size(V3,1),1);
Ysurf = V3;
Zsurf = r3;

% --- Stable LC surface
surf(Xsurf(:,LPC3:end), Ysurf(:,LPC3:end), Zsurf(:,LPC3:end), ...
    'FaceColor',[0 0 1], 'EdgeColor','none', 'FaceAlpha',0.08); hold on

% --- Unstable LC surface
surf(Xsurf(:,1:LPC3-1), Ysurf(:,1:LPC3-1), Zsurf(:,1:LPC3-1), ...
    'FaceColor',[1 0 0], 'EdgeColor','none', 'FaceAlpha',0.08);

grid off

% % Stable LCs
% r_min = min(r,[],1); r_max = max(r,[],1);
% plot(Md(LPC:end),r_min(LPC:end),'b-','LineWidth',1)%, hold on
% plot(Md(LPC:end),r_max(LPC:end),'b-','LineWidth',1)
% 
% % Unstable LCs
% plot(Md(1:LPC),r_min(1:LPC),'r--','LineWidth',1)
% plot(Md(1:LPC),r_max(1:LPC),'r--','LineWidth',1)
% 
% % LPC
% plot(Md(LPC),r_min(LPC),msty,'MarkerSize',ms)
% plot(Md(LPC),r_max(LPC),msty,'MarkerSize',ms)
% 
% % Stable EPs
% plot(Md2(1:H1),r2(1:H1),'-','Color',ggreen,'LineWidth',1)
% plot(Md2(H2:end),r2(H2:end),'-','Color',ggreen,'LineWidth',1)
% 
% % Unstable EPs
% plot(Md2(H1:H2),r2(H1:H2),'--','Color',ggreen,'LineWidth',1)
% 
% % Hopf points
% plot(Md2([H1,H2]),r2([H1,H2]),msty,'MarkerSize',ms)
% 
% % Axes
% xlabel('$M_{D1}$'); ylabel('$r$');
% xlims = [0,1]; xticks(0:0.2:1);
% ylims = [-0.005,0.16]; yticks(0:0.04:0.16);
% 
% xlim(xlims),ylim(ylims),xtickangle(0)

% Axes
xlabel('$M_{D1}$'); ylabel(''); zlabel('$r$')

xlims = [0,1]; xticks(0:0.2:1);
ylims = [-100,-43]; yticks([]);
zlims = [-0.005,0.16]; zticks(0:0.04:0.16);
xlim(xlims), ylim(ylims), zlim(zlims), xtickangle(0)

view([0,0])

% Better transparency rendering
set(gcf,'Renderer','opengl')

camlight headlight
lighting gouraud

% --- [Md,frequency] ------------------------------------------------------
nexttile

% Stable LCs
plot(Md(LPC:end),freq(LPC:end),'b-','LineWidth',1), hold on

% Unstable LCs
plot(Md(1:LPC),freq(1:LPC),'r--','LineWidth',1)

% LPC
plot(Md(LPC),freq(LPC),msty,'MarkerSize',ms)

% Stable EPs
plot(Md2([1,H1]),[0,0],'-','Color',ggreen,'LineWidth',1)
plot(Md2([H2,end]),[0,0],'-','Color',ggreen,'LineWidth',1)

% Unstable EPs
plot(Md2([H1,H2]),[0,0],'--','Color',ggreen,'LineWidth',1)

% Hopf points
plot(Md2([H1,H2]),[0,0],msty,'MarkerSize',ms)
plot(Md2([H1,H2]),freq([1,end]),msty,'MarkerSize',ms)

% Axes
xlabel('$M_{D1}$'); ylabel('Interburst frequency [Hz]');
xlims = [0,1]; xticks(0:0.2:1);
ylims = [-0.1,7]; yticks(0:2:6);

xlim(xlims),ylim(ylims),xtickangle(0)





%% FIG. 10: LC manifold w.r.t. Md [right]
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LC_Md_Delta_0.5_right.mat'])

V = x(2:4:end-4,[1:3:124,124:2:223]);
r = x(1:4:end-5,[1:3:124,124:2:223]);
Md = x(end,[1:3:124,124:2:223]);

V_dense = x(2:4:end-4,:);
Md_dense = x(end,:);
r_dense = x(1:4:end-5,:);
freq = 1000 ./ x(end-1,:);

% EP
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\EP_Md_Delta_0.5_a.mat'])
H1 = 962; 

V2 = x(2,1:2000);
r2 = x(1,1:2000);
Md2 = x(end,1:2000);

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.5])
tiledlayout(2,4,'TileSpacing','compact','Padding','compact')


% --- [Md,v,r] ------------------------------------------------------------
nexttile([2,3]), hold on

% =================================
% CREATE SURFACE FROM LIMIT CYCLES
% =================================

% Surface coordinates:
%   rows    -> phase along LC
%   columns -> continuation step

Xsurf = repmat(Md(:)',size(V,1),1);
Ysurf = V;
Zsurf = r;

% --- Transparent LC manifold surface
surf(Xsurf, Ysurf, Zsurf, ...
    'FaceColor',[0 0 1], 'EdgeColor','none', 'FaceAlpha',0.10);

% =================
% ORIGINAL CURVES
% =================

% Stable LCs
for i = 1:4:length(Md)
    plot3(Md(i)*ones(1,size(r,1)),V(:,i),r(:,i),'b-','LineWidth',0.5)
end

% Stable EPs
plot3(Md2(1:H1),V2(1:H1),r2(1:H1),'-','Color',ggreen,'LineWidth',1)

% Unstable EPs
plot3(Md2(H1:end),V2(H1:end),r2(H1:end),'--','Color',ggreen,'LineWidth',1)

% Hopf points
plot3(Md2(H1),V2(H1),r2(H1),msty,'MarkerSize',ms)

% Axes
xlabel('$M_{D1}$'); ylabel('$v$'); zlabel('$r$')

xlims = [4.6,5.66]; xticks(4.6:0.2:5.66);
ylims = [-130,175]; yticks(-130:100:175);
zlims = [0,2.25]; zticks(0:0.5:2);

xlim(xlims), ylim(ylims), zlim(zlims), ytickangle(0)

view([-15,7])

% Better transparency rendering
set(gcf,'Renderer','opengl')

camlight headlight
lighting gouraud


% --- [Md,r] --------------------------------------------------------------
nexttile

Xsurf = repmat(Md_dense(:)',size(V,1),1);
Ysurf = V_dense;
Zsurf = r_dense;

% --- Transparent LC manifold surface
surf(Xsurf, Ysurf, Zsurf, ...
    'FaceColor',[0 0 1], 'EdgeColor','none', 'FaceAlpha',0.10); hold on

grid off

% Better transparency rendering
set(gcf,'Renderer','opengl')

camlight headlight
lighting gouraud

% % ================
% % ORIGINAL CURVES
% % ================
% 
% % Stable LCs
% r_min = min(r_dense,[],1); r_max = max(r_dense,[],1);
% v_aux = repelem(-100,length(Md_dense));
% plot3(Md_dense,v_aux,r_min,'b-','LineWidth',1)
% plot3(Md_dense,v_aux,r_max,'b-','LineWidth',1)
% 
% % Stable EPs
% plot3(Md2(1:H1),V2(1:H1),r2(1:H1),'-','Color',ggreen,'LineWidth',1)
% 
% % Unstable EPs
% plot3(Md2(H1:end),V2(H1:end),r2(H1:end),'--','Color',ggreen,'LineWidth',1)
% 
% % Hopf points
% plot3(Md2(H1),V2(H1),r2(H1),msty,'MarkerSize',ms)

% Axes
xlabel('$M_{D1}$'); ylabel(''); zlabel('$r$');
xlims = [4.6,5.66]; xticks(4.6:0.2:5.66);
ylims = [-130,175]; yticks([]);
zlims = [-0.05,2.25]; zticks(0:0.5:2);

xlim(xlims), ylim(ylims), zlim(zlims), xtickangle(0)

view([0,0])


% --- [Md,frequency] ------------------------------------------------------
nexttile

% Stable LCs
plot(Md_dense,freq,'b-','LineWidth',1), hold on

% Stable EPs
plot(Md2([1,H1]),[0,0],'-','Color',ggreen,'LineWidth',1)

% Unstable EPs
plot(Md2([H1,end]),[0,0],'--','Color',ggreen,'LineWidth',1)

% Hopf points
plot(Md2(H1),0,msty,'MarkerSize',ms)
plot(Md2(H1),freq(1),msty,'MarkerSize',ms)


% Axes
xlabel('$M_{D1}$'); ylabel('Oscillation frequency [Hz]');
xlims = [4.6,5.66]; xticks(4.6:0.2:5.66);
ylims = [-10,350]; yticks(0:100:350);

xlim(xlims),ylim(ylims),xtickangle(0)


%% NOT USED: FREQUENCY ESTIMATE (single method)
lwd = 1.5;

% Data loading ------------------------------------------------------------
Md_max = 10;
Md_max = 20;

load(['freq_data_FP_',num2str(Md_max),'.mat']); 
load(['freq_data_LS_',num2str(Md_max),'.mat']);
load(['freq_data_EP_',num2str(Md_max),'.mat']); 
% load(['freq_data_LC_',num2str(Md_max),'.mat']);

% Left bubble -------------------------------------------------------------
L_Md_vals    = data.L_freq(1,:,:);
L_Delta_vals = data.L_freq(2,:,:);
L_freq       = data.L_freq(3,:,:);
% L_freq(L_freq==0) = NaN;

L_Md_from = data.L_bubble.Md_from;
L_Md_to = data.L_bubble.Md_to;
L_Md = [L_Md_from,flip(L_Md_to)];

L_Delta_from = data.L_bubble.Delta_from;
L_Delta_to = data.L_bubble.Delta_to;
L_Delta = [L_Delta_from,flip(L_Delta_to)];

% Right bubble ------------------------------------------------------------
R_Md_vals    = data.R_freq(1,:,:);
R_Delta_vals = data.R_freq(2,:,:);
R_freq       = data.R_freq(3,:,:);
% R_freq(R_freq==0) = NaN;

R_Md_from = data.R_bubble.Md_from;
R_Md_to = data.R_bubble.Md_to; 
ind = R_Md_to < max(R_Md_to); R_Md_to = R_Md_to(ind);
R_Md = [R_Md_from,flip(R_Md_to)];

R_Delta_from = data.R_bubble.Delta_from;
R_Delta_to = data.R_bubble.Delta_to; R_Delta_to = R_Delta_to(ind);
R_Delta = [R_Delta_from,flip(R_Delta_to)];


% Plotting ================================================================
ms = 3;             % marker size
cmap = parula(256); %[[1, 1, 1]; [1, 1, 1]; parula(254)];

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[10,6.5])
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');


% BURSTING ----------------------------------------------------------------
ax1 = nexttile([1,1]);
scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq(:), 'filled'), hold on
% plot(L_Md,L_Delta,'k-','LineWidth',1)

% Colormap
colormap(cmap)
% freq_lim = [floor(min(L_freq(:),[],'omitnan')) ...
%             ceil(max(L_freq(:),[],'omitnan'))];
freq_lim = [0 8];
cb = colorbar; clim(freq_lim);
cb.Label.String = 'Interburst frequency [Hz]';
cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';

% Axes
xlabel('$M_{D1}$', 'Interpreter', 'latex');
ylabel('$\Delta_{v_{\theta}}$', 'Interpreter', 'latex');

Md_lims = [0 0.8]; Md_ticks = 0:0.2:0.8;
Delta_lims = [0 6]; Delta_ticks = 0:1:6;
xlim(Md_lims); xticks(Md_ticks); xtickangle(0)
ylim(Delta_lims); yticks(Delta_ticks)   

xaxis = get(gca, 'XAxis'); xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis'); yaxis.TickLabelInterpreter = 'latex';
set(gca, 'Layer', 'top', 'box', 'on')


% SPIKING -----------------------------------------------------------------
ax2 = nexttile([1,1]);
scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq(:), 'filled'), hold on
% plot(R_Md,R_Delta,'k-','LineWidth',1)

% Colormap
colormap(cmap)
% freq_lim = [floor(min(R_freq,[],'all')) ceil(max(R_freq,[],'all'))];
freq_lim = [0 500];
cb = colorbar; clim(freq_lim);
cb.Label.String = 'Oscillation frequency [Hz]'; 
cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';

% Axes
xlabel('$M_{D1}$', 'Interpreter', 'latex');
ylabel('$\Delta_{v_{\theta}}$', 'Interpreter', 'latex');

Md_lims = [4 Md_max]; Md_ticks = linspace(4,Md_max,5);
Delta_lims = [0 6]; Delta_ticks = 0:1:6;
xlim(Md_lims); xticks(Md_ticks); xtickangle(0)
ylim(Delta_lims); yticks(Delta_ticks)

xaxis = get(gca, 'XAxis'); xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis'); yaxis.TickLabelInterpreter = 'latex';
set(gca, 'Layer', 'top', 'box', 'on')


% BURSTING & SPIKING ------------------------------------------------------

ax3 = nexttile([1,2]);
scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq(:), 'filled'), hold on
% plot(L_Md,L_Delta,'k-','LineWidth',1)

scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq(:), 'filled')
% plot(R_Md,R_Delta,'k-','LineWidth',1)

% Colormap
colormap(cmap)
% freq_lim = [floor(min([L_freq(:); R_freq(:)],[],'all')), ...
%             ceil(max([L_freq(:); R_freq(:)],[],'all'))];
freq_lim = [0 500];
cb = colorbar; clim(freq_lim);
cb.Label.String = 'Frequency [Hz]'; 
cb.TickLabelInterpreter = 'latex'; cb.Label.Interpreter = 'latex';
    
% Axes
xlabel('$M_{D1}$', 'Interpreter', 'latex');
ylabel('$\Delta_{v_{\theta}}$', 'Interpreter', 'latex');

Md_lims = [0 Md_max]; Md_ticks = linspace(0,Md_max,11);
Delta_lims = [0 6]; Delta_ticks = 0:1:6;
xlim(Md_lims); xticks(Md_ticks); xtickangle(0)
ylim(Delta_lims); yticks(Delta_ticks)

xaxis = get(gca, 'XAxis'); xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis'); yaxis.TickLabelInterpreter = 'latex';
set(gca, 'Layer', 'top', 'box', 'on')


% Bifurcation curves ------------------------------------------------------
% Hopf
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1a.mat'])
Hx1 = x(end-2,:);
Hy1 = x(end-1,:);
GH = 118;

load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
Hx2 = x(end-2,:);
Hy2 = x(end-1,:);

plot(ax1,[flip(Hx1(1:GH)),Hx2],[flip(Hy1(1:GH)),Hy2],'k--','LineWidth',lwd);
plot(ax3,[flip(Hx1(1:GH)),Hx2],[flip(Hy1(1:GH)),Hy2],'k--','LineWidth',lwd);
plot(ax1,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
plot(ax3,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
plot(ax1,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)
plot(ax3,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)

%
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
Hx1 = x(end-2,:);
Hy1 = x(end-1,:);
GH = 1503;

plot(ax2,Hx1(1:GH),Hy1(1:GH),'k--','LineWidth',lwd);
plot(ax3,Hx1(1:GH),Hy1(1:GH),'k--','LineWidth',lwd);
plot(ax2,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
plot(ax3,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
plot(ax2,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)
plot(ax3,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)


% LPC
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1a.mat'])
LPCx = x(end-1,:);
LPCy = x(end,:);

load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1b.mat'])
LPCx = [LPCx, x(end-1,:)];
LPCy = [LPCy, x(end,:)];

plot(ax1,LPCx,LPCy,'k-.','LineWidth',lwd);
plot(ax3,LPCx,LPCy,'k-.','LineWidth',lwd);

%
load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2a.mat'])
LPCx = x(end-1,:);
LPCy = x(end,:);

load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2b.mat'])
LPCx = [flip(LPCx), x(end-1,:)];
LPCy = [flip(LPCy), x(end,:)];

plot(ax2,LPCx,LPCy,'k-.','LineWidth',lwd);
plot(ax3,LPCx,LPCy,'k-.','LineWidth',lwd);


%% FIG. 3: FREQUENCY ESTIMATION (all methods)
% Md_max = 10;
Md_max = 20;

methods = {'FP','LS','EP'};
n_methods = numel(methods);

ms = 4;
cmap = parula(256);
lwd = 1.5;

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[11,6.25])
tlo = tiledlayout(2,3,'TileSpacing','compact','Padding','tight');

for i = 1:n_methods
    
    load(['freq_data_',methods{i},'_',num2str(Md_max),'.mat']);
    
    % --- Left bubble ---
    L_Md_vals    = data.L_freq(1,:,:);
    L_Delta_vals = data.L_freq(2,:,:);
    L_freq       = data.L_freq(3,:,:);

    % --- Right bubble ---
    R_Md_vals    = data.R_freq(1,:,:);
    R_Delta_vals = data.R_freq(2,:,:);
    R_freq       = data.R_freq(3,:,:);

    ind = data.R_bubble.Md_to < max(data.R_bubble.Md_to);
    
    % ================= BURSTING =================
    ax1 = nexttile(i);
    scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq(:), 'filled'); hold on

    colormap(cmap); clim([0 8]);

    if i == n_methods
        cb = colorbar;
        cb.Label.String = 'Interburst frequency [Hz]';
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';
    end

    % xlabel('$M_{D1}$','Interpreter','latex')
    % ylabel('$\Delta_{v_{\theta}}$','Interpreter','latex')

    xlims = [0 0.8]; xlim(xlims); xticks(0:0.2:0.8)
    ylims = [0 6]; ylim(ylims); yticks(0:1:6)

    set(gca,'TickLabelInterpreter','latex','Layer','top','box','on',...
        'DataAspectRatio',[range(xlims),range(ylims),1])


    % ================= OSCILLATIONS =================
    ax2 = nexttile(i + 3);
    scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq(:), 'filled'); hold on

    colormap(cmap); clim([0 500]);

    if i == n_methods
        cb = colorbar;
        cb.Label.String = 'Oscillation frequency [Hz]';
        cb.Label.Interpreter = 'latex';
        cb.TickLabelInterpreter = 'latex';
    end

    % xlabel('$M_{D1}$','Interpreter','latex')
    % ylabel('$\Delta_{v_{\theta}}$','Interpreter','latex')

    xlims = [4 Md_max]; xlim(xlims); xticks(linspace(4,Md_max,5))
    ylims = [0 6]; ylim(ylims); yticks(0:1:6)

    set(gca,'TickLabelInterpreter','latex','Layer','top','box','on',...
        'DataAspectRatio',[range(xlims),range(ylims),1])


    % ================= BIFURCATION CURVES =================

    % ---- Hopf (left bubble) ----
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1a.mat'])
    Hx1 = x(end-2,:); Hy1 = x(end-1,:);
    GH = 118;

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
    Hx2 = x(end-2,:); Hy2 = x(end-1,:);

    plot(ax1,[flip(Hx1(1:GH)),Hx2],[flip(Hy1(1:GH)),Hy2],'k--','LineWidth',lwd);
    plot(ax1,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
    plot(ax1,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)

    % ---- Hopf (right bubble) ----
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
    Hx = x(end-2,:); Hy = x(end-1,:);
    GH = 1503;

    plot(ax2,Hx(1:GH),Hy(1:GH),'k--','LineWidth',lwd);
    plot(ax2,Hx(GH:end),Hy(GH:end),'k-','LineWidth',lwd);
    plot(ax2,Hx(GH),Hy(GH),'k.','MarkerSize',12)

    % ---- LPC (left) ----
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1a.mat'])
    LPCx = x(end-1,:); LPCy = x(end,:);

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1b.mat'])
    LPCx = [LPCx, x(end-1,:)]; LPCy = [LPCy, x(end,:)];

    plot(ax1,LPCx,LPCy,'k-.','LineWidth',lwd);

    % ---- LPC (right) ----
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2a.mat'])
    LPCx = x(end-1,:); LPCy = x(end,:);

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2b.mat'])
    LPCx = [flip(LPCx), x(end-1,:)]; LPCy = [flip(LPCy), x(end,:)];

    plot(ax2,LPCx,LPCy,'k-.','LineWidth',lwd);

end

% ===== Shared labels =====
xlabel(tlo,'$M_{D1}$','Interpreter','latex','FontSize',14)
ylabel(tlo,'$\Delta_{v_{\theta}}$','Interpreter','latex','FontSize',14)



%% NOT USED: COMPARISON OF FREQUENCY ESTIMATIONS
Md_max = 20;

methods = {'FP','LS'};
n_methods = numel(methods);

ms = 4;
cmap = bluewhitered(256);
lwd = 1.5;

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.25])

tlo = tiledlayout(2,n_methods+1,'TileSpacing','compact','Padding','tight');

% === Reference (EP) ===
load(['freq_data_EP_',num2str(Md_max),'.mat']); dataEP = data;

for i = 1:n_methods
    
    load(['freq_data_',methods{i},'_',num2str(Md_max),'.mat']);
    dataRef = data;

    % --- Left bubble ---
    L_Md_vals    = dataRef.L_freq(1,:,:);
    L_Delta_vals = dataRef.L_freq(2,:,:);
    L_freq_method = dataRef.L_freq(3,:,:);
    L_freq_EP     = dataEP.L_freq(3,:,:);
    L_freq_rel    = NaN(size(L_freq_EP));
    L_valid       = isfinite(L_freq_method) & isfinite(L_freq_EP) & (L_freq_EP > 0);
    L_freq_rel(L_valid) = 100 * (L_freq_method(L_valid)-L_freq_EP(L_valid)) ./ L_freq_EP(L_valid);

    % --- Right bubble ---
    R_Md_vals    = dataRef.R_freq(1,:,:);
    R_Delta_vals = dataRef.R_freq(2,:,:);
    R_freq_method = dataRef.R_freq(3,:,:);
    R_freq_EP     = dataEP.R_freq(3,:,:);
    R_freq_rel    = NaN(size(R_freq_EP));
    R_valid       = isfinite(R_freq_method) & isfinite(R_freq_EP) & (R_freq_EP > 0);
    R_freq_rel(R_valid) = 100 * (R_freq_method(R_valid)-R_freq_EP(R_valid)) ./ R_freq_EP(R_valid);

    % ================= BURSTING =================
    ax1 = nexttile(i);
    scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq_rel(:), 'filled'); hold on
    colormap(cmap); clim([-1 1]);

    xlims = [0 0.8]; xlim(xlims); xticks(0:0.2:0.8)
    ylims = [0 6]; ylim(ylims); yticks(0:1:6)

    set(gca,'TickLabelInterpreter','latex','Layer','top','box','on',...
        'DataAspectRatio',[range(xlims),range(ylims),1])


    % ================= OSCILLATIONS =================
    ax2 = nexttile(i + (n_methods+1));
    scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq_rel(:), 'filled'); hold on
    colormap(cmap); clim(30*[-1 1]);

    xlims = [4 Md_max]; xlim(xlims); xticks(linspace(4,Md_max,5))
    ylims = [0 6]; ylim(ylims); yticks(0:1:6)

    set(gca,'TickLabelInterpreter','latex','Layer','top','box','on',...
        'DataAspectRatio',[range(xlims),range(ylims),1])

    % ================= BIFURCATION CURVES =================
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1a.mat'])
    Hx1 = x(end-2,:); Hy1 = x(end-1,:); GH = 118;
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
    Hx2 = x(end-2,:); Hy2 = x(end-1,:);

    plot(ax1,[flip(Hx1(1:GH)),Hx2],[flip(Hy1(1:GH)),Hy2],'k--','LineWidth',lwd);
    plot(ax1,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
    plot(ax1,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
    Hx = x(end-2,:); Hy = x(end-1,:); GH = 1503;

    plot(ax2,Hx(1:GH),Hy(1:GH),'k--','LineWidth',lwd);
    plot(ax2,Hx(GH:end),Hy(GH:end),'k-','LineWidth',lwd);
    plot(ax2,Hx(GH),Hy(GH),'k.','MarkerSize',12)

end

% ================= 1D CUT (RIGHT COLUMN) =================
Delta_selected = 21; % approximately Delta = 0.5
ggreen = [0 0.8 0];
lwd2 = 0.75;

load(['freq_data_FP_',num2str(Md_max),'.mat']); dataFP = data;
load(['freq_data_LS_',num2str(Md_max),'.mat']); dataLS = data;
load(['freq_data_LC_',num2str(Md_max),'.mat']); dataLC = data;

% ---- BURSTING CUT ----
axc1 = nexttile(n_methods+1);
plot(dataFP.L_freq(1,:,Delta_selected), dataFP.L_freq(3,:,Delta_selected),'r','LineWidth',lwd2), hold on
plot(dataLS.L_freq(1,:,Delta_selected), dataLS.L_freq(3,:,Delta_selected),'--','Color',ggreen,'LineWidth',lwd2)
plot(dataEP.L_freq(1,:,Delta_selected), dataEP.L_freq(3,:,Delta_selected),'b','LineWidth',lwd2)
% plot(dataLC.L_freq(1,:,Delta_selected), dataLC.L_freq(3,:,Delta_selected),'r:','LineWidth',lwd2)

plot([0.148,0.643],[5.712,6.288],'k.','MarkerSize',12,'HandleVisibility','off')

xlims = [0 0.8]; ylims = [0 7];
xlim(xlims), ylim(ylims)

ylabel('$f$ [Hz]','Interpreter','latex','FontSize',14)

set(gca,'TickLabelInterpreter','latex','Layer','top','box','on', ...
    'DataAspectRatio',[range(xlims),range(ylims),1])

% ---- OSCILLATIONS CUT ----
axc2 = nexttile(2*(n_methods+1));
plot(dataFP.R_freq(1,:,Delta_selected), dataFP.R_freq(3,:,Delta_selected),'r','LineWidth',lwd2), hold on
plot(dataLS.R_freq(1,:,Delta_selected), dataLS.R_freq(3,:,Delta_selected),'--','Color',ggreen,'LineWidth',lwd2)
plot(dataEP.R_freq(1,:,Delta_selected), dataEP.R_freq(3,:,Delta_selected),'b','LineWidth',lwd2)
% plot(dataLC.R_freq(1,:,Delta_selected), dataLC.R_freq(3,:,Delta_selected),'r:','LineWidth',lwd2)

plot(4.582,323.8,'k.','MarkerSize',12,'HandleVisibility','off')

xlims = [4 Md_max]; ylims = [0 500];
xlim(xlims), ylim(ylims)

ylabel('$f$ [Hz]','Interpreter','latex','FontSize',14)

legend({'FP','LS','EP','LC'},'Location','southwest')
set(gca,'TickLabelInterpreter','latex','Layer','top','box','on', ...
    'DataAspectRatio',[range(xlims),range(ylims),1])

% ===== SHARED LABELS =====
xlabel(tlo,'$M_{D1}$','Interpreter','latex','FontSize',14)
ylabel(tlo,'$\Delta_{v_{\theta}}$','Interpreter','latex','FontSize',14)


%% FIG. 4: COMPARISON OF FREQUENCY ESTIMATIONS
Md_max = 20;

methods = {'FP','LS'};
n_methods = numel(methods);

ms = 4;
cmap = bluewhitered(256);
lwd = 1.5;

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.25])

tlo = tiledlayout(2,n_methods+1,'TileSpacing','compact','Padding','tight');

% === Reference (EP) ===
load(['freq_data_EP_',num2str(Md_max),'.mat']); dataEP = data;

for i = 1:n_methods
    
    load(['freq_data_',methods{i},'_',num2str(Md_max),'.mat']);
    dataRef = data;

    % --- Left bubble ---
    L_Md_vals    = dataRef.L_freq(1,:,:);
    L_Delta_vals = dataRef.L_freq(2,:,:);
    L_freq_method = dataRef.L_freq(3,:,:);
    L_freq_EP     = dataEP.L_freq(3,:,:);
    L_freq_rel    = NaN(size(L_freq_EP));
    L_valid       = isfinite(L_freq_method) & isfinite(L_freq_EP) & (L_freq_EP > 0);
    L_freq_rel(L_valid) = 100 * (L_freq_method(L_valid)-L_freq_EP(L_valid)) ./ L_freq_EP(L_valid);

    % --- Right bubble ---
    R_Md_vals    = dataRef.R_freq(1,:,:);
    R_Delta_vals = dataRef.R_freq(2,:,:);
    R_freq_method = dataRef.R_freq(3,:,:);
    R_freq_EP     = dataEP.R_freq(3,:,:);
    R_freq_rel    = NaN(size(R_freq_EP));
    R_valid       = isfinite(R_freq_method) & isfinite(R_freq_EP) & (R_freq_EP > 0);
    R_freq_rel(R_valid) = 100 * (R_freq_method(R_valid)-R_freq_EP(R_valid)) ./ R_freq_EP(R_valid);

    % ================= BURSTING =================
    ax1 = nexttile(i);
    scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq_rel(:), 'filled'); hold on
    colormap(cmap); clim(30*[-1 1]); 
    
    xlims = [0 0.8]; xlim(xlims); xticks(0:0.2:0.8)
    ylims = [0 6]; ylim(ylims); yticks(0:1:6)

    set(gca,'TickLabelInterpreter','latex','Layer','top','box','on',...
        'DataAspectRatio',[range(xlims),range(ylims),1])

    % ================= OSCILLATIONS =================
    ax2 = nexttile(i + (n_methods+1));
    scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq_rel(:), 'filled'); hold on
    colormap(cmap); clim(30*[-1 1]);

    xlims = [4 Md_max]; xlim(xlims); xticks(linspace(4,Md_max,5))
    ylims = [0 6]; ylim(ylims); yticks(0:1:6)

    set(gca,'TickLabelInterpreter','latex','Layer','top','box','on',...
        'DataAspectRatio',[range(xlims),range(ylims),1])

    % ===== COLORBARS (jen v pravém sloupci daného řádku) =====
    if i == n_methods
        cb1 = colorbar(ax1,'eastoutside');
        cb1.Label.String = 'Relative difference [\%]';
        cb1.Label.Interpreter = 'latex';
        cb1.TickLabelInterpreter = 'latex';
        cb1.Ticks = -30:10:30;
        cb1.TickLabels = {'$<-30$','$-20$','$-10$','$0$','$10$','$20$','$>30$'}; 

        cb2 = colorbar(ax2,'eastoutside');
        cb2.Label.String = 'Relative difference [\%]';
        cb2.Label.Interpreter = 'latex';
        cb2.TickLabelInterpreter = 'latex';
        cb2.Ticks = -30:10:30;
        cb2.TickLabels = {'$<-30$','$-20$','$-10$','$0$','$10$','$20$','$>30$'}; 
    end

    % ================= BIFURCATION CURVES =================
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1a.mat'])
    Hx1 = x(end-2,:); Hy1 = x(end-1,:); GH = 118;
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
    Hx2 = x(end-2,:); Hy2 = x(end-1,:);

    plot(ax1,[flip(Hx1(1:GH)),Hx2],[flip(Hy1(1:GH)),Hy2],'k--','LineWidth',lwd);
    plot(ax1,Hx1(GH:end),Hy1(GH:end),'k-','LineWidth',lwd);
    plot(ax1,Hx1(GH),Hy1(GH),'k.','MarkerSize',12)

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\H_DeltaMd_1b.mat'])
    Hx = x(end-2,:); Hy = x(end-1,:); GH = 1503;

    plot(ax2,Hx(1:GH),Hy(1:GH),'k--','LineWidth',lwd);
    plot(ax2,Hx(GH:end),Hy(GH:end),'k-','LineWidth',lwd);
    plot(ax2,Hx(GH),Hy(GH),'k.','MarkerSize',12)

    % ---- LPC (left) ----
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1a.mat'])
    LPCx = x(end-1,:); LPCy = x(end,:);

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_1b.mat'])
    LPCx = [LPCx, x(end-1,:)]; LPCy = [LPCy, x(end,:)];

    plot(ax1,LPCx,LPCy,'k-.','LineWidth',lwd);

    % ---- LPC (right) ----
    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2a.mat'])
    LPCx = x(end-1,:); LPCy = x(end,:);

    load([Path,'\Dopamine_NMM_fast4\diagram_DeltaMd\LPC_DeltaMd_2b.mat'])
    LPCx = [flip(LPCx), x(end-1,:)]; LPCy = [flip(LPCy), x(end,:)];

    plot(ax2,LPCx,LPCy,'k-.','LineWidth',lwd);

end

% ================= 1D CUT (RIGHT COLUMN) =================
Delta_selected = 21; % approximately Delta = 0.5
ggreen = [0 0.8 0];
lwd2 = 0.75;

load(['freq_data_FP_',num2str(Md_max),'.mat']); dataFP = data;
load(['freq_data_LS_',num2str(Md_max),'.mat']); dataLS = data;
load(['freq_data_LC_',num2str(Md_max),'.mat']); dataLC = data;

% ---- BURSTING CUT ----
axc1 = nexttile(n_methods+1);
plot(dataFP.L_freq(1,:,Delta_selected), dataFP.L_freq(3,:,Delta_selected),'r','LineWidth',lwd2), hold on
plot(dataLS.L_freq(1,:,Delta_selected), dataLS.L_freq(3,:,Delta_selected),'--','Color',ggreen,'LineWidth',lwd2)
plot(dataEP.L_freq(1,:,Delta_selected), dataEP.L_freq(3,:,Delta_selected),'b','LineWidth',lwd2)

plot([0.148,0.643],[5.712,6.288],'k.','MarkerSize',12,'HandleVisibility','off')

xlims = [0 0.8]; ylims = [0 7];
xlim(xlims), xticks(0:0.2:0.8)
ylim(ylims), yticks(0:1:7)

ylabel('$f$ [Hz]','Interpreter','latex')

set(gca,'TickLabelInterpreter','latex','Layer','top','box','on', ...
    'DataAspectRatio',[range(xlims),range(ylims),1])

% ---- OSCILLATIONS CUT ----
axc2 = nexttile(2*(n_methods+1));
plot(dataFP.R_freq(1,:,Delta_selected), dataFP.R_freq(3,:,Delta_selected),'r','LineWidth',lwd2), hold on
plot(dataLS.R_freq(1,:,Delta_selected), dataLS.R_freq(3,:,Delta_selected),'--','Color',ggreen,'LineWidth',lwd2)
plot(dataEP.R_freq(1,:,Delta_selected), dataEP.R_freq(3,:,Delta_selected),'b','LineWidth',lwd2)

plot(4.582,323.8,'k.','MarkerSize',12,'HandleVisibility','off')

xlims = [4 Md_max]; ylims = [0 500];
xlim(xlims), xticks(4:4:Md_max)
ylim(ylims), yticks(0:100:500)

ylabel('$f$ [Hz]','Interpreter','latex')

legend({'FP','LS','EP','LC'},'Location','southwest','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex','Layer','top','box','on', ...
    'DataAspectRatio',[range(xlims),range(ylims),1])

% ===== SHARED LABELS =====
xlabel(tlo,'$M_{D1}$','Interpreter','latex','FontSize',14)
ylabel(tlo,'$\Delta_{v_{\theta}}$','Interpreter','latex','FontSize',14)




