clear
close all
warning('off')

% MatCont procedures initialization
init

% Parallel pool initialization
if isempty(gcp('nocreate'))
    parpool;
end

% LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');


% %% PRECOMPUTED CURVES DELINEATING THE BUBBLES
% % Load the Hopf curve
% PATH = './Systems/Dopamine_NMM_fast4/diagram_DeltaMd/';
% curve_Hopf1 = load(strcat(PATH,'H_DeltaMd_1a.mat'));
% curve_Hopf2 = load(strcat(PATH,'H_DeltaMd_1b.mat'));
% 
% EP = [curve_Hopf1.x(1:end-3,end:-1:2), curve_Hopf2.x(1:end-3,:)];
% Md = [curve_Hopf1.x(end-2,end:-1:2), curve_Hopf2.x(end-2,:)];
% Delta = [curve_Hopf1.x(end-1,end:-1:2), curve_Hopf2.x(end-1,:)];


%% HOPF CURVE CONTINUATION
% (1) convergence to a stable equilibrium point (EP)
% (2) continuation of the EP toward the Hopf bifurcation (HB)
% (3) continuation of the HB curves in both directions
% (4) extraction of HB branches delineating the left (bursting) and right
%     (spiking) bubble

tstamp = datetime;

% System of ODEs
hls = Dopamine_NMM_fast4;

% Free parameters
Md_idx    = 4;
Delta_idx = 9;
par_Md    = 4;
par_Delta = 5;

% Modifiable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default (fixed) parameters
par_C      = 50;
par_Iext   = 0;
par_ujump  = 45;
par_B      = 0.2;
par_gammaa = 1.1;
par_Vr     = -60;
par_Vt     = -43;   % Vt > Vr
par_cexc   = 0;
par_Ta     = 6;
par_K      = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pars = [par_C, par_Iext, par_ujump, par_Md, par_B, par_gammaa, par_Vr, ...
    par_Vt, par_Delta, par_cexc, par_Ta, par_K];

odefun = @(t,y,par_Md,par_Delta) hls{2}(t, y, pars(1), pars(2), ...
    pars(3), par_Md, pars(5), pars(6), pars(7), pars(8), par_Delta, ...
    pars(10), pars(11), pars(12));

% Options setup
optODE = odeset('RelTol',1e-5,'MaxStep',0.05);

% Convergence to the (stable) EP
IC    = [0.006 -70 -50 0.06]; % r, V, u, Sa
[t,y] = ode45(@(t,y) odefun(t,y,par_Md,par_Delta),[0 1000], IC, optODE);
EP    = y(end,:)';

% figure, plot(t,y(:,1))

% Direction of continuation
apEP    = Delta_idx;   % EP continuation parameter: Delta
[x0,v0] = init_EP_EP(@Dopamine_NMM_fast4,EP,pars,apEP);

optEP = contset;
optEP = contset(optEP,'MaxNumPoints',20);
optEP = contset(optEP,'InitStepsize',0.01);
optEP = contset(optEP,'MinStepsize',0.001);
optEP = contset(optEP,'MaxStepsize',0.1);
% optEP = contset(optEP,'VarTolerance',1e-06);
% optEP = contset(optEP,'FunTolerance',1e-06);
% optEP = contset(optEP,'TestTolerance',1e-06);
optEP = contset(optEP,'Adapt',3);
optEP = contset(optEP,'Singularities',1);
optEP = contset(optEP,'Eigenvalues',1);
[x,~,~,~,~] = cont(@equilibrium,x0,[],optEP);

if x(end,1) <= x(end,end)
    optEP = contset(optEP,'Backward',1);  % wrong direction => backward continuation
end

% EP continuation
optEP = contset(optEP,'MaxNumPoints',2000);
[x,v,s,h,f] = cont(@equilibrium,x0,[],optEP);

% Select Hopf bifurcation point
Hopf_idx = find(strcmp({s.label},'H'),1);
if isempty(Hopf_idx), Hopf_idx = find(strcmp({s.label},'H '),1); end % sometimes 'H ' instead of 'H'

if isempty(Hopf_idx)
    error('No Hopf bifurcation detected!');
else
    idx = s(Hopf_idx).index;
end

% Initialization of Hopf curve continuation
pars(apEP) = x(end,idx);    % update parameter value found at Hopf point
apH = [Md_idx Delta_idx];   % Hopf continuation parameters
[xH,vH] = init_H_H(@Dopamine_NMM_fast4,x(1:end-1,idx),pars',apH);

% Forward continuation
optH = contset;
optH = contset(optH,'MaxNumPoints',10000);
optH = contset(optH,'InitStepsize',0.01);
optH = contset(optH,'MinStepsize',0.001);
optH = contset(optH,'MaxStepsize',0.1);
optH = contset(optH,'Singularities',0);
optH = contset(optH,'Backward',0);
[xHf,vHf,sHf,hHf,fHf] = cont(@hopf,xH,vH,optH);

% Backward continuation
optH = contset(optH,'Backward',1);
[xHb,vHb,sHb,hHb,fHb] = cont(@hopf,xH,vH,optH);

% Hopf curve
EP    = [xHb(1:end-3,end:-1:2), xHf(1:end-3,:)];
Md    = [xHb(end-2,end:-1:2), xHf(end-2,:)];
Delta = [xHb(end-1,end:-1:2), xHf(end-1,:)];

% figure, plot(Md, Delta,'.')

% Selection of positive Md and Delta
ind_positive = Md > 0 & Delta > 0;
EP    = EP(:,ind_positive);
Md    = Md(ind_positive); 
Delta = Delta(ind_positive);
 
% figure, plot(Md, Delta,'.')

% Modifiable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Md threshold separating bursting and spiking bubbles
thrs = 2;
Md_max = 20;

% Grid size
L_grid_size = [300,200];
R_grid_size = [400,200];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Left bubble for Md (bursting)
L_range = [0 thrs];
L_bubble = extract_branch(EP, Md, Delta, L_range, L_grid_size, 'L');

% Right bubble for Md (spiking)
R_range = [thrs Md_max];
R_bubble = extract_branch(EP, Md, Delta, R_range, R_grid_size, 'R');

% Visualization of these curves
fig = figure; 
plot(L_bubble.Md_from, L_bubble.Delta_from,'.', ...
    L_bubble.Md_to, L_bubble.Delta_to,'.', ...
    R_bubble.Md_from, R_bubble.Delta_from,'.', ...
    R_bubble.Md_to, R_bubble.Delta_to,'.')
xlim([min([L_range, R_range]), max([L_range, R_range])])
title('Are the bubbles separated correctly? If so, press any key.')

% pause
close(fig)


%% SETTING

% Modifiable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method of frequency estimation
% opts.method = 'LS';     % SIM & Lomb–Scargle periodogram [~560 s, ~1500 s]
% opts.method = 'FP';     % SIM & findpeaks [~260 s]
opts.method = 'EP';     % Equilibrium Point continuation [~40 s]
% opts.method = 'LC';     % Limit Cycle continuation [~30 min]


% ---- SIM options ----
% Simulation time span
opts.SIM.T_span.L = [0 2000];
opts.SIM.T_span.R = [0 100];

% Maximum frequency [Hz]
opts.SIM.f_max.L = 50;
opts.SIM.f_max.R = 500;

% RHS and solver setting
opts.SIM.odefun = odefun;
opts.SIM.optODE = odeset('RelTol',1e-3);


% ---- EP options ----
% Maximum number of continuation steps
opts.EP.MaxNumPts.L = 2000;
opts.EP.MaxNumPts.R = 2000;

% Maximum continuation step size
opts.EP.MaxStepsize.L = 0.1;
opts.EP.MaxStepsize.R = 0.75;

% Other parameters
opts.EP.Md_idx    = Md_idx;
opts.EP.Delta_idx = Delta_idx;
opts.EP.pars      = pars;


% ---- LC options ----
% Maximum number of continuation steps
opts.LC.MaxNumPts.L = 300;
opts.LC.MaxNumPts.R = 300;

% Maximum continuation step size
opts.LC.MaxStepsize.L = 10;
opts.LC.MaxStepsize.R = 40;

% Other parameters
opts.LC.Md_idx    = Md_idx;
opts.LC.Delta_idx = Delta_idx;
opts.LC.pars      = pars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% COMPUTATION

t_method = tic;

[L_freq, L_cost] = compute_frequency_grid(L_bubble,opts);  % left bubble
[R_freq, R_cost] = compute_frequency_grid(R_bubble,opts);  % right bubble

method_wall_time = toc(t_method);

cost = combine_costs(L_cost, R_cost);
cost.method_wall_time = method_wall_time;
cost.method = opts.method;

print_cost_report(cost);

% Saving
data = struct('method', opts.method, 'pars', pars, 'L_bubble', L_bubble, ...
    'R_bubble', R_bubble, 'L_freq', L_freq, 'R_freq', R_freq, ...
    'opts', opts, 'cost', cost);

save(['freq_data_',opts.method,'_new.mat'],'data')

warning('on')


%% PLOTTING ===============================================================
% Data loading ------------------------------------------------------------
Md_max = 10;
Md_max = 20;

% load(['freq_data_FP_',num2str(Md_max),'.mat']); 
% load(['freq_data_LS_',num2str(Md_max),'.mat']);
% load(['freq_data_EP_',num2str(Md_max),'.mat']); 
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
ms = 4;             % marker size
cmap = parula(256); %[[1, 1, 1]; [1, 1, 1]; parula(254)];

fig = figure;
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[13,6.5])
tlo = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');


% BURSTING ----------------------------------------------------------------
nexttile([1,1])
scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq(:), 'filled'), hold on
plot(L_Md,L_Delta,'k-','LineWidth',1)

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
xlim(Md_lims); xticks(Md_ticks)
ylim(Delta_lims); yticks(Delta_ticks)

xaxis = get(gca, 'XAxis'); xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis'); yaxis.TickLabelInterpreter = 'latex';
set(gca, 'Layer', 'top', 'box', 'on')


% SPIKING -----------------------------------------------------------------
nexttile([1,1])
scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq(:), 'filled'), hold on
plot(R_Md,R_Delta,'k-','LineWidth',1)

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
xlim(Md_lims); xticks(Md_ticks);
ylim(Delta_lims); yticks(Delta_ticks)

xaxis = get(gca, 'XAxis'); xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis'); yaxis.TickLabelInterpreter = 'latex';
set(gca, 'Layer', 'top', 'box', 'on')


% BURSTING & SPIKING ------------------------------------------------------

nexttile([1,2])
scatter(L_Md_vals(:), L_Delta_vals(:), ms, L_freq(:), 'filled'), hold on
plot(L_Md,L_Delta,'k-','LineWidth',1)

scatter(R_Md_vals(:), R_Delta_vals(:), ms, R_freq(:), 'filled')
plot(R_Md,R_Delta,'k-','LineWidth',1)

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
xlim(Md_lims); xticks(Md_ticks)
ylim(Delta_lims); yticks(Delta_ticks)

xaxis = get(gca, 'XAxis'); xaxis.TickLabelInterpreter = 'latex';
yaxis = get(gca, 'YAxis'); yaxis.TickLabelInterpreter = 'latex';
set(gca, 'Layer', 'top', 'box', 'on')

savefig(fig, 'freq.fig');

if exist('tstamp','var')
    fprintf('\nTOTAL ELAPSED TIME = %.1f secs\n', seconds(datetime-tstamp))
end


%% INTERBURST (AND SPIKING) FREQUENCY ESTIMATION

function [freq_matrix, cost] = compute_frequency_grid(bubble, opts)
% COMPUTE_FREQUENCY_GRID
% ------------------------------------------------------------
% Computes grid-based frequency estimation over a 2D parameter space
% (Md, Delta) inside a given bubble (bursting or spiking).
%
% This function provides a common interface for multiple frequency
% estimation methods:
%   - LS: Lomb–Scargle periodogram (frequency-domain)
%   - FP: peak detection (time-domain)
%   - EP: equilibrium continuation + eigenvalues
%   - LC: limit cycle continuation
%
% In addition, this modified version collects basic computational-cost
% characteristics in the output structure cost.

b          = bubble.spec;
EP_from    = bubble.EP_from;
Md_from    = bubble.Md_from;
Md_to      = bubble.Md_to;
Md_grid    = bubble.Md_grid;
Delta_from = bubble.Delta_from;
method     = opts.method;

N_Delta = length(Delta_from);
N_Md    = length(Md_grid);

freq_matrix = NaN(4, N_Md, N_Delta);

cost_sections(N_Delta) = empty_cost_section();
t_grid = tic;

% --- select method options ---
switch method
    case {'LS','FP'}
        m = opts.SIM;
    case {'EP','LC'}
        m = opts.(method);
    otherwise
        error('Unknown method!')
end

% --- computation ---
parfor i = 1:N_Delta

    fprintf('Computing iteration %s: %d/%d\n', b, i, N_Delta);

    section_cost = empty_cost_section();
    section_cost.method = method;
    section_cost.bubble = b;
    section_cost.slice_index = i;

    t_section = tic;

    IC = EP_from(:,i);

    % --- local buffer required for parfor ---
    freq_section = NaN(4, N_Md);

    % --- select Md points inside bubble ---
    idx_inside = Md_grid > Md_from(i) & Md_grid < Md_to(i);
    Md_inside  = Md_grid(idx_inside);

    section_cost.n_total = N_Md;
    section_cost.n_inside = numel(Md_inside);

    if isempty(Md_inside)
        freq_matrix(:,:,i) = freq_section;

        section_cost.wall_s = toc(t_section);
        cost_sections(i) = section_cost;

        continue
    end

    % ============================================================
    % EP METHOD
    if strcmp(method,'EP')

        par = m.pars;
        par(m.Delta_idx) = Delta_from(i);

        [f_dom, f_sec, section_cost_method] = compute_frequency_section_EP( ...
            IC, par, Md_from(i), Md_to(i), Md_inside, ...
            m.Md_idx, m.MaxNumPts.(b), m.MaxStepsize.(b));

    % ============================================================
    % LC METHOD
    elseif strcmp(method,'LC')

        par = m.pars;
        par(m.Delta_idx) = Delta_from(i);
        
        if Md_from(i) < 2 % left
            [f_dom, f_sec, section_cost_method] = compute_frequency_section_LC( ...
                IC, par, Md_from(i), Md_to(i), Md_inside, ...
                m.Md_idx, m.MaxNumPts.(b), m.MaxStepsize.(b));
        else % right
            [f_dom, f_sec, section_cost_method] = compute_frequency_section_LC_bwd( ...
                par, Md_from(i), Md_to(i), Md_inside, ...
                m.Md_idx, m.MaxNumPts.(b), m.MaxStepsize.(b));
        end

    % ============================================================
    % SIMULATION METHODS (LS / FP)
    else

        [f_dom, f_sec, section_cost_method] = compute_frequency_section_SIM( ...
            @(t,y,par_Md) m.odefun(t,y,par_Md,Delta_from(i)), ...
            m.T_span.(b), m.optODE, Md_inside, IC, m.f_max.(b), method);

    end

    % --- fill slice ---
    idx = find(idx_inside);

    freq_section(1,idx) = Md_inside;
    freq_section(2,idx) = Delta_from(i);
    freq_section(3,idx) = f_dom;
    freq_section(4,idx) = f_sec;

    section_cost.wall_s = toc(t_section);
    section_cost = add_cost_sections(section_cost, section_cost_method);

    freq_matrix(:,:,i) = freq_section;
    cost_sections(i) = section_cost;

end

cost = summarize_cost_sections(cost_sections);
cost.grid_wall_s = toc(t_grid);
cost.method = method;
cost.bubble = b;

end


%% EP CONTINUATION

function [f_dom, f_sec, cost] = compute_frequency_section_EP( ...
    EP, pars, Md_from, Md_to, Md_grid, Md_idx, MaxNumPoints, MaxStepsize)

% COMPUTE_FREQUENCY_SECTION_EP
% ------------------------------------------------------------
% Frequency estimation via Equilibrium Point (EP) continuation and eigenvalues.
%
% IDEA:
%   Continue equilibrium from Hopf point -> extract eigenvalues ->
%   frequency ~ Im(lambda)/(2*pi)

f_dom = NaN(size(Md_grid));
f_sec = NaN(size(Md_grid));

cost = empty_cost_section();
cost.method = 'EP';
cost.n_requested = numel(Md_grid);

par         = pars;
par(Md_idx) = Md_from;

ap = Md_idx;

[x0,~] = init_EP_EP(@Dopamine_NMM_fast4,EP,par,ap);

optEP = contset;
optEP = contset(optEP,'MaxNumPoints',20);
optEP = contset(optEP,'InitStepsize',0.1);
optEP = contset(optEP,'MinStepsize',0.01);
optEP = contset(optEP,'MaxStepsize',MaxStepsize);
optEP = contset(optEP,'Singularities',0);
optEP = contset(optEP,'Eigenvalues',1);

% --- choose direction ---
try
    [x,~,~,~,~] = cont(@equilibrium,x0,[],optEP);
    cost.n_cont_runs = cost.n_cont_runs + 1;
    cost.n_cont_points = cost.n_cont_points + size(x,2);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

if size(x,2) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

Backward = x(end,1) >= x(end,end); % continuation from left to right

% --- full continuation ---
optEP = contset(optEP,'MaxNumPoints',MaxNumPoints);
optEP = contset(optEP,'Backward',Backward);

try
    [x,~,~,~,f] = cont(@equilibrium,x0,[],optEP);
    cost.n_cont_runs = cost.n_cont_runs + 1;
    cost.n_cont_points = cost.n_cont_points + size(x,2);
    cost.n_jac_eig = cost.n_jac_eig + size(f,2);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

if size(x,2) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

idx_inside = x(end,:) < Md_to;
if sum(idx_inside) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

Md_cont = x(end,idx_inside);

[f_dom_cont, f_sec_cont] = dominant_frequencies_EP(f(:,idx_inside));

% --- interpolate to fixed grid ---
f_dom = interp1(Md_cont, f_dom_cont, Md_grid, 'linear');
f_sec = interp1(Md_cont, f_sec_cont, Md_grid, 'linear');

cost.n_valid_dom = sum(~isnan(f_dom));
cost.n_valid_sec = sum(~isnan(f_sec));

end 


%% LC CONTINUATION

function [f_dom, f_sec, cost] = compute_frequency_section_LC( ...
    EP, pars, Md_from, Md_to, Md_grid, Md_idx, MaxNumPoints, MaxStepsize)

% COMPUTE_FREQUENCY_SECTION_LC
% ------------------------------------------------------------
% Frequency estimation via Limit Cycle (LC) continuation.
%
% IDEA:
%   Continue limit cycles from Hopf -> extract period T ->
%   frequency = 1/T

f_dom = NaN(size(Md_grid));
f_sec = NaN(size(Md_grid));

cost = empty_cost_section();
cost.method = 'LC';
cost.n_requested = numel(Md_grid);

par         = pars;
par(Md_idx) = Md_from;

try
    [x0,v0] = init_H_LC(@Dopamine_NMM_fast4,EP,par,Md_idx,1e-6,150,4);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

InitStepsize = MaxStepsize;

for i = 1:5

    optLC = contset;
    optLC = contset(optLC,'MaxNumPoints',MaxNumPoints);
    optLC = contset(optLC,'InitStepsize',InitStepsize);
    optLC = contset(optLC,'MinStepsize',0.01);
    optLC = contset(optLC,'MaxStepsize',MaxStepsize);
    % optLC = contset(optLC,'VarTolerance',1e-05);
    % optLC = contset(optLC,'FunTolerance',1e-05);
    % optLC = contset(optLC,'TestTolerance',1e-04);
    optLC = contset(optLC,'Multipliers',1);
    optLC = contset(optLC,'Singularities',0);
    
    try
        [x,~,~,~,~] = cont(@limitcycle,x0,v0,optLC);
        cost.n_cont_runs = cost.n_cont_runs + 1;
        cost.n_cont_points = cost.n_cont_points + size(x,2);
        cost.n_lc_points = cost.n_lc_points + size(x,2);
    catch
        cost.n_fail = cost.n_fail + 1;
        return
    end
    
    if size(x,2) > 2 
        break % OK
    end
    
    InitStepsize = 2 * InitStepsize;
end

if size(x,2) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

Md_cont = x(end,:);
[~, idx_LPC] = max(Md_cont);

idx_before_LPC = 1:length(Md_cont) <= idx_LPC;
idx_inside     = (Md_cont < Md_to) & idx_before_LPC;

if ~any(idx_inside)
    cost.n_fail = cost.n_fail + 1;
    return
end

Md_cont = x(end,idx_inside);
T       = x(end-1,idx_inside);

f_cont = 1000 ./ T;     % ms -> s

% --- interpolate to fixed grid ---
if length(Md_cont) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

f_dom = interp1(Md_cont, f_cont, Md_grid, 'linear', NaN);

cost.n_valid_dom = sum(~isnan(f_dom));
cost.n_valid_sec = sum(~isnan(f_sec));

end


function [f_dom, f_sec, cost] = compute_frequency_section_LC_bwd(pars, Md_from, ...
    Md_to, Md_grid, Md_idx, MaxNumPoints, MaxStepsize)

% COMPUTE_FREQUENCY_SECTION_LC_BWD
% ------------------------------------------------------------
% Frequency estimation via backward Limit Cycle (LC) continuation.

f_dom = NaN(size(Md_grid));
f_sec = NaN(size(Md_grid));

cost = empty_cost_section();
cost.method = 'LC';
cost.n_requested = numel(Md_grid);

par         = pars;
par(Md_idx) = Md_to;

hls = Dopamine_NMM_fast4;

odefun = @(t,y) hls{2}(t, y, par(1), par(2), par(3), par(4), par(5), ...
    par(6), par(7), par(8), par(9), par(10), par(11), par(12));

% Options setup
optODE = odeset('AbsTol',1e-7,'RelTol',1e-7,'MaxStep',0.01);

% Convergence to the (stable) LC
IC     = [ 0.016 3.974 521.684 24.083 ]; % r, V, u, Sa

try
    cost.n_ode_solves = cost.n_ode_solves + 1;
    [t,y]  = ode45(@(t,y) odefun(t,y),[0 1000], IC, optODE);
    cost.n_time_points = cost.n_time_points + length(t);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

if isempty(y)
    cost.n_fail = cost.n_fail + 1;
    return
end

new_IC = y(end,:)';

% Initialization
T_guess = 6;

try
    cost.n_ode_solves = cost.n_ode_solves + 1;
    [t,y] = ode45(@(t,y) odefun(t,y), [0 T_guess], new_IC, optODE);
    cost.n_time_points = cost.n_time_points + length(t);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

try
    t = t-t(1);
    [x0,v0] = initOrbLC(@Dopamine_NMM_fast4,t,y,par,Md_idx,250,4,1e-4);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

optLC = contset;
optLC = contset(optLC,'MaxNumPoints',10);
optLC = contset(optLC,'InitStepsize',0.01);
optLC = contset(optLC,'MinStepsize',0.001);
optLC = contset(optLC,'MaxStepsize',MaxStepsize);
optLC = contset(optLC,'VarTolerance',1e-03);
optLC = contset(optLC,'FunTolerance',1e-03);
optLC = contset(optLC,'TestTolerance',1e-03);
optLC = contset(optLC,'Multipliers',0);
optLC = contset(optLC,'Singularities',0);

% Probe continuation
try
    [x,~,~,~,~] = cont(@limitcycle,x0,v0,optLC);
    cost.n_cont_runs = cost.n_cont_runs + 1;
    cost.n_cont_points = cost.n_cont_points + size(x,2);
    cost.n_lc_points = cost.n_lc_points + size(x,2);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

if size(x,2) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

% Full continuation
optLC = contset(optLC,'MaxNumPoints',MaxNumPoints);
if x(end,1) <= x(end,end) % we want to continuate from right to left
    optLC = contset(optLC,'Backward',1);  % wrong direction => backward continuation
end

try
    [x,~,~,~,~] = cont(@limitcycle,x0,v0,optLC);
    cost.n_cont_runs = cost.n_cont_runs + 1;
    cost.n_cont_points = cost.n_cont_points + size(x,2);
    cost.n_lc_points = cost.n_lc_points + size(x,2);
catch
    cost.n_fail = cost.n_fail + 1;
    return
end

if size(x,2) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

Md_cont = x(end,:);
[~, idx_LPC] = min(Md_cont);

idx_before_LPC = 1:length(Md_cont) <= idx_LPC;
idx_inside     = (Md_cont > Md_from) & idx_before_LPC; % 'from' is on the left

if ~any(idx_inside)
    cost.n_fail = cost.n_fail + 1;
    return
end

Md_cont = x(end,idx_inside);
T       = x(end-1,idx_inside);

f_cont = 1000 ./ T;     % ms -> s

% --- interpolate to fixed grid ---
if length(Md_cont) < 2
    cost.n_fail = cost.n_fail + 1;
    return
end

f_dom = interp1(Md_cont, f_cont, Md_grid, 'linear', NaN);

cost.n_valid_dom = sum(~isnan(f_dom));
cost.n_valid_sec = sum(~isnan(f_sec));

end


%% SIMULATION METHODS

function [f_dom, f_sec, cost] = compute_frequency_section_SIM(odefun, T_span, optODE, Md_vals, EP, f_max, method)

% COMPUTE_FREQUENCY_SECTION_SIM
% ------------------------------------------------------------
% Frequency estimation via direct orbit simulation.
%
% METHODS:
%   LS        - Lomb–Scargle periodogram
%   FP        - findpeaks
%
% This modified version also counts:
%   - number of ODE solves
%   - number of RHS evaluations
%   - number of stored time points
%   - number of valid frequency estimates

N = length(Md_vals);

f_dom = NaN(1,N);
f_sec = NaN(1,N);

cost = empty_cost_section();
cost.method = method;
cost.n_requested = N;

rhs_count = 0;
ode_solve_count = 0;
n_time_points = 0;

% Wrapper for counting RHS evaluations
    function dy = counted_odefun(t,y,par_Md)
        rhs_count = rhs_count + 1;
        dy = odefun(t,y,par_Md);
    end

% --- transient ---
try
    ode_solve_count = ode_solve_count + 1;
    [~,y] = ode45(@(t,y) counted_odefun(t,y,Md_vals(1)), 10*T_span, EP, optODE);
catch
    cost.n_fail = cost.n_fail + 1;
    cost.n_rhs = rhs_count;
    cost.n_ode_solves = ode_solve_count;
    return
end

if isempty(y)
    cost.n_fail = cost.n_fail + 1;
    cost.n_rhs = rhs_count;
    cost.n_ode_solves = ode_solve_count;
    return
end

IC = y(end,:)';

frequencies = linspace(0, f_max, 2000); % 500
minPeakDistance_ms = 1000 / f_max;

for i = 1:N
    
    % --- simulation ---
    try
        ode_solve_count = ode_solve_count + 1;
        [t,y] = ode45(@(t,y) counted_odefun(t,y,Md_vals(i)), T_span, IC, optODE);
    catch
        cost.n_fail = cost.n_fail + 1;
        continue
    end

    if length(t) < 2
        cost.n_fail = cost.n_fail + 1;
        continue
    end

    n_time_points = n_time_points + length(t);

    IC = y(end,:)'; % IC for the next iteration

    signal = y(:,1);
    
    % --- frequency computation ---
    if strcmp(method,'LS')
        [f_dom(i), f_sec(i)] = dominant_frequencies_LS(signal, t, frequencies);
    else
        [f_dom(i), f_sec(i)] = dominant_frequencies_FP(signal, t, minPeakDistance_ms);
    end
end

cost.n_rhs = rhs_count;
cost.n_ode_solves = ode_solve_count;
cost.n_time_points = n_time_points;
cost.n_valid_dom = sum(~isnan(f_dom));
cost.n_valid_sec = sum(~isnan(f_sec));

end


%% EIGENVALUES -> FREQUENCY

function [f_dom, f_sec] = dominant_frequencies_EP(eigvals)

% DOMINANT_FREQUENCIES_EP
% ------------------------------------------------------------
% Convert eigenvalues -> frequencies using eigenvalue branch tracking.
%   This function identifies the relevant complex eigenvalue pair at 
%   the first continuation point and tracks the same eigenvalue branch 
%   along the continuation.
%
% INPUT:
%   eigvals   - matrix of eigenvalues, each column = spectrum
%
% OUTPUT:
%   f_dom     - tracked dominant frequency [Hz]
%   f_sec     - tracked secondary frequency [Hz], if initialized

n_cols = size(eigvals, 2);

f_dom = nan(1, n_cols);
f_sec = nan(1, n_cols);

% Tolerances
imag_tol = 1e-10;   % eigenvalue is considered complex if |Im(lambda)| > imag_tol
freq_tol = 1e-8;    % used only for identifying distinct modes at initialization

if n_cols == 0
    return
end

% --- initialization from the first spectrum ---
lam0 = eigvals(:,1);

% Prefer eigenvalues with positive imaginary part to avoid conjugate duplicates.
idx_complex = find(imag(lam0) > imag_tol);

% Fallback: if for some reason only negative-imaginary representatives are present
if isempty(idx_complex)
    idx_complex = find(abs(imag(lam0)) > imag_tol);
end

if isempty(idx_complex)
    % No complex mode at the starting point -> no frequency estimate
    return
end

% Sort candidate complex eigenvalues by closeness of real part to imaginary axis
[~, order] = sort(abs(real(lam0(idx_complex))));
idx_complex = idx_complex(order);

% Build a list of distinct initial modes.
% This removes possible conjugate duplicates if both signs of Im(lambda) are present.
mode_idx = [];
mode_freqs = [];

for k = 1:length(idx_complex)

    idx_k = idx_complex(k);
    omega_k = abs(imag(lam0(idx_k)));
    f_k = omega_k / (2*pi) * 1000;  % ms^{-1} -> Hz

    if isempty(mode_freqs)
        mode_idx(end+1) = idx_k; %#ok<AGROW>
        mode_freqs(end+1) = f_k; %#ok<AGROW>
    else
        rel_diff = abs(f_k - mode_freqs) ./ max(abs(mode_freqs), eps);

        if all(rel_diff > freq_tol)
            mode_idx(end+1) = idx_k; %#ok<AGROW>
            mode_freqs(end+1) = f_k; %#ok<AGROW>
        end
    end
end

if isempty(mode_idx)
    return
end

% Dominant mode = first complex mode closest to the imaginary axis
lambda_dom_ref = lam0(mode_idx(1));
has_dom = true;

% Secondary mode = second distinct complex mode, if present
has_sec = length(mode_idx) >= 2;
if has_sec
    lambda_sec_ref = lam0(mode_idx(2));
else
    lambda_sec_ref = NaN;
end

% --- track modes along continuation ---
for j = 1:n_cols

    lam = eigvals(:,j);

    % ---------------------------------------------------------------------
    % Dominant tracked mode
    % ---------------------------------------------------------------------
    if has_dom

        % Match eigenvalue closest to the previously tracked representative
        [~, idx_dom] = min(abs(lam - lambda_dom_ref));
        lambda_dom = lam(idx_dom);

        % Update reference for next column
        lambda_dom_ref = lambda_dom;

        % If the tracked mode is still complex, convert to frequency.
        % If it has become real, do NOT switch to another complex pair.
        if abs(imag(lambda_dom)) > imag_tol
            f_dom(j) = abs(imag(lambda_dom)) / (2*pi) * 1000;
        else
            f_dom(j) = 0; % or f_dom(j) = NaN;
        end

    end

    % ---------------------------------------------------------------------
    % Secondary tracked mode
    % ---------------------------------------------------------------------
    if has_sec

        % Match eigenvalue closest to the previously tracked secondary mode
        [~, idx_sec] = min(abs(lam - lambda_sec_ref));
        lambda_sec = lam(idx_sec);

        % Update reference for next column
        lambda_sec_ref = lambda_sec;

        % Again: if this tracked mode becomes real, do not replace it by
        % another complex pair.
        if abs(imag(lambda_sec)) > imag_tol
            f_sec(j) = abs(imag(lambda_sec)) / (2*pi) * 1000;
        else
            f_sec(j) = 0; % or f_sec(j) = NaN;
        end

    end

end

end


% function [f_dom, f_sec] = dominant_frequencies_EP(eigvals)
% 
% % DOMINANT_FREQUENCIES_EP
% % ------------------------------------------------------------
% % Convert eigenvalues -> frequencies
% %
% % IDEA:
% %   f = |Im(lambda)| / (2*pi)
% %   scaling: ms -> Hz
% 
% n_cols = size(eigvals, 2);
% f_dom = nan(1, n_cols);
% f_sec = nan(1, n_cols);
% 
% for j = 1:n_cols
%     lam = eigvals(:, j);
% 
%     % sort by closeness of real part to imaginary axis
%     [~, idx] = sort(abs(real(lam)));
% 
%     idx_dom = idx(1);
% 
%     % secondary index: choose next eigenvalue with nonzero imaginary part if possible
%     idx_sec = NaN;
%     for k = 3:length(idx)
%         if abs(imag(lam(idx(k)))) > 0
%             idx_sec = idx(k);
%             break
%         end
%     end
% 
%     omega_dom = abs(imag(lam(idx_dom)));
%     f_dom(j) = omega_dom / (2*pi) * 1000;
% 
%     if ~isnan(idx_sec)
%         omega_sec = abs(imag(lam(idx_sec)));
%         f_sec(j) = omega_sec / (2*pi) * 1000;
%     end
% end
% 
% end



%% LOMB–SCARGLE

function [f_dom, f_sec] = dominant_frequencies_LS(y, t, frequencies)

% DOMINANT_FREQUENCIES_LS
% ------------------------------------------------------------
% Frequency estimation via Lomb–Scargle periodogram.
%
% Secondary frequency is selected such that it is NOT
% a harmonic (integer multiple) of the dominant frequency.

f_dom = NaN;
f_sec = NaN;

% --- Lomb–Scargle periodogram ---
try
    power = plomb(y, t/1000, frequencies);
catch
    return
end

[pks,locs] = findpeaks(power);
if isempty(pks), return, end

% --- sort peaks by power ---
[~,idx] = sort(pks,'descend');
peak_freqs = frequencies(locs(idx));

% --- dominant frequency ---
f_dom = peak_freqs(1);

% --- parameters for harmonic rejection ---
tol = 0.05;   % tolerance (5%)

% --- find first non-harmonic frequency ---
for k = 2:length(peak_freqs)

    f_candidate = peak_freqs(k);

    ratio = f_candidate / f_dom;
    nearest_int = round(ratio);

    % check if close to integer multiple
    if abs(ratio - nearest_int) > tol
        f_sec = f_candidate;
        return
    end

end

% if all peaks are harmonics f_sec remains NaN

end


%% FINDPEAKS

function [f_dom, f_sec] = dominant_frequencies_FP(y, t, minPeakDistance_ms)

% DOMINANT_FREQUENCIES_FP
% ------------------------------------------------------------
% Frequency estimation via peak detection in time domain

f_dom = NaN;
f_sec = NaN;

signal = y - mean(y);
signal = -signal;

[~,locs] = findpeaks(signal,'MinPeakHeight',0.9*max(signal));

if numel(locs) < 2, return, end

periods = diff(t(locs));
f_dom = 1000 / mean(periods);

end


%% BRANCH EXTRACTION

function branch = extract_branch(EP, Md, Delta, selected_range, grid_size, branch_id)
%EXTRACT_BRANCH Extracts starting points and interpolated ending points.

% --- select range ---
ind_selected = Md > selected_range(1) & Md < selected_range(2);

EP = EP(:, ind_selected);
Md = Md(ind_selected);
Delta = Delta(ind_selected);

if isempty(Md)
    error('\nNo branch found within interval [%d,%d]!\n', selected_range(1), selected_range(2));
end

% --- ensure increasing Md ---
if Md(1) > Md(end)
    EP    = flip(EP, 2);
    Md    = flip(Md);
    Delta = flip(Delta);
end

% --- maximum Delta ---
[Delta_max, Delta_max_ind] = max(Delta);

% --- starting branch ---
EP_from_orig    = EP(:, 1:Delta_max_ind);
Md_from_orig    = Md(1:Delta_max_ind);
Delta_from_orig = Delta(1:Delta_max_ind);

Delta_from = linspace(0, Delta_max, grid_size(2));  % new
Md_from = interp1(Delta_from_orig, Md_from_orig, Delta_from, 'spline');
EP_from = interp1(Delta_from_orig, EP_from_orig', Delta_from, 'spline')';

% --- ending branch (original) ---
Md_to_orig    = Md(Delta_max_ind+1:end);
Delta_to_orig = Delta(Delta_max_ind+1:end);

% --- interpolate ending branch ---
if ~isempty(Md_to_orig)
    Md_to = interp1(Delta_to_orig, Md_to_orig, Delta_from, 'spline');
    Md_to(Md_to > selected_range(2)) = selected_range(2);
else
    Md_to = Md_from(end) * ones(size(Delta_from));  % no ending branch
end

Delta_to = Delta_from;

[Md_min,Md_max] = bounds([Md_from,Md_to]);
Md_grid = linspace(Md_min,Md_max,grid_size(1));

branch = struct(...
    'spec',branch_id, ...
    'range',selected_range, ...
    'grid_size',grid_size, ...
    'EP_from',EP_from, ...
    'Md_from',Md_from, ...
    'Md_to',Md_to, ...
    'Md_grid',Md_grid, ...
    'Delta_from',Delta_from, ...
    'Delta_to',Delta_to, ...
    'Delta_max',Delta_max ...
    );

end


%% COST-PROFILING HELPER FUNCTIONS

function cost = empty_cost_section()

cost = struct();

cost.method = '';
cost.bubble = '';
cost.slice_index = NaN;

% Timing
cost.wall_s = 0;              % sum of individual slice wall times
cost.grid_wall_s = 0;         % wall time of one bubble computation
cost.method_wall_time = 0;    % wall time of L + R together

% Parameter samples
cost.n_total = 0;             % total grid points considered per section
cost.n_inside = 0;            % points inside the bubble
cost.n_requested = 0;         % requested points passed to method

% Frequency estimates
cost.n_valid_dom = 0;         % valid dominant frequency estimates
cost.n_valid_sec = 0;         % valid secondary frequency estimates

% Simulation costs
cost.n_rhs = 0;               % RHS evaluations
cost.n_ode_solves = 0;        % number of ode45 calls
cost.n_time_points = 0;       % stored ODE time points

% Continuation costs
cost.n_cont_runs = 0;         % number of continuation calls
cost.n_cont_points = 0;       % continuation points
cost.n_jac_eig = 0;           % Jacobian/eigenvalue spectra
cost.n_lc_points = 0;         % LC continuation points

% Failures
cost.n_fail = 0;

end


function c = add_cost_sections(c, a)

fields = fieldnames(a);

for k = 1:numel(fields)
    fld = fields{k};

    if isfield(c, fld) && isnumeric(c.(fld)) && isnumeric(a.(fld))
        c.(fld) = c.(fld) + a.(fld);
    elseif isfield(c, fld) && isempty(c.(fld))
        c.(fld) = a.(fld);
    end
end

end


function cost = summarize_cost_sections(cost_sections)

cost = empty_cost_section();

for i = 1:numel(cost_sections)
    cost = add_cost_sections(cost, cost_sections(i));
end

if ~isempty(cost_sections)
    cost.method = cost_sections(1).method;
    cost.bubble = cost_sections(1).bubble;
end

end


function cost = combine_costs(varargin)

cost = empty_cost_section();

for i = 1:nargin
    cost = add_cost_sections(cost, varargin{i});
end

end


function print_cost_report(cost)

coverage_dom = 100 * cost.n_valid_dom / max(cost.n_inside, 1);
coverage_sec = 100 * cost.n_valid_sec / max(cost.n_inside, 1);

fprintf('\n');
fprintf('============================================================\n');
fprintf('COMPUTATIONAL COST REPORT\n');
fprintf('============================================================\n');
fprintf('Method                         : %s\n', cost.method);
fprintf('Total wall time [s]            : %.3f\n', cost.method_wall_time);
fprintf('Sum of slice wall times [s]    : %.3f\n', cost.wall_s);
fprintf('Grid wall time [s]             : %.3f\n', cost.grid_wall_s);
fprintf('\n');

fprintf('Parameter samples total        : %d\n', cost.n_total);
fprintf('Parameter samples inside       : %d\n', cost.n_inside);
fprintf('Valid dominant estimates       : %d\n', cost.n_valid_dom);
fprintf('Valid secondary estimates      : %d\n', cost.n_valid_sec);
fprintf('Coverage dominant [%%]          : %.2f\n', coverage_dom);
fprintf('Coverage secondary [%%]         : %.2f\n', coverage_sec);
fprintf('\n');

fprintf('ODE solves                     : %d\n', cost.n_ode_solves);
fprintf('RHS evaluations                : %d\n', cost.n_rhs);
fprintf('Stored ODE time points         : %d\n', cost.n_time_points);
fprintf('\n');

fprintf('Continuation runs              : %d\n', cost.n_cont_runs);
fprintf('Continuation points            : %d\n', cost.n_cont_points);
fprintf('Jacobian/eigenvalue spectra    : %d\n', cost.n_jac_eig);
fprintf('Limit-cycle points             : %d\n', cost.n_lc_points);
fprintf('\n');

fprintf('Failed sections / operations   : %d\n', cost.n_fail);
fprintf('============================================================\n');

fprintf('\nCompact table row:\n');
fprintf('method | wall_s | worker_s | N_inside | N_valid | RHS_eval | cont_pts | eig_specs | coverage\n');
fprintf('%-8s | %10.2f | %10.2f | %8d | %8d | %10d | %10d | %10d | %8.2f\n', ...
    cost.method, ...
    cost.method_wall_time, ...
    cost.wall_s, ...
    cost.n_inside, ...
    cost.n_valid_dom, ...
    cost.n_rhs, ...
    cost.n_cont_points, ...
    cost.n_jac_eig, ...
    coverage_dom);

fprintf('\n');

end