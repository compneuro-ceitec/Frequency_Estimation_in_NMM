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


%% USER SETTINGS
target_band = [5 25];        % [Hz]
target_freq_mode = 'sec';    % 'dom' = dominant frequency, 'sec' = secondary frequency

% --- par1 & par2 ---
par1_range = [0 20]; par1_idx = 4; par1_label = '$M_{D1}$';                 % Md
par2_range = [0 6];  par2_idx = 9; par2_label = '$\Delta_{v_\theta}$';      % Delta

% --- par3 ---
par3_range = [1 200];   par3_idx = 1;  par3_label = '$C$';                  % C
% par3_range = [0 70];    par3_idx = 2;  par3_label = '$I_{ext}$';            % Iext
% par3_range = [0 120];   par3_idx = 3;  par3_label = '$u_{jump}$';           % ujump
% par3_range = [0 1];     par3_idx = 5;  par3_label = '$B_{M_{D1}}$';         % B
% par3_range = [1 3];     par3_idx = 6;  par3_label = '$\gamma$';             % gamma
% par3_range = [-80 -50]; par3_idx = 7;  par3_label = '$v_r$';                % Vr
% par3_range = [-55 -35]; par3_idx = 8;  par3_label = '$\bar{v}_{\theta}$';   % Vt
% par3_range = [1 10];    par3_idx = 11; par3_label = '$\tau_{s_a}$';         % Ta

n_par3 = 15;            % one point per core
n_refine_iters = 1;     % number of boundary-refinement iterations
par3_grid = linspace(par3_range(1), par3_range(2), n_par3);

% --- names / indices of parameters in pars vector ---
cfg = struct();
cfg.par1_idx = par1_idx;
cfg.par2_idx = par2_idx;
cfg.par3_idx = par3_idx;

cfg.par1_range       = par1_range;
cfg.par2_range       = par2_range;
cfg.par3_range       = par3_range;
cfg.target_band      = target_band;
cfg.target_freq_mode = target_freq_mode;

cfg.init_time_step       = 1;
cfg.max_points_ep        = 100;
cfg.max_points_hopf_test = 10;
cfg.max_points_hopf_full = 50;

% --- fixed/default parameter vector ---
cfg.pars = [...
    50, ...   % C
    0,  ...   % Iext
    45, ...   % ujump
    0,  ...   % Md
    0.2,...   % B
    1.1,...   % gammaa
    -60,...   % Vr                                                 < Vt (!)
    -43,...   % Vt
    0,  ...   % Delta
    0,  ...   % cexc
    6,  ...   % Ta
    10  ...   % K
];

% --- initial equilibrium guess ---
cfg.x0_eq = [0.006 -70 -50 0.06]';


%% COMPUTATION
tic;

% --- preallocation ---
slice_results = cell(n_par3,1);

% --- parfor loop ---
parfor i = 1:n_par3
    par3_value = par3_grid(i);
    slice_results{i} = compute_hopf_slice(par3_value, cfg);
end

% --- collect summary ---
summary = struct( ...
    'par3', num2cell(nan(n_par3,1)), ...
    'f_min', num2cell(nan(n_par3,1)), ...
    'f_max', num2cell(nan(n_par3,1)), ...
    'f_mean', num2cell(nan(n_par3,1)), ...
    'f_std', num2cell(nan(n_par3,1)), ...
    'f2_min', num2cell(nan(n_par3,1)), ...
    'f2_max', num2cell(nan(n_par3,1)), ...
    'f2_mean', num2cell(nan(n_par3,1)), ...
    'f2_std', num2cell(nan(n_par3,1)), ...
    'distance', num2cell(nan(n_par3,1)), ...
    'feasible', num2cell(false(n_par3,1)), ...
    'num_points', num2cell(zeros(n_par3,1)), ...
    'num_branches', num2cell(zeros(n_par3,1)), ...
    'status', repmat({''}, n_par3,1), ...
    'branches', repmat({[]}, n_par3,1));

all_points = [];        % [par1, par2, par3, f_dom, f_sec, branch_id]
all_branch_stats = [];  % [par3, branch_id, f_min, f_mean, f_max, f_std, f2_min, f2_mean, f2_max, f2_std, distance, feasible]

for i = 1:n_par3
    s = slice_results{i};

    summary(i).par3       = s.par3;
    summary(i).f_min      = s.f_min;
    summary(i).f_max      = s.f_max;
    summary(i).f_mean     = s.f_mean;
    summary(i).f_std      = s.f_std;
    summary(i).f2_min     = s.f2_min;
    summary(i).f2_max     = s.f2_max;
    summary(i).f2_mean    = s.f2_mean;
    summary(i).f2_std     = s.f2_std;
    summary(i).distance   = s.distance;
    summary(i).feasible   = s.feasible;
    summary(i).num_points = s.num_points;
    summary(i).status     = s.status;
    summary(i).branches   = s.branches;

    if isempty(s.points)
        summary(i).num_branches = 0;
    else
        summary(i).num_branches = numel(unique(s.points(:,end)));
        all_points = [all_points; s.points];
    end

    if ~isempty(s.branches)
        for j = 1:numel(s.branches)
            bj = s.branches(j);
            all_branch_stats = [all_branch_stats; ...
                s.par3, bj.branch_id, ...
                bj.f_min, bj.f_mean, bj.f_max, bj.f_std, ...
                bj.f2_min, bj.f2_mean, bj.f2_max, bj.f2_std, ...
                bj.distance, bj.feasible];
        end
    end
end

% --- sorting ---
[~,ord] = sort([summary.par3]);
summary = summary(ord);


%% REFINEMENT
[left_bracket, right_bracket, refine_history] = ...
    refine_par3_bounds(summary, cfg, n_par3, n_refine_iters);


%% SUMMARY
fprintf('\n=== COARSE GRID COMPUTATION SUMMARY ===\n');
fprintf('Target band: [%.2f, %.2f] Hz, target frequency mode: %s\n\n', ...
    cfg.target_band(1), cfg.target_band(2), cfg.target_freq_mode);

rowfmt = '%-15s # %-14s | %-6s | %-6s # %-14s | %-6s | %-6s # %-8s | %-8s | %-8s | %-20s\n';

% --- table header ---
header1 = sprintf(['%-15s # %-32s # %-32s # \n'], ...
    '', '       DOMINANT FREQUENCY', '       SECONDARY FREQUENCY');

header2 = sprintf(rowfmt, ...
    'item', '[f_min, f_max]', 'f_mean', 'f_std', '[f_min, f_max]', ...
    'f_mean', 'f_std', 'distance', 'feasible', 'n_points', 'status');

fprintf('%s\n', repelem('-', length(header2)-1));
fprintf('%s', header1);
fprintf('%s', header2);
fprintf('%s\n', repelem('=', length(header2)-1));

for i = 1:numel(summary)

    % --- slice row ---
    slice_name = sprintf('par3 = %.3g', summary(i).par3);

    fprintf(rowfmt, ...
        slice_name, ...
        fmt_interval(summary(i).f_min, summary(i).f_max), ...
        fmt_num1(summary(i).f_mean), ...
        fmt_num1(summary(i).f_std), ...
        fmt_interval(summary(i).f2_min, summary(i).f2_max), ...
        fmt_num1(summary(i).f2_mean), ...
        fmt_num1(summary(i).f2_std), ...
        fmt_num1(summary(i).distance), ...
        fmt_logical(summary(i).feasible), ...
        sprintf(' %d / %d', summary(i).num_points, summary(i).num_branches), ...
        fmt_string(summary(i).status));

    % --- branch rows ---
    if ~isempty(summary(i).branches)
        for j = 1:numel(summary(i).branches)
            bj = summary(i).branches(j);

            if isfield(bj,'f') && ~isempty(bj.f)
                branch_n_points = numel(bj.f);
            else
                branch_n_points = [];
            end

            fprintf(rowfmt, ...
                sprintf('     branch %d', bj.branch_id), ...
                fmt_interval(bj.f_min, bj.f_max), ...
                fmt_num1(bj.f_mean), ...
                fmt_num1(bj.f_std), ...
                fmt_interval(bj.f2_min, bj.f2_max), ...
                fmt_num1(bj.f2_mean), ...
                fmt_num1(bj.f2_std), ...
                fmt_num1(bj.distance), ...
                fmt_logical(bj.feasible), ...
                fmt_int(branch_n_points), ...
                '');
        end
    end

    % --- separator between slices ---
    fprintf('%s\n', repelem('-', length(header2)-1));
end

feas = [summary.feasible];
p3   = [summary.par3];

idx = find(feas);
if isempty(idx)
    fprintf('\nNo feasible par3 point found on coarse grid.\n\n');
else
    breaks = [1, find(diff(idx) > 1) + 1, length(idx)+1];
    fprintf('\nFeasible coarse intervals in par3:    ');
    for b = 1:length(breaks)-1
        block = idx(breaks(b):breaks(b+1)-1);
        fprintf('[%.4g, %.4g]\n\n', p3(block(1)), p3(block(end)));
    end
end

% REFINEMENT
if ~isempty(left_bracket) || ~isempty(right_bracket)
    fprintf('=== REFINEMENT SUMMARY ===\n')
end

if ~isempty(left_bracket)
    fprintf('Lower-bound bracket in par3:   [%.4g, %.4g]\n', ...
        left_bracket(1), left_bracket(2));
end

if ~isempty(right_bracket)
    fprintf('Upper-bound bracket in par3:   [%.4g, %.4g]\n', ...
        right_bracket(1), right_bracket(2));
end

if ~isempty(left_bracket) && ~isempty(right_bracket)
    fprintf('Estimated par3 interval:       [%.4g, %.4g]\n', ...
        mean(left_bracket), mean(right_bracket));

    m = max(diff(left_bracket), diff(right_bracket));
    fprintf(['Upper bound on bracket width:  %.4g ', ...
        '(~%.2f%% of par3range span)\n\n'], m, 100 * m / range(cfg.par3_range));
end

toc;


%% VISUALIZATION
titles_on = true;

% SCATTER (all_points: [par1, par2, par3, f_dom, f_sec, branch_id]) -------
if ~isempty(all_points)
    fig = figure;
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[11,4.5])

    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % --- dominant frequency plot ---
    ax1 = nexttile;
    scatter3(ax1, all_points(:,1), all_points(:,2), all_points(:,3), ...
             20, all_points(:,4), 'filled');

    xlabel(ax1, par1_label);
    ylabel(ax1, par2_label);
    zlabel(ax1, par3_label);
    if titles_on
        title(ax1, 'DOMINANT FREQUENCY');
    end

    box(ax1, 'on');
    grid(ax1, 'on');
    xlim(ax1, cfg.par1_range);
    ylim(ax1, cfg.par2_range);
    zlim(ax1, cfg.par3_range);
    xtickangle(ax1, 0);
    ytickangle(ax1, 0);
    ztickangle(ax1, 0);
    view(ax1, -15, 6);

    cb1 = colorbar(ax1);
    cb1.Label.String = 'Frequency [Hz]';
    cb1.TickLabelInterpreter = 'latex';
    cb1.Label.Interpreter = 'latex';

    % --- secondary frequency plot ---
    ax2 = nexttile;
    scatter3(ax2, all_points(:,1), all_points(:,2), all_points(:,3), ...
             20, all_points(:,5), 'filled');

    xlabel(ax2, par1_label);
    ylabel(ax2, par2_label);
    zlabel(ax2, par3_label);
    if titles_on
        title(ax2, 'SECONDARY FREQUENCY');
    end

    box(ax2, 'on');
    grid(ax2, 'on');
    xlim(ax2, cfg.par1_range);
    ylim(ax2, cfg.par2_range);
    zlim(ax2, cfg.par3_range);
    xtickangle(ax2, 0);
    ytickangle(ax2, 0);
    ztickangle(ax2, 0);
    view(ax2, -15, 6);

    cb2 = colorbar(ax2);
    cb2.Label.String = 'Frequency [Hz]';
    cb2.TickLabelInterpreter = 'latex';
    cb2.Label.Interpreter = 'latex';
end

% DIAGNOSTIC PLOT ---------------------------------------------------------
if ~isempty(all_branch_stats)
    p3_branch       = all_branch_stats(:,1);
    fmin_branch     = all_branch_stats(:,3);
    fmean_branch    = all_branch_stats(:,4);
    fmax_branch     = all_branch_stats(:,5);
    fstd_branch     = all_branch_stats(:,6); 
    f2min_branch    = all_branch_stats(:,7);
    f2mean_branch   = all_branch_stats(:,8);
    f2max_branch    = all_branch_stats(:,9);
    f2std_branch    = all_branch_stats(:,10); 
    dist_branch     = all_branch_stats(:,11);
    feasible_branch = logical(all_branch_stats(:,12));

    fig_diag = figure;
    set(fig_diag,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[11,4.5])

    tl_diag = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % --- dominant frequency diagnostic ---
    axd1 = nexttile;
    hold(axd1, 'on');

    % Variant A: highlight target only in the relevant panel
    if strcmpi(cfg.target_freq_mode, 'dom')
        draw_feasible_band_patch(axd1, summary, cfg, left_bracket, right_bracket);

        h = yline(axd1, cfg.target_band(1), '--k');
        h.Label = '$f^*_{\min}$';
        h.Interpreter = 'latex';

        h = yline(axd1, cfg.target_band(2), '--k');
        h.Label = '$f^*_{\max}$';
        h.Interpreter = 'latex';
    end

    lower_err = fmean_branch - fmin_branch;
    upper_err = fmax_branch - fmean_branch;

    errorbar(axd1, p3_branch, fmean_branch, lower_err, upper_err, '.', ...
        'LineStyle','none', 'MarkerSize',15);

    xlabel(axd1, par3_label);
    ylabel(axd1, 'Frequency [Hz]');
    if titles_on
        title(axd1, 'DOMINANT FREQUENCY');
    end

    box(axd1, 'on');
    grid(axd1, 'off');
    xlim(axd1, cfg.par3_range);

    % --- secondary frequency diagnostic ---
    axd2 = nexttile;
    hold(axd2, 'on');

    % Variant A: highlight target only in the relevant panel
    if strcmpi(cfg.target_freq_mode, 'sec')
        draw_feasible_band_patch(axd2, summary, cfg, left_bracket, right_bracket);

        h = yline(axd2, cfg.target_band(1), '--k');
        h.Label = '$f^*_{\min}$';
        h.Interpreter = 'latex';

        h = yline(axd2, cfg.target_band(2), '--k');
        h.Label = '$f^*_{\max}$';
        h.Interpreter = 'latex';
    end

    lower_err2 = f2mean_branch - f2min_branch;
    upper_err2 = f2max_branch - f2mean_branch;

    errorbar(axd2, p3_branch, f2mean_branch, lower_err2, upper_err2, '.', ...
        'LineStyle','none', 'MarkerSize',15);

    xlabel(axd2, par3_label);
    ylabel(axd2, 'Frequency [Hz]');
    if titles_on
        title(axd2, 'SECONDARY FREQUENCY');
    end

    box(axd2, 'on');
    grid(axd2, 'off');
    xlim(axd2, cfg.par3_range);

    % --- unify y-axis limits across both panels ---
    y_all = [ ...
        fmin_branch(:); fmax_branch(:); ...
        f2min_branch(:); f2max_branch(:) ...
    ];

    % include target band so it remains visible when active
    y_all = [y_all; cfg.target_band(:)];
    y_all = y_all(~isnan(y_all));

    if isempty(y_all)
        ylim_shared = [0 1];
    else
        y_min = min(y_all);
        y_max = max(y_all);

        if y_min == y_max
            ylim_shared = [y_min, y_min + 1];
        else
            pad = 0.05 * (y_max - y_min);
            ylim_shared = [max([0, y_min - pad]), y_max + pad];
        end
    end

    ylim(axd1, ylim_shared);
    ylim(axd2, ylim_shared);
end

% OVERLAP DISTANCE --------------------------------------------------------
if ~isempty(all_branch_stats)
    figure;
    hold on;
    scatter(p3_branch(~feasible_branch), dist_branch(~feasible_branch), ...
        15, 'r', 'filled');
    scatter(p3_branch(feasible_branch), dist_branch(feasible_branch), ...
        15, 'g', 'filled');
    yline(0, '--k');
    xlabel(par3_label);
    ylabel('Overlap distance');
    if titles_on
        title(['Branch-dependent overlap distance vs ', par3_label, ...
            ' (', upper(cfg.target_freq_mode), ')']);
    end
    legend('infeasible', 'feasible', 'Location', 'best');

    box on; grid off;
    xlim(cfg.par3_range);
end


%% SAVE
% save('hopf_scan_par3.mat', 'cfg', 'par3_grid', ...
%     'slice_results', 'summary', 'all_points', 'all_branch_stats', ...
%     'left_bracket', 'right_bracket', 'refine_history');

warning('on')


%% COMPUTE HOPF SLICE
function slice = compute_hopf_slice(par3_value, cfg)
% COMPUTE_HOPF_SLICE Compute one Hopf slice for fixed par3.
%
% INPUT:
%   par3_value                  - value of parameter par3
%   cfg                         - structure with slice settings:
%     .par1_idx / .par2_idx     - index of parameter par1 (par2) in cfg.pars
%     .par3_idx                 - index of parameter par3 in cfg.pars
%     .par1_range / .par2_range - allowed interval [min max] for par1 (par2)
%     .par3_range               - scanned interval [min max] for par3
%     .target_band              - target frequency band [f_min, f_max] [Hz]
%     .target_freq_mode         - which frequency is compared to target_band:
%                                 'dom' or 'sec'
%     .init_time_step           - initial continuation step size
%     .max_points_ep            - max. number of equilibrium-continuation points
%     .max_points_hopf_test     - max. number of test Hopf-continuation points
%     .max_points_hopf_full     - max. number of full Hopf-continuation points
%     .pars                     - default parameter vector
%                                 [C, Iext, ujump, Md, B, gamma, Vr, Vt, Delta, cexc, Ta, K]
%     .x0_eq                    - initial equilibrium guess
%
% OUTPUT:
%   slice                       - structure representing the computed slice:
%     .par3                     - current par3 value
%     .f_min / .f2_min          - minimum dominant (secondary) frequency [Hz]
%     .f_max / .f2_max          - maximum dom. (sec.) frequency [Hz]
%     .f_mean / .f2_mean        - mean dom. (sec.) frequency [Hz]
%     .f_std / .f2_std          - std. dev. of dom. (sec.) frequency [Hz]
%     .distance                 - minimum overlap distance to target_band
%     .feasible                 - logical flag, true if some branch is feasible
%     .num_points               - number of valid Hopf points in the box
%     .status                   - computation status string
%                                 ('init', 'ok', 'no_hopf_on_EP',
%                                  'no_valid_points_in_box', 'failed: ...')
%     .points                   - valid points, columns [par1, par2, par3, f_dom, f_sec, branch_id]
%     .branches                 - structure array of valid Hopf branches:
%       .branch_id              - integer identifier
%       .par1 / .par2           - branch values of par1 (par2)
%       .f / .f2                - dom. (sec.) frequencies [Hz]
%       .f_min / .f2_min        - minimum dom. (sec.) frequency [Hz]
%       .f_max / .f2_max        - maximum dom. (sec.) frequency [Hz]
%       .f_mean / .f2_mean      - mean dom. (sec.) frequency [Hz]
%       .f_std / .f2_std        - std. dev. of dom. (sec.) frequency [Hz]
%       .distance               - overlap distance of selected branch interval to target_band
%       .feasible               - logical flag for branch feasibility

    slice            = struct();
    slice.par3       = par3_value;
    slice.f_min      = NaN;
    slice.f_max      = NaN;
    slice.f_mean     = NaN;
    slice.f_std      = NaN;
    slice.f2_min     = NaN;
    slice.f2_max     = NaN;
    slice.f2_mean    = NaN;
    slice.f2_std     = NaN;
    slice.distance   = Inf;
    slice.feasible   = false;
    slice.num_points = 0;
    slice.status     = 'init';
    slice.points     = [];   % [par1, par2, par3, f_dom, f_sec, branch_id]
    slice.branches   = [];

    try
        % SET PARAMETERS
        pars = cfg.pars;
        pars(cfg.par3_idx) = par3_value;

        par1_idx = cfg.par1_idx;
        par2_idx = cfg.par2_idx;

        hls = Dopamine_NMM_fast4;
        odefun = @(t,y) hls{2}(t, y, pars(1), pars(2), pars(3), ...
            pars(4), pars(5), pars(6), pars(7), pars(8), pars(9), ...
            pars(10), pars(11), pars(12));

        % EQUILIBRIUM
        opts = optimoptions('fsolve', ...
            'Display','off', ...
            'FunctionTolerance',1e-12, ...
            'StepTolerance',1e-12, ...
            'OptimalityTolerance',1e-12);

        EP = fsolve(@(Y) odefun(0,Y), cfg.x0_eq, opts);

        % EP CONTINUATION IN par1
        apEP = par1_idx;
        [x0, v0] = init_EP_EP(@Dopamine_NMM_fast4, EP, pars, apEP);

        opt = contset;
        opt = contset(opt, 'InitStepsize', cfg.init_time_step);
        opt = contset(opt, 'MinStepsize', 0.001);
        opt = contset(opt, 'MaxStepsize', 10);
        opt = contset(opt, 'MaxNumPoints', cfg.max_points_ep);
        opt = contset(opt, 'Adapt', 3);
        opt = contset(opt, 'Singularities', 1);
        opt = contset(opt, 'Eigenvalues', 1);

        [x_ep, v_ep, s_ep, h_ep, f_ep] = cont(@equilibrium, x0, [], opt); 

        idx_hopf = find_hopf_point(s_ep);

        if isempty(idx_hopf)
            slice.status = 'no_hopf_on_EP';
            return;
        end

        % LOOP OVER HOPF POINTS
        valid_points = [];
        valid_f2_points = [];
        branch_counter = 0;
        valid_branches = struct( ...
            'branch_id', {}, ...
            'par1', {}, ...
            'par2', {}, ...
            'f', {}, ...
            'f2', {}, ...
            'f_min', {}, ...
            'f_max', {}, ...
            'f_mean', {}, ...
            'f_std', {}, ...
            'f2_min', {}, ...
            'f2_max', {}, ...
            'f2_mean', {}, ...
            'f2_std', {}, ...
            'distance', {}, ...
            'feasible', {});

        for k = 1:length(idx_hopf)

            idx = idx_hopf(k);
        
            pars(apEP) = x_ep(end, idx);
        
            % =============================================================
            % Either compute the Hopf branch, or use only the detected 
            % Hopf point from EP continuation.
            % =============================================================
            if cfg.max_points_hopf_full <= 1

                % ---------------------------------------------------------
                % FAST MODE: only the detected Hopf point on EP continuation
                % ---------------------------------------------------------
                par1 = x_ep(end, idx);
                par2 = pars(par2_idx);
        
                [f_dom, f_sec] = dominant_frequencies_EP(f_ep(:, idx));
        
            else
                
                % ---------------------------------------------------------
                % FULL MODE: H-H continuation
                % ---------------------------------------------------------
                apH = [par1_idx, par2_idx];
        
                [xH, vH] = init_H_H(@Dopamine_NMM_fast4, ...
                                    x_ep(1:end-1, idx), pars', apH);
        
                % --- short test continuation to choose direction toward par2 > 0 ---
                opt_test = contset;
                opt_test = contset(opt_test, 'MaxNumPoints', cfg.max_points_hopf_test);
                opt_test = contset(opt_test, 'InitStepsize', cfg.init_time_step);
                opt_test = contset(opt_test, 'MinStepsize', 0.001);
                opt_test = contset(opt_test, 'MaxStepsize', 10);
                opt_test = contset(opt_test, 'Singularities', 0);
        
                % forward
                opt_f = contset(opt_test, 'Backward', 0);
                [xHf,~,~,~,~] = cont(@hopf, xH, vH, opt_f);
        
                par2_f = xHf(end-1,:);
                n_f = length(par2_f);
                w_f = max(1, n_f-5):n_f;
                cond_f = (mean(par2_f(w_f)) > 0);
        
                % backward
                opt_b = contset(opt_test, 'Backward', 1);
                [xHb,~,~,~,~] = cont(@hopf, xH, vH, opt_b);
        
                par2_b = xHb(end-1,:);
                n_b = length(par2_b);
                w_b = max(1, n_b-5):n_b;
                cond_b = (mean(par2_b(w_b)) > 0);
        
                % continuation toward positive values of par2
                if cond_f
                    selected_dir = 0;
                elseif cond_b
                    selected_dir = 1;
                else
                    continue
                end
        
                % --- full continuation ---
                opt_full = contset(opt_test, 'MaxNumPoints', cfg.max_points_hopf_full);
                opt_full = contset(opt_full, 'Backward', selected_dir);
                opt_full = contset(opt_full, 'Eigenvalues', 1);
        
                [xH_full, ~, ~, ~, fH_full] = ...
                    cont(@hopf, xH, vH, opt_full);
        
                [f_dom, f_sec] = dominant_frequencies_EP( ...
                    fH_full(:,1:size(xH_full,2)));
        
                par1 = xH_full(end-2,:);
                par2 = xH_full(end-1,:);
        
            end
        
            % =============================================================
            % Common processing for both FULL and FAST mode
            % =============================================================
        
            % restrict to desired par1_range x par2_range box
            idx_box = (par1 >= cfg.par1_range(1)) & (par1 <= cfg.par1_range(2)) & ...
                      (par2 >= cfg.par2_range(1)) & (par2 <= cfg.par2_range(2)) & ...
                      (f_dom > 0);
        
            par1_sel = par1(idx_box);
            par2_sel = par2(idx_box);
            f_sel    = f_dom(idx_box);
            f2_sel   = f_sec(idx_box);
        
            if isempty(par1_sel)
                continue
            end
        
            branch_counter = branch_counter + 1;
        
            pts = [par1_sel(:), ...
                   par2_sel(:), ...
                   par3_value * ones(numel(par1_sel),1), ...
                   f_sel(:), ...
                   f2_sel(:), ...
                   branch_counter * ones(numel(par1_sel),1)];
        
            [f_min, f_max, f_mean, f_std]     = freq_stats(f_sel);
            [f2_min, f2_max, f2_mean, f2_std] = freq_stats(f2_sel);
        
            branch_struct = struct();
            branch_struct.branch_id = branch_counter;
            branch_struct.par1      = par1_sel(:);
            branch_struct.par2      = par2_sel(:);
            branch_struct.f         = f_sel(:);
            branch_struct.f2        = f2_sel(:);
            branch_struct.f_min     = f_min;
            branch_struct.f_max     = f_max;
            branch_struct.f_mean    = f_mean;
            branch_struct.f_std     = f_std;
            branch_struct.f2_min    = f2_min;
            branch_struct.f2_max    = f2_max;
            branch_struct.f2_mean   = f2_mean;
            branch_struct.f2_std    = f2_std;
        
            % choose which frequency interval is compared to target_band
            switch lower(string(cfg.target_freq_mode))
                case "dom"
                    f_low_target  = branch_struct.f_min;
                    f_high_target = branch_struct.f_max;
                case "sec"
                    f_low_target  = branch_struct.f2_min;
                    f_high_target = branch_struct.f2_max;
                otherwise
                    error('Invalid cfg.target_freq_mode. Use ''dom'' or ''sec''.');
            end
        
            branch_struct.distance = overlap_distance(f_low_target, ...
                f_high_target, cfg.target_band);
            branch_struct.feasible = (branch_struct.distance <= 0);
        
            valid_branches(end+1) = branch_struct;
            valid_points = [valid_points; pts];
            valid_f2_points = [valid_f2_points; f2_sel(:)];
        
        end

        if isempty(valid_points)
            slice.status = 'no_valid_points_in_box';
            return;
        end

        % FREQUENCY INTERVAL
        branch_distances = [valid_branches.distance];
        slice_distance = min(branch_distances);

        [f_min, f_max, f_mean, f_std]     = freq_stats(valid_points(:,4));
        [f2_min, f2_max, f2_mean, f2_std] = freq_stats(valid_points(:,5));

        slice.points     = valid_points;
        slice.branches   = valid_branches;
        slice.num_points = size(valid_points,1);
        slice.f_min      = f_min;
        slice.f_max      = f_max;
        slice.f_mean     = f_mean;
        slice.f_std      = f_std;
        slice.f2_min     = f2_min;
        slice.f2_max     = f2_max;
        slice.f2_mean    = f2_mean;
        slice.f2_std     = f2_std;
        slice.distance   = slice_distance;
        slice.feasible   = any([valid_branches.feasible]);
        slice.status     = 'ok';

    catch ME
        slice.status = ['failed: ' ME.identifier];
    end
end


%% OVERLAP DISTANCE
function d = overlap_distance(f_min_slice, f_max_slice, target_band)
% OVERLAP_DISTANCE Compute overlap distance between [f_min_slice, f_max_slice] and target_band.
%
% INPUT:
%   f_min_slice / f_max_slice - lower / upper bound of slice interval [Hz]
%   target_band               - target interval [f_target_min, f_target_max] [Hz]
%
% OUTPUT:
%   d                         - overlap distance; d <= 0 means overlap

    f_target_min = target_band(1);
    f_target_max = target_band(2);

    d = max([f_target_min - f_max_slice, ...
             f_min_slice - f_target_max]);
end


%% REFINEMENT
function [left_bracket, right_bracket, refine_history] = ...
    refine_par3_bounds(summary, cfg, n_par3, n_refine_iters)
% REFINE_PAR3_BOUNDS Refine lower / upper feasible boundary in par3.
%
% INPUT:
%   summary                         - coarse-grid summary structure
%     .par3 / .feasible             - tested par3 values / feasibility flags
%   cfg                             - settings passed to compute_hopf_slice
%   n_par3                          - number of inner points per refinement step
%   n_refine_iters                  - number of refinement iterations
%
% OUTPUT:
%   left_bracket / right_bracket    - lower / upper boundary bracket [a, b]
%   refine_history                  - structure with refinement summaries:
%     .left / .right                - cell arrays of intermediate summaries

    % INITIAL BRACKETS
    left_info  = find_left_boundary_info(summary);
    right_info = find_right_boundary_info(summary);

    if isempty(left_info)
        left_bracket = [];
    else
        left_bracket = left_info.bracket;
    end

    if isempty(right_info)
        right_bracket = [];
    else
        right_bracket = right_info.bracket;
    end

    % PREALLOCATE HISTORY
    refine_history = struct();
    refine_history.left = cell(n_refine_iters,1);
    refine_history.right = cell(n_refine_iters,1);

    % ITERATIVE REFINEMENT
    for it = 1:n_refine_iters

        % LEFT BOUNDARY
        if ~isempty(left_info) && diff(left_info.bracket) > 0
            a = left_info.bracket(1);
            b = left_info.bracket(2);

            left_grid_full  = linspace(a, b, n_par3+2);
            left_grid_inner = left_grid_full(2:end-1);

            left_summary_inner = evaluate_par3_grid_summary(left_grid_inner, cfg);

            left_summary_full = merge_boundary_and_inner_summary( ...
                left_info.left_summary, left_summary_inner, left_info.right_summary);

            left_info = find_left_boundary_info(left_summary_full);
            left_bracket = left_info.bracket;
            refine_history.left{it} = left_summary_full;
        end

        % RIGHT BOUNDARY
        if ~isempty(right_info) && diff(right_info.bracket) > 0
            a = right_info.bracket(1);
            b = right_info.bracket(2);

            right_grid_full  = linspace(a, b, n_par3+2);
            right_grid_inner = right_grid_full(2:end-1);

            right_summary_inner = evaluate_par3_grid_summary(right_grid_inner, cfg);

            right_summary_full = merge_boundary_and_inner_summary( ...
                right_info.left_summary, right_summary_inner, right_info.right_summary);

            right_info = find_right_boundary_info(right_summary_full);
            right_bracket = right_info.bracket;
            refine_history.right{it} = right_summary_full;
        end
    end
end


function summary = evaluate_par3_grid_summary(par3_grid, cfg)
% EVALUATE_PAR3_GRID_SUMMARY Evaluate compute_hopf_slice on a par3 grid.
%
% INPUT:
%   par3_grid                   - vector of tested par3 values
%   cfg                         - settings passed to compute_hopf_slice
%
% OUTPUT:
%   summary                     - structure array sorted by par3:
%     .par3                     - tested par3 value
%     .f_min / .f2_min          - minimum dom. (sec.) frequency [Hz]
%     .f_max / .f2_max          - maximum dom. (sec.) frequency [Hz]
%     .f_mean / .f2_mean        - mean dom. (sec.) frequency [Hz]
%     .f_std / .f2_std          - std. dev. of dom. (sec.) frequency [Hz]
%     .distance                 - minimum overlap distance to target_band
%     .feasible                 - logical feasibility flag
%     .num_points               - number of valid points
%     .num_branches             - number of valid branches
%     .status                   - computation status string
%     .branches                 - branch structure array

    % PREALLOCATE
    n_pts = numel(par3_grid);
    slice_results_local = cell(n_pts,1);

    % EVALUATION
    parfor i = 1:n_pts
        par3_value = par3_grid(i);
        slice_results_local{i} = compute_hopf_slice(par3_value, cfg);
    end

    % COLLECT SUMMARY
    summary = struct( ...
        'par3', num2cell(nan(n_pts,1)), ...
        'f_min', num2cell(nan(n_pts,1)), ...
        'f_max', num2cell(nan(n_pts,1)), ...
        'f_mean', num2cell(nan(n_pts,1)), ...
        'f_std', num2cell(nan(n_pts,1)), ...
        'f2_min', num2cell(nan(n_pts,1)), ...
        'f2_max', num2cell(nan(n_pts,1)), ...
        'f2_mean', num2cell(nan(n_pts,1)), ...
        'f2_std', num2cell(nan(n_pts,1)), ...
        'distance', num2cell(nan(n_pts,1)), ...
        'feasible', num2cell(false(n_pts,1)), ...
        'num_points', num2cell(zeros(n_pts,1)), ...
        'num_branches', num2cell(zeros(n_pts,1)), ...
        'status', repmat({''}, n_pts,1), ...
        'branches', repmat({[]}, n_pts,1));

    for i = 1:n_pts
        s = slice_results_local{i};

        summary(i).par3       = s.par3;
        summary(i).f_min      = s.f_min;
        summary(i).f_max      = s.f_max;
        summary(i).f_mean     = s.f_mean;
        summary(i).f_std      = s.f_std;
        summary(i).f2_min     = s.f2_min;
        summary(i).f2_max     = s.f2_max;
        summary(i).f2_mean    = s.f2_mean;
        summary(i).f2_std     = s.f2_std;
        summary(i).distance   = s.distance;
        summary(i).feasible   = s.feasible;
        summary(i).num_points = s.num_points;
        summary(i).status     = s.status;
        summary(i).branches   = s.branches;

        if isempty(s.points)
            summary(i).num_branches = 0;
        else
            summary(i).num_branches = numel(unique(s.points(:,end)));
        end
    end

    % SORT BY par3
    [~,ord] = sort([summary.par3]);
    summary = summary(ord);
end


function summary_full = merge_boundary_and_inner_summary(left_summary, ...
    inner_summary, right_summary)
% MERGE_BOUNDARY_AND_INNER_SUMMARY Merge boundary points with inner refinement points.
%
% INPUT:
%   left_summary / right_summary - already computed boundary summaries
%   inner_summary                - newly computed inner summaries
%
% OUTPUT:
%   summary_full                 - merged summary sorted by par3

    % MERGE
    if isempty(inner_summary)
        summary_full = [left_summary; right_summary];
    else
        summary_full = [left_summary; inner_summary; right_summary];
    end

    % SORT BY par3
    [~,ord] = sort([summary_full.par3]);
    summary_full = summary_full(ord);
end


function info = find_left_boundary_info(summary)
% FIND_LEFT_BOUNDARY_INFO Find first transition to feasible region.
%
% INPUT:
%   summary                      - summary structure with .par3 / .feasible
%
% OUTPUT:
%   info                         - boundary-info structure:
%     .bracket                   - lower boundary bracket [a, b]
%     .left_summary              - left bracket endpoint summary
%     .right_summary             - right bracket endpoint summary

    % FIND LEFT BOUNDARY
    info = [];

    feas = [summary.feasible];
    p3   = [summary.par3];

    if ~any(feas)
        return;
    end

    first_feasible = find(feas, 1, 'first');

    if first_feasible == 1
        info.bracket = [p3(1), p3(1)];
        info.left_summary  = summary(1);
        info.right_summary = summary(1);
    else
        info.bracket = [p3(first_feasible-1), p3(first_feasible)];
        info.left_summary  = summary(first_feasible-1);
        info.right_summary = summary(first_feasible);
    end
end


function info = find_right_boundary_info(summary)
% FIND_RIGHT_BOUNDARY_INFO Find last transition from feasible region.
%
% INPUT:
%   summary                      - summary structure with .par3 / .feasible
%
% OUTPUT:
%   info                         - boundary-info structure:
%     .bracket                   - upper boundary bracket [a, b]
%     .left_summary              - left bracket endpoint summary
%     .right_summary             - right bracket endpoint summary

    % FIND RIGHT BOUNDARY
    info = [];

    feas = [summary.feasible];
    p3   = [summary.par3];

    if ~any(feas)
        return;
    end

    last_feasible = find(feas, 1, 'last');

    if last_feasible == numel(p3)
        info.bracket = [p3(end), p3(end)];
        info.left_summary  = summary(end);
        info.right_summary = summary(end);
    else
        info.bracket = [p3(last_feasible), p3(last_feasible+1)];
        info.left_summary  = summary(last_feasible);
        info.right_summary = summary(last_feasible+1);
    end
end


%% DOMINANT / SECONDARY FREQUENCIES
function [f_dom, f_sec] = dominant_frequencies_EP(eigvals)
% DOMINANT_FREQUENCIES_EP Extract dominant and secondary oscillation frequencies.
%
% INPUT:
%   eigvals       - matrix of eigenvalues; one column per continuation point
%
% OUTPUT:
%   f_dom / f_sec - dom. (sec.) frequency at each point [Hz]

    n_cols = size(eigvals, 2);
    f_dom = nan(1, n_cols);
    f_sec = nan(1, n_cols);

    for j = 1:n_cols
        lam = eigvals(:, j);

        % sort by closeness of real part to imaginary axis
        [~, idx] = sort(abs(real(lam)));

        idx_dom = idx(1);

        % secondary index: choose next eigenvalue with nonzero imaginary part if possible
        idx_sec = NaN;
        for k = 3:length(idx)
            if abs(imag(lam(idx(k)))) > 0
                idx_sec = idx(k);
                break
            end
        end

        omega_dom = abs(imag(lam(idx_dom)));
        f_dom(j) = omega_dom / (2*pi) * 1000;

        if ~isnan(idx_sec)
            omega_sec = abs(imag(lam(idx_sec)));
            f_sec(j) = omega_sec / (2*pi) * 1000;
        end
    end
end


%% FIND HOPF POINTS
function idx = find_hopf_point(s)
% FIND_HOPF_POINT Return indices of Hopf bifurcation points from MatCont output.
%
% INPUT:
%   s           - MatCont singularity structure
%
% OUTPUT:
%   idx         - indices of entries labelled 'H'

    if isempty(s)
        idx = [];
        return
    end

    mask = strcmp(strtrim({s.label}), 'H');
    idx  = [s(mask).index];
end


%% FREQUENCY STATISTICS
function [f_min, f_max, f_mean, f_std] = freq_stats(f)
% FREQ_STATS Compute basic statistics of a frequency vector.
%
% INPUT:
%   f               - vector of frequencies [Hz]
%
% OUTPUT:
%   f_min / f_max   - minimum / maximum valid frequency [Hz]
%   f_mean / f_std  - mean / std. dev. of valid frequency [Hz]

    f = f(:);
    f = f(~isnan(f));

    if isempty(f)
        f_min  = NaN;
        f_max  = NaN;
        f_mean = NaN;
        f_std  = NaN;
    else
        f_min  = min(f);
        f_max  = max(f);
        f_mean = mean(f);
        f_std  = std(f);
    end
end


%% DRAW FEASIBLE BAND PATCH
function draw_feasible_band_patch(ax, summary, cfg, left_bracket, right_bracket)
% DRAW_FEASIBLE_BAND_PATCH Draw feasible par3 x target-band rectangle(s).

    y1 = cfg.target_band(1);
    y2 = cfg.target_band(2);

    % first try refined interval (if available)
    if ~isempty(left_bracket) && ~isempty(right_bracket)
        x1 = mean(left_bracket);
        x2 = mean(right_bracket);

        patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], ...
              0.8*[1 1 1], ...
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.3);

    else
        % fallback to coarse feasible intervals (possibly multiple blocks)
        feas = [summary.feasible];
        p3   = [summary.par3];
        idx = find(feas);

        if ~isempty(idx)
            breaks = [1, find(diff(idx) > 1) + 1, length(idx)+1];

            for b = 1:length(breaks)-1
                block = idx(breaks(b):breaks(b+1)-1);

                x1 = p3(block(1));
                x2 = p3(block(end));

                patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], ...
                      0.8*[1 1 1], ...
                      'EdgeColor', 'none', ...
                      'FaceAlpha', 0.3);
            end
        end
    end
end


%% FMT CELL
function out = fmt_interval(a, b)
% FMT_INTERVAL Format interval [a, b] for table output.
%
% INPUT:
%   a / b           - lower / upper bound
%
% OUTPUT:
%   out             - formatted character vector

    if isempty(a) || isempty(b) || any(isnan([a b]))
        out = '';
    else
        out = sprintf(' %5.1f %6.1f', a, b);
    end
end

function out = fmt_num1(x)
% FMT_NUM1 Format one scalar number for table output.
%
% INPUT:
%   x               - scalar value
%
% OUTPUT:
%   out             - formatted character vector

    if isempty(x) || any(isnan(x))
        out = '';
    else
        out = sprintf(' %5.1f', x);
    end
end

function out = fmt_int(x)
% FMT_INT Format one integer for table output.
%
% INPUT:
%   x               - integer-like scalar
%
% OUTPUT:
%   out             - formatted character vector

    if isempty(x) || any(isnan(x))
        out = '';
    else
        out = sprintf(' %d', round(x));
    end
end

function out = fmt_logical(x)
% FMT_LOGICAL Format one logical for table output.
%
% INPUT:
%   x               - logical / logical-like scalar
%
% OUTPUT:
%   out             - formatted character vector

    if isempty(x) || any(isnan(x))
        out = '';
    else
        out = sprintf('    %d', logical(x));
    end
end

function out = fmt_string(x)
% FMT_STRING Convert string-like input to char for table output.
%
% INPUT:
%   x               - string / char / empty value
%
% OUTPUT:
%   out             - character vector

    if isempty(x)
        out = '';
    else
        out = char(string(x));
    end
end