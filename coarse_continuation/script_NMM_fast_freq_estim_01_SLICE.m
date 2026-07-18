clear
close all
warning('off')

% MatCont procedures initialization
init

% LaTeX
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');


%% USER SETTINGS

% --- par1 & par2 ---
par1_range = [0 20]; par1_idx = 4; par1_label = '$M_{D1}$';                 % Md
par2_range = [0 6];  par2_idx = 9; par2_label = '$\Delta_{v_\theta}$';      % Delta

% --- alternative choices for par1 ---
% par1_range = [1 200];   par1_idx = 1;  par1_label = '$C$';                  % C
% par1_range = [0 70];    par1_idx = 2;  par1_label = '$I_{ext}$';            % Iext
% par1_range = [0 120];   par1_idx = 3;  par1_label = '$u_{jump}$';           % ujump
% par1_range = [0 1];     par1_idx = 5;  par1_label = '$B_{M_{D1}}$';         % B
% par1_range = [1 3];     par1_idx = 6;  par1_label = '$\gamma$';             % gamma
% par1_range = [-80 -50]; par1_idx = 7;  par1_label = '$v_r$';                % Vr
% par1_range = [-55 -35]; par1_idx = 8;  par1_label = '$\bar{v}_{\theta}$';   % Vt
% par1_range = [1 10];    par1_idx = 11; par1_label = '$\tau_{s_a}$';         % Ta

% --- names / indices of parameters in pars vector ---
cfg = struct();
cfg.par1_idx = par1_idx;
cfg.par2_idx = par2_idx;

cfg.par1_range = par1_range;
cfg.par2_range = par2_range;

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

slice = compute_hopf_slice(cfg);


%% SUMMARY
fprintf('\n=== SINGLE SLICE COMPUTATION SUMMARY ===\n');

% --- table header ---
header1 = sprintf(['%-15s # %-32s # %-32s # \n'], ...
    '', '       DOMINANT FREQUENCY', '       SECONDARY FREQUENCY');
header2 = sprintf([ ...
    '%-15s # %-14s | %-6s | %-6s # %-14s | %-6s | %-6s # %-8s | %-20s\n'], ...
    'item', '[f_min, f_max]', 'f_mean', 'f_std', '[f_min, f_max]', ...
    'f_mean', 'f_std', 'n_points', 'status');

fprintf('%s\n', repelem('-', length(header2)-1));
fprintf('%s', header1);
fprintf('%s', header2);
fprintf('%s\n', repelem('=', length(header2)-1));

% --- slice row ---
fprintf([ ...
    '%-15s # %-14s | %-6s | %-5s # %-14s | %-6s | %-5s # %-8s | %-20s\n'], ...
    'slice', ...
    fmt_interval(slice.f_min, slice.f_max), ...
    fmt_num1(slice.f_mean), ...
    fmt_num1(slice.f_std), ...
    fmt_interval(slice.f2_min, slice.f2_max), ...
    fmt_num1(slice.f2_mean), ...
    fmt_num1(slice.f2_std), ...
    sprintf(' %d / %d', slice.num_points, numel(slice.branches)), ...
    fmt_string(slice.status));

% --- branch rows ---
if ~isempty(slice.branches)
    for j = 1:numel(slice.branches)
        bj = slice.branches(j);

        if isfield(bj,'f') && ~isempty(bj.f)
            branch_n_points = numel(bj.f);
        else
            branch_n_points = [];
        end

        fprintf([ ...
            '%-15s # %-14s | %-6s | %-6s # %-14s | %-6s | %-6s # %-8s | %-20s\n'], ...
            sprintf('     branch %d', bj.branch_id), ...
            fmt_interval(bj.f_min, bj.f_max), ...
            fmt_num1(bj.f_mean), ...
            fmt_num1(bj.f_std), ...
            fmt_interval(bj.f2_min, bj.f2_max), ...
            fmt_num1(bj.f2_mean), ...
            fmt_num1(bj.f2_std), ...
            fmt_int(branch_n_points), ...
            '');
    end
end

fprintf('%s\n\n', repelem('-', length(header2)-1));

toc;


%% VISUALIZATION
titles_on = true;

% SCATTER (slice.points: [par1, par2, f_dom, f_sec, branch_id]) ----------
if ~isempty(slice.points)

    fig = figure;
    set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[11,4.5])

    tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    % --- dominant frequency plot ---
    ax1 = nexttile;
    scatter(ax1, slice.points(:,1), slice.points(:,2), 20, slice.points(:,3), 'filled');

    xlabel(ax1, par1_label);
    ylabel(ax1, par2_label);
    if titles_on
        title(ax1, 'DOMINANT FREQUENCY');
    end

    box(ax1, 'on');
    grid(ax1, 'on');
    xlim(ax1, cfg.par1_range);
    ylim(ax1, cfg.par2_range);
    xtickangle(ax1, 0);
    ytickangle(ax1, 0);

    cb1 = colorbar(ax1);
    cb1.Label.String = 'Frequency [Hz]';
    cb1.TickLabelInterpreter = 'latex';
    cb1.Label.Interpreter = 'latex';

    % --- secondary frequency plot ---
    ax2 = nexttile;
    scatter(ax2, slice.points(:,1), slice.points(:,2), 20, slice.points(:,4), 'filled');

    xlabel(ax2, par1_label);
    ylabel(ax2, par2_label);
    if titles_on
        title(ax2, 'SECONDARY FREQUENCY');
    end

    box(ax2, 'on');
    grid(ax2, 'on');
    xlim(ax2, cfg.par1_range);
    ylim(ax2, cfg.par2_range);
    xtickangle(ax2, 0);
    ytickangle(ax2, 0);

    cb2 = colorbar(ax2);
    cb2.Label.String = 'Frequency [Hz]';
    cb2.TickLabelInterpreter = 'latex';
    cb2.Label.Interpreter = 'latex';
end

warning('on')


%% COMPUTE HOPF SLICE
function slice = compute_hopf_slice(cfg)
% COMPUTE_HOPF_SLICE Compute one Hopf slice in the (par1, par2) plane.
%
% INPUT:
%   cfg                         - structure with slice settings:
%     .par1_idx / .par2_idx     - index of parameter par1 (par2) in cfg.pars
%     .par1_range / .par2_range - allowed interval [min max] for par1 (par2)
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
%     .f_min / .f2_min          - minimum dominant (secondary) frequency [Hz]
%     .f_max / .f2_max          - maximum dom. (sec.) frequency [Hz]
%     .f_mean / .f2_mean        - mean dom. (sec.) frequency [Hz]
%     .f_std / .f2_std          - std. dev. of dom. (sec.) frequency [Hz]
%     .num_points               - number of valid Hopf points in the box
%     .status                   - computation status string
%                                 ('init', 'ok', 'no_hopf_on_EP',
%                                  'no_valid_points_in_box', 'failed: ...')
%     .points                   - valid points, columns [par1, par2, f_dom, f_sec, branch_id]
%     .branches                 - structure array of valid Hopf branches:
%       .branch_id              - integer identifier
%       .par1 / .par2           - branch values of par1 (par2)
%       .f / .f2                - dom. (sec.) frequencies [Hz]
%       .f_min / .f2_min        - minimum dom. (sec.) frequency [Hz]
%       .f_max / .f2_max        - maximum dom. (sec.) frequency [Hz]
%       .f_mean / .f2_mean      - mean dom. (sec.) frequency [Hz]
%       .f_std / .f2_std        - std. dev. of dom. (sec.) frequency [Hz]

    slice            = struct();
    slice.f_min      = NaN;
    slice.f_max      = NaN;
    slice.f_mean     = NaN;
    slice.f_std      = NaN;
    slice.f2_min     = NaN;
    slice.f2_max     = NaN;
    slice.f2_mean    = NaN;
    slice.f2_std     = NaN;
    slice.num_points = 0;
    slice.status     = 'init';
    slice.points     = [];   % [par1, par2, f_dom, f_sec, branch_id]
    slice.branches   = [];

    try
        % SET PARAMETERS
        pars = cfg.pars;

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
            'f2_std', {});

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
        
            valid_branches(end+1) = branch_struct;
            valid_points = [valid_points; pts];
            valid_f2_points = [valid_f2_points; f2_sel(:)];
        
        end

        if isempty(valid_points)
            slice.status = 'no_valid_points_in_box';
            return;
        end

        % FREQUENCY INTERVAL
        [f_min, f_max, f_mean, f_std]     = freq_stats(valid_points(:,3));
        [f2_min, f2_max, f2_mean, f2_std] = freq_stats(valid_f2_points);

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
        slice.status     = 'ok';

    catch ME
        slice.status = ['failed: ' ME.identifier];
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
