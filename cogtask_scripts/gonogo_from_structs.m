function metrics_by_period = gonogo_from_structs
% gonogo_from_structs.m
% 
% Reads ALL .mat files under <root_dir>\Data (recursively).
% Infers period from filename token "ses-#" and aggregates metrics
% by (participant × period). Outputs figures and CSV to <root_dir>\Results.
%
% Each .mat is expected to contain:
%   results (1xN struct) with fields: stimulus, response(0/1), RT, ISI, isGo(0/1), correct(0/1)
%   subject (char/string) [optional; otherwise parsed from filename]
%   session (ignored; rely on "ses-#" in filename)
%
% Metrics per participant × period:
%   H, FA (log-linear/Hautus correction), zH, zFA, d', criterion c,
%   mean RT (correct Go only; weighted across files).

%% PATH CONFIGURATION
% Set USE_HARDCODED_PATHS=true for reproducible batch runs.
% Leave it false to choose the project root interactively when the function runs.
USE_HARDCODED_PATHS = false;

% Optional hardcoded paths. Edit these for your machine if needed.
HARD_CODED_ROOT_DIR = '';   % project folder containing Data/ and receiving Results/
HARD_CODED_DATA_DIR = '';   % optional; leave blank to use fullfile(root_dir, 'Data')
HARD_CODED_OUT_DIR = '';    % optional; leave blank to use fullfile(root_dir, 'Results')

if USE_HARDCODED_PATHS
    root_dir = char(HARD_CODED_ROOT_DIR);
    if isempty(root_dir) && isempty(char(HARD_CODED_DATA_DIR))
        error('Set HARD_CODED_ROOT_DIR or HARD_CODED_DATA_DIR, or use interactive folder selection.');
    end
else
    root_dir = uigetdir(pwd, 'Select Go/No-Go project root folder containing Data/');
    if isequal(root_dir, 0)
        disp('Folder selection cancelled.');
        return;
    end
end

if ~isempty(char(HARD_CODED_DATA_DIR))
    data_dir = char(HARD_CODED_DATA_DIR);
else
    data_dir = fullfile(root_dir, 'Data');
end

if ~isempty(char(HARD_CODED_OUT_DIR))
    out_dir = char(HARD_CODED_OUT_DIR);
else
    out_dir = fullfile(root_dir, 'Results');
end
if ~exist(data_dir,'dir'); error('Data folder not found: %s', data_dir); end
if ~exist(out_dir,'dir');  mkdir(out_dir); end

task_out_dir = fullfile(out_dir, 'GoNoGo');
csv_dir      = fullfile(task_out_dir, 'CSV');
stats_dir    = fullfile(task_out_dir, 'Stats');
fig_dir      = fullfile(task_out_dir, 'Figures');
part_fig_dir = fullfile(fig_dir, 'Participants');
group_fig_dir= fullfile(fig_dir, 'Group');

dirs_to_make = {task_out_dir, csv_dir, stats_dir, fig_dir, part_fig_dir, group_fig_dir};
for id = 1:numel(dirs_to_make)
    if ~exist(dirs_to_make{id}, 'dir')
        mkdir(dirs_to_make{id});
    end
end

csv_out = fullfile(csv_dir, 'gonogo_metrics_by_period.csv');
csv_out_filelevel = fullfile(csv_dir, 'gonogo_file_level_counts.csv');
csv_out_inventory = fullfile(csv_dir, 'gonogo_file_inventory.csv');

%% EDITABLE: ses-# → period mapping
% Adjust these if your numbering differs.
ses_to_period = containers.Map( ...
    {1, 2, 3, 4, 5}, ...
    {'Baseline TMSEEG', 'Dosing Session', '1-week follow-up','2-week follow-up','1-month follow-up'} ...
);

% Canonical plotting order (clinical timeline)
PERIOD_ORDER = [ ...
    "Baseline fMRI", ...
    "Baseline TMSEEG", ...
    "Dosing Session", ...
    "1-week follow-up", ...
    "2-week follow-up", ...
    "1-month follow-up" ...
];

%% EXPECTED FIELDS INSIDE `results`
fld.results  = 'results';
fld.subject  = 'subject';   % optional
fld.isGo     = 'isGo';
fld.response = 'response';
fld.correct  = 'correct';
fld.RT       = 'RT';

%% FIND FILES
% Search all .mat files and identify Go/No-Go files by their contents.
% This avoids missing files if the filename does not contain CCPT.
files = dir(fullfile(data_dir, '**', '*.mat'));

if isempty(files)
    error('No .mat files found under: %s', data_dir);
end

%% PER-FILE COUNTS
rows = [];
inventory_rows = [];
skipped_map = 0;
skipped_nores = 0;
skipped_not_task = 0;

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);

    inv.file = string(fpath);
    inv.filename = string(files(k).name);
    inv.detected_task = false;
    inv.mapped = false;
    inv.period = "";
    inv.participant = "";
    inv.session_num = NaN;
    inv.reason = "";

    [pid, sesnum] = parseFilename(files(k).name);
    pid_from_path = inferParticipantFromPath(fpath);
    if (pid == "unknown" || strlength(string(pid)) > 4) && pid_from_path ~= ""
        pid = pid_from_path;
    end

    period_from_path = inferClinicalPeriodFromPath(fpath);
    if period_from_path ~= ""
        period = period_from_path;
    elseif ~isnan(sesnum) && isKey(ses_to_period, sesnum)
        period = string(ses_to_period(sesnum));
    else
        period = "";
    end

    try
        S = load(fpath);
    catch ME
        inv.reason = "Could not load MAT file: " + string(ME.message);
        inventory_rows = [inventory_rows; inv]; %#ok<AGROW>
        continue;
    end

    if isfield(S, fld.results)
        R = S.(fld.results);
    else
        R = pickFirstStructArray(S);
    end

    if ~isstruct(R) || isempty(R)
        skipped_nores = skipped_nores + 1;
        inv.reason = "No results struct";
        inventory_rows = [inventory_rows; inv]; %#ok<AGROW>
        continue;
    end

    Rt = struct2table(R);
    required_fields = {fld.isGo, fld.response, fld.correct, fld.RT};
    if ~all(ismember(required_fields, Rt.Properties.VariableNames))
        skipped_not_task = skipped_not_task + 1;
        inv.reason = "Results struct does not contain Go/No-Go fields";
        inventory_rows = [inventory_rows; inv]; %#ok<AGROW>
        continue;
    end
    inv.detected_task = true;

    if period == ""
        skipped_map = skipped_map + 1;
        inv.reason = "Could not infer period from folder path or filename";
        inventory_rows = [inventory_rows; inv]; %#ok<AGROW>
        fprintf('Skipping Go/No-Go file (cannot infer period): %s\n', fpath);
        continue;
    end
    inv.mapped = true;
    inv.period = period;
    inv.participant = string(pid);
    inv.session_num = double(sesnum);
    inv.reason = "Processed";

    if isfield(S, fld.subject) && ~isempty(S.(fld.subject))
        pid = string(S.(fld.subject));
        inv.participant = string(pid);
    end

    isGo      = logical(Rt.(fld.isGo));
    responded = logical(Rt.(fld.response));
    correct   = logical(Rt.(fld.correct));
    RT        = Rt.(fld.RT);

    nGo    = sum(isGo == 1);
    nNoGo  = sum(isGo == 0);
    hits   = sum(isGo & correct);
    fas    = sum(~isGo & responded);
    rt_sum = sum(RT(isGo & correct), 'omitnan');

    rows = [rows; makeFileRow(pid, period, nGo, nNoGo, hits, fas, rt_sum, fpath, sesnum)]; %#ok<AGROW>
    inventory_rows = [inventory_rows; inv]; %#ok<AGROW>
end

if ~isempty(inventory_rows)
    writetable(struct2table(inventory_rows), csv_out_inventory);
end

file_metrics = struct2table(rows);
if isempty(file_metrics)
    error('No valid files with mappable ses-# were processed. Skipped (map): %d, (no results): %d', skipped_map, skipped_nores);
end

%% ---------- AGGREGATE BY (participant × period) ----------
grp = findgroups(file_metrics.participant, file_metrics.period);
agg.nGo    = splitapply(@nansum, file_metrics.n_go,    grp);
agg.nNoGo  = splitapply(@nansum, file_metrics.n_nogo,  grp);
agg.hits   = splitapply(@nansum, file_metrics.hits,    grp);
agg.fas    = splitapply(@nansum, file_metrics.fas,     grp);
agg.rt_sum = splitapply(@nansum, file_metrics.rt_sum,  grp);

[pIDs, periods] = splitapply(@firstPair, file_metrics.participant, file_metrics.period, grp); %#ok<ASGLU>

metrics_by_period = table(pIDs, periods, ...
    agg.nGo, agg.nNoGo, agg.hits, agg.fas, agg.rt_sum, ...
    'VariableNames', {'participant','period','n_go','n_nogo','hits','fas','rt_sum'});

% SDT rates (Hautus), z, d', c, mean RT (weighted)
metrics_by_period.H   = arrayfun(@(h,n) safeRate(h, n),  metrics_by_period.hits,  metrics_by_period.n_go);
metrics_by_period.FA  = arrayfun(@(f,n) safeRate(f, n),  metrics_by_period.fas,   metrics_by_period.n_nogo);
metrics_by_period.zH  = safeNormInv(metrics_by_period.H);
metrics_by_period.zFA = safeNormInv(metrics_by_period.FA);
metrics_by_period.dprime    = metrics_by_period.zH - metrics_by_period.zFA;
metrics_by_period.criterion = -0.5 * (metrics_by_period.zH + metrics_by_period.zFA);
metrics_by_period.mean_rt   = metrics_by_period.rt_sum ./ metrics_by_period.hits;
metrics_by_period.mean_rt(metrics_by_period.hits==0) = NaN;

% ROBUST ORDERING (fix for your error)
% Make categorical, then reorder only to the intersection with PERIOD_ORDER.
metrics_by_period.period = categorical(string(metrics_by_period.period)); % ensure valid, nonblank
present = categories(metrics_by_period.period);
desired = cellstr(PERIOD_ORDER(:)');      % cellstr for older MATLAB compatibility
keep_order = intersect(desired, present, 'stable');
metrics_by_period.period = reordercats(metrics_by_period.period, keep_order);
metrics_by_period = sortrows(metrics_by_period, {'participant','period'});

% Save CSV
try
    writetable(file_metrics, csv_out_filelevel);
    writetable(metrics_by_period, csv_out);
catch
    warning('Could not write one or more Go/No-Go CSV files. Check whether a file is open.');
end
% Save group summary and basic longitudinal statistics
numeric_vars = {'H','FA','dprime','criterion','mean_rt'};
group_summary = makeGroupSummary(metrics_by_period, numeric_vars);
writetable(group_summary, fullfile(stats_dir, 'gonogo_group_summary_by_period.csv'));

baseline_stats = baselinePairedStats(metrics_by_period, numeric_vars, "Baseline TMSEEG");
writetable(baseline_stats, fullfile(stats_dir, 'gonogo_paired_ttests_vs_baseline.csv'));
fprintf('Group summary saved: %s\n', fullfile(stats_dir, 'gonogo_group_summary_by_period.csv'));
fprintf('Paired stats saved: %s\n', fullfile(stats_dir, 'gonogo_paired_ttests_vs_baseline.csv'));


%% PLOTTING
participants = unique(metrics_by_period.participant, 'stable');

% Per-participant time series
for ip = 1:numel(participants)
    sub = metrics_by_period(metrics_by_period.participant==participants(ip), :);
    sub = sortrows(sub, 'period');

    % Handle participants missing some periods gracefully
    x = 1:height(sub);
    xlbls = string(sub.period);

    f = figure('Color','k','Position',[220 80 950 900]);
    tl = tiledlayout(5,1,'TileSpacing','compact','Padding','compact');
    title(tl, sprintf('Participant %s — Go/No-Go over periods', participants(ip)));

    nexttile; plot(x, sub.H, '-o','LineWidth',1.5); grid on; ylabel('Hit Rate'); ylim([0 1]);
    set(gca,'XTick',x,'XTickLabel',xlbls,'XTickLabelRotation',45);

    nexttile; plot(x, sub.FA, '-o','LineWidth',1.5); grid on; ylabel('FA Rate'); ylim([0 1]);
    set(gca,'XTick',x,'XTickLabel',xlbls,'XTickLabelRotation',45);

    nexttile; plot(x, sub.dprime, '-o','LineWidth',1.5); grid on; ylabel("d'");
    set(gca,'XTick',x,'XTickLabel',xlbls,'XTickLabelRotation',45);

    nexttile; plot(x, sub.criterion, '-o','LineWidth',1.5); grid on; ylabel('c');
    set(gca,'XTick',x,'XTickLabel',xlbls,'XTickLabelRotation',45);

    nexttile; plot(x, sub.mean_rt, '-o','LineWidth',1.5); grid on; ylabel('Mean RT');
    set(gca,'XTick',x,'XTickLabel',xlbls,'XTickLabelRotation',45); xlabel('Period');

    exportgraphics( ...
        f, ...
        fullfile(part_fig_dir, sprintf('participant_%s_timeseries.png', participants(ip))), ...
        'Resolution', 300);
    close(f);
end

% Group trend plots (use only periods actually present)
plotGroupTrend(metrics_by_period, 'H',         'Hit Rate',      keep_order, group_fig_dir);
plotGroupTrend(metrics_by_period, 'FA',        'False Alarm',   keep_order, group_fig_dir);
plotGroupTrend(metrics_by_period, 'dprime',    "d'",            keep_order, group_fig_dir);
plotGroupTrend(metrics_by_period, 'criterion', 'Criterion (c)', keep_order, group_fig_dir);
plotGroupTrend(metrics_by_period, 'mean_rt',   'Mean RT',       keep_order, group_fig_dir);

fprintf('Done.\nResults (CSV & figures) -> %s\n', out_dir);
if skipped_map>0 || skipped_nores>0
    fprintf('Note: skipped %d file(s) (unmapped ses-#) and %d file(s) (no results struct).\n', skipped_map, skipped_nores);
end
end

%% HELPERS
function [pid, sesnum] = parseFilename(name)
% Parse "Subject_0005_ses-5_CCPT_11.37h_10102025.mat"
    pid = "unknown";
    sesnum = NaN;
    try
        % subject id (digits after "Subject_" or first run of 3+ digits)
        m = regexp(name,'Subject[_-]?(\d+)|(\d{3,})','tokens','once');
        if ~isempty(m)
            g = m(~cellfun('isempty',m));
            if ~isempty(g); pid = string(g{1}); end
        end
        % ses-#
        s = regexp(name,'ses[-_]?(\d+)','tokens','once');
        if ~isempty(s); sesnum = str2double(s{1}); end
    catch
        % leave defaults
    end
end

function mustHave(T, name, fname)
    if ~ismember(name, T.Properties.VariableNames)
        error('Missing field "%s" in %s', name, fname);
    end
end

function row = makeFileRow(pid, period, nGo, nNoGo, hits, fas, rt_sum, fpath, sesnum)
    row.participant = string(pid);
    row.period      = string(period);
    row.sesnum      = double(sesnum);
    row.n_go        = double(nGo);
    row.n_nogo      = double(nNoGo);
    row.hits        = double(hits);
    row.fas         = double(fas);
    row.rt_sum      = double(rt_sum);
    row.file        = string(fpath);
end

function z = safeNormInv(p)
    try
        z = norminv(p,0,1);
    catch
        z = -sqrt(2) * erfcinv(2*p);
    end
end

function r = safeRate(k, n)
% Log-linear (Hautus) adjusted rate: (k + 0.5) / (n + 1)
    if n==0
        r = NaN;
    else
        r = (k + 0.5) / (n + 1);
        epsv = 1e-12;
        r = min(max(r, epsv), 1 - epsv);
    end
end

function [firstA, firstB] = firstPair(A, B)
    firstA = A(1);
    firstB = B(1);
end

function R = pickFirstStructArray(S)
% Find first non-empty struct array in a loaded .mat
    R = [];
    fns = fieldnames(S);
    for i = 1:numel(fns)
        val = S.(fns{i});
        if isstruct(val) && ~isempty(val)
            R = val; return;
        end
    end
end

function plotGroupTrend(M, varname, ylab, ORDER_CELLSTR, out_dir)
    % ORDER_CELLSTR is a cellstr of desired period order present in the data.
    if isempty(ORDER_CELLSTR)
        % fall back to natural order present in M
        ORDER_CELLSTR = categories(categorical(string(M.period)));
    end

    % Order by provided sequence
    P = categorical(string(M.period));
    P = reordercats(P, ORDER_CELLSTR);
    M.period = P;
    M = sortrows(M, {'period','participant'});

    sessions = categories(M.period);
    parts    = unique(M.participant,'stable');

    % Build Y(period × participant)
    Y = NaN(numel(sessions), numel(parts));
    for ip = 1:numel(parts)
        sub = M(M.participant==parts(ip), {'period', varname});
        [~,loc] = ismember(string(sub.period), sessions);
        Y(loc, ip) = sub.(varname);
    end

    mu  = mean(Y, 2, 'omitnan');
    sem = std(Y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Y),2));

    x = 1:numel(sessions);
    f = figure('Color','k','Position',[220 100 1000 450]); hold on; box on; grid on;

    for ip = 1:numel(parts)
        plot(x, Y(:,ip), '-', 'LineWidth', 0.8, 'HandleVisibility','off');
    end

    if numel(x) >= 2
        px = [x, fliplr(x)];
        py = [mu-sem; flipud(mu+sem)]';
        patch(px, py, [0.7 0.7 0.7], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility','off');
    end
    plot(x, mu, '-o', 'LineWidth', 2.5, 'DisplayName','Group mean');

    set(gca,'XTick',x,'XTickLabel',sessions,'XTickLabelRotation',45);
    xlabel('Period'); ylabel(ylab); title(sprintf('Group trend — %s', ylab));
    legend('Location','best');

    exportgraphics( ...
        f, ...
        fullfile(out_dir, sprintf('group_trend_%s.png', varname)), ...
        'Resolution', 300);
    close(f);
end


function group_summary = makeGroupSummary(T, numeric_vars)
    rows = [];
    periods = categories(categorical(string(T.period)));

    for ip = 1:numel(periods)
        period = string(periods{ip});
        maskP = string(T.period) == period;

        for iv = 1:numel(numeric_vars)
            varname = numeric_vars{iv};
            vals = T.(varname)(maskP);
            vals = vals(~isnan(vals));

            row.period = period;
            row.metric = string(varname);
            row.n = double(numel(vals));
            row.mean = mean(vals, 'omitnan');
            row.sd = std(vals, 0, 'omitnan');
            row.sem = row.sd ./ sqrt(row.n);
            row.median = median(vals, 'omitnan');
            row.min = min(vals);
            row.max = max(vals);

            rows = [rows; row]; %#ok<AGROW>
        end
    end

    if isempty(rows)
        group_summary = table();
    else
        group_summary = struct2table(rows);
    end
end

function stats_tbl = baselinePairedStats(T, numeric_vars, baseline_label)
    rows = [];

    periods = categories(categorical(string(T.period)));
    periods = string(periods);
    compare_periods = periods(periods ~= baseline_label);

    participants = unique(string(T.participant), 'stable');

    for iv = 1:numel(numeric_vars)
        varname = numeric_vars{iv};

        for iper = 1:numel(compare_periods)
            comp_label = compare_periods(iper);

            base_vals = [];
            comp_vals = [];

            for isub = 1:numel(participants)
                pid = participants(isub);

                bmask = string(T.participant) == pid & string(T.period) == baseline_label;
                cmask = string(T.participant) == pid & string(T.period) == comp_label;

                if any(bmask) && any(cmask)
                    b = mean(T.(varname)(bmask), 'omitnan');
                    c = mean(T.(varname)(cmask), 'omitnan');

                    if ~isnan(b) && ~isnan(c)
                        base_vals(end+1,1) = b; %#ok<AGROW>
                        comp_vals(end+1,1) = c; %#ok<AGROW>
                    end
                end
            end

            row.metric = string(varname);
            row.baseline_period = string(baseline_label);
            row.comparison_period = string(comp_label);
            row.n_pairs = double(numel(base_vals));
            row.mean_baseline = mean(base_vals, 'omitnan');
            row.mean_comparison = mean(comp_vals, 'omitnan');
            row.mean_difference = mean(comp_vals - base_vals, 'omitnan');

            if numel(base_vals) >= 2
                [~,p,ci,stats] = ttest(comp_vals, base_vals);
                row.t_stat = stats.tstat;
                row.df = stats.df;
                row.p_value = p;
                row.ci_low = ci(1);
                row.ci_high = ci(2);
            else
                row.t_stat = NaN;
                row.df = NaN;
                row.p_value = NaN;
                row.ci_low = NaN;
                row.ci_high = NaN;
            end

            rows = [rows; row]; %#ok<AGROW>
        end
    end

    if isempty(rows)
        stats_tbl = table();
    else
        stats_tbl = struct2table(rows);
    end
end


function period = inferClinicalPeriodFromPath(fpath)
    s = lower(string(fpath));
    period = "";

    if contains(s, "baseline tmseeg") || contains(s, "baseline_tms") || contains(s, "baseline tms")
        period = "Baseline TMSEEG";
    elseif contains(s, "baseline fmri") || contains(s, "baseline_fmri")
        period = "Baseline fMRI";
    elseif contains(s, "dosing session") || contains(s, "dosing")
        period = "Dosing Session";
    elseif contains(s, "1-week follow-up") || contains(s, "1 week follow-up") || contains(s, "1-week") || contains(s, "1_week")
        period = "1-week follow-up";
    elseif contains(s, "2-week follow-up") || contains(s, "2 week follow-up") || contains(s, "2-week") || contains(s, "2_week")
        period = "2-week follow-up";
    elseif contains(s, "1-month follow-up") || contains(s, "1 month follow-up") || contains(s, "1-month") || contains(s, "1_month")
        period = "1-month follow-up";
    end
end

function pid = inferParticipantFromPath(fpath)
    pid = "";
    parts = regexp(char(fpath), '[\\/]', 'split');
    for i = 1:numel(parts)
        tok = regexp(parts{i}, '^0*(\d{1,4})$', 'tokens', 'once');
        if ~isempty(tok)
            pid = string(tok{1});
            return;
        end
    end
end
