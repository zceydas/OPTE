function metrics = nback_from_structs
% nback_from_structs.m
% -------------------------------------------------------------
% Reads ALL .mat files under <root_dir>\Data (recursively).
% Each .mat contains a numeric matrix (unnamed columns).
%   Col  6 : N-back level (e.g., 0,1,2,3)
%   Col 10 : correctness code (1=Hit, 2=False Alarm, 3=Miss, 4=Correct Rejection)
%   Col 11 : Reaction time (RT)
%
% File naming: Subject_XXXX_ses-#_*.mat  (ses # mapped to clinical period)
%
% Aggregates by participant × period × level:
%   counts (Hits,FAs,Misses,CRs),
%   SDT rates with Hautus correction, z-scores, d′,
%   mean RT for Hits and FAs.
%
% Outputs:
%   - CSV: <root_dir>\Results\nback_metrics_by_period_level.csv
%   - Participant plots: one figure per participant, lines per level
%   - Group trend plots: mean±SEM across participants, lines per level
% -------------------------------------------------------------

%% ---------- PATH CONFIGURATION ----------
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
    root_dir = uigetdir(pwd, 'Select N-back project root folder containing Data/');
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
nback_dir = fullfile(out_dir, 'NBack');
csv_dir = fullfile(nback_dir, 'CSV');
stats_dir = fullfile(nback_dir, 'Stats');
fig_dir = fullfile(nback_dir, 'Figures');
participant_fig_dir = fullfile(fig_dir, 'Participant_TimeSeries');
group_fig_dir = fullfile(fig_dir, 'Group_Trends');

if ~exist(data_dir,'dir'); error('Data folder not found: %s', data_dir); end
if ~exist(csv_dir,'dir'); mkdir(csv_dir); end
if ~exist(stats_dir,'dir'); mkdir(stats_dir); end
if ~exist(participant_fig_dir,'dir'); mkdir(participant_fig_dir); end
if ~exist(group_fig_dir,'dir'); mkdir(group_fig_dir); end

csv_out_metrics = fullfile(csv_dir, 'nback_metrics_by_period_level.csv');
csv_out_all = fullfile(csv_dir, 'nback_all_participant_metrics.csv');
csv_out_filelevel = fullfile(csv_dir, 'nback_file_level_counts.csv');
csv_out_filemap = fullfile(csv_dir, 'nback_filemap.csv');
csv_out_group_summary = fullfile(stats_dir, 'nback_group_summary_by_period_level.csv');
csv_out_paired_stats = fullfile(stats_dir, 'nback_paired_ttests_vs_baseline_by_level.csv');

% Subject ID formatting
SUBJECT_PAD_WIDTH = 4;  % e.g., 5 -> "0005"

% Session → period mapping (edit if needed)
ses_to_period = containers.Map( ...
    {1,                   2,               3,                  4,                   5}, ...
    {'Baseline TMSEEG',   'Dosing Session','1-week follow-up', '2-week follow-up',  '1-month follow-up'} ...
);
PERIOD_ORDER = unique(string(values(ses_to_period)),'stable');

%% ---------- FIND FILES ----------
% Search by task folder OR filename. This is more robust than relying only on
% filenames, because some files may sit inside an N-back folder even if their
% filename is not perfectly standardized.
all_mat = dir(fullfile(data_dir, '**', '*.mat'));
task_mask = false(numel(all_mat),1);
for k = 1:numel(all_mat)
    fpath_tmp = lower(fullfile(all_mat(k).folder, all_mat(k).name));
    task_mask(k) = contains(fpath_tmp, 'n-back') || contains(fpath_tmp, 'nback');
end
files = all_mat(task_mask);
files = files(~contains({files.name}, '_all', 'IgnoreCase', true));

if isempty(files)
    error('No N-back .mat files found under: %s', data_dir);
end

%% ---------- PREFLIGHT: file → subject/session/period map ----------
frows = [];
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    [pid_raw, sesnum] = parseFilename(files(k).name);
    pid_fmt = normalizeSubject(pid_raw, SUBJECT_PAD_WIDTH);

    mapped = true; reason = "";
    period = "";
    period_from_path = inferClinicalPeriodFromPath(fpath);
    if period_from_path ~= ""
        period = period_from_path;
    elseif ~isnan(sesnum) && isKey(ses_to_period, sesnum)
        period = string(ses_to_period(sesnum));
    else
        mapped = false;
        if isnan(sesnum)
            reason = "No session number found in filename and no period found in folder path";
        else
            reason = "Session number not in ses_to_period map and no period found in folder path";
        end
    end

    pid_from_path = inferParticipantFromPath(fpath);
    if (pid_raw == "unknown" || strlength(string(pid_raw)) > 4) && pid_from_path ~= ""
        pid_raw = pid_from_path;
        pid_fmt = normalizeSubject(pid_raw, SUBJECT_PAD_WIDTH);
    end

    frows = [frows; struct( ...
        'file', fpath, ...
        'subject_raw', string(pid_raw), ...
        'subject', string(pid_fmt), ...
        'session_num', double(sesnum), ...
        'period', string(period), ...
        'mapped', logical(mapped), ...
        'reason', string(reason)) ...
    ];
end

filemap = struct2table(frows);
writetable(filemap, csv_out_filemap);
fprintf('Preflight map saved: %s\n', csv_out_filemap);

% Summary
n_total  = height(filemap);
n_mapped = sum(filemap.mapped);
n_unmapped = n_total - n_mapped;
fprintf('Found %d files: %d mapped, %d unmapped.\n', n_total, n_mapped, n_unmapped);
if n_unmapped > 0
    um = filemap(~filemap.mapped, {'file','reason'});
    disp(um(1:min(10, height(um)), :));
end

% Keep mapped rows
filemap = filemap(filemap.mapped, :);
if isempty(filemap)
    error('All files were unmapped. Check nback_filemap.csv and adjust parsing or ses_to_period.');
end

%% ---------- EXTRACT PER-FILE STATS (BY LEVEL) ----------
rows = [];
skipped_nodata = 0;  % no usable "results" matrix

for k = 1:height(filemap)
    fpath = filemap.file{k};
    pid   = filemap.subject{k};
    period= filemap.period{k};

    % Load and obtain numeric matrix "M" with >= 11 columns
    S = load(fpath);
    M = coerceResultsToNumeric(S, 11);  % <---- conversion-aware
    if isempty(M)
        skipped_nodata = skipped_nodata + 1;
        fprintf('Skipping (no usable "results" with >=11 cols): %s\n', fpath);
        continue;
    end

    % Column extraction (1-based)
    level = M(:,6);
    corr  = M(:,10);      % 1=Hit, 2=FA, 3=Miss, 4=CR
    RT    = M(:,11);      % ms; -99 means missing

    % Valid rows: known correctness + finite level
    valid = ismember(corr, [1 2 3 4]) & ~isnan(level);
    if ~any(valid); continue; end
    level = level(valid);  corr = corr(valid);  RT = RT(valid);

    % Exclude RT==-99 from RT means (treat as missing)
    RT(RT == -99) = NaN;

    % Per-level counts for this file
    ulev = unique(level(:)');  % row vector
    for L = ulev
        mask = (level == L);

        hits = sum(corr(mask) == 1);
        fas  = sum(corr(mask) == 2);
        miss = sum(corr(mask) == 3);
        cr   = sum(corr(mask) == 4);

        rt_hit_sum = sum(RT(mask & corr==1), 'omitnan');
        rt_fa_sum  = sum(RT(mask & corr==2), 'omitnan');

        rows = [rows; makeFileRow(pid, period, L, hits, fas, miss, cr, rt_hit_sum, rt_fa_sum)];
    end
end

file_tbl = struct2table(rows);
if isempty(file_tbl)
    error('No valid N-back rows were produced. (Unmapped files: %d, No-data files: %d)', n_unmapped, skipped_nodata);
end

%% ---------- AGGREGATE BY (participant × period × level) ----------
grp = findgroups(file_tbl.participant, file_tbl.period, file_tbl.level);
A.hits       = splitapply(@nansum, file_tbl.hits,        grp);
A.fas        = splitapply(@nansum, file_tbl.fas,         grp);
A.miss       = splitapply(@nansum, file_tbl.miss,        grp);
A.cr         = splitapply(@nansum, file_tbl.cr,          grp);
A.rt_hit_sum = splitapply(@nansum, file_tbl.rt_hit_sum,  grp);
A.rt_fa_sum  = splitapply(@nansum, file_tbl.rt_fa_sum,   grp);
[pIDs, periods, levels] = splitapply(@firstTriple, file_tbl.participant, file_tbl.period, file_tbl.level, grp);

metrics = table(pIDs, periods, levels, A.hits, A.fas, A.miss, A.cr, A.rt_hit_sum, A.rt_fa_sum, ...
    'VariableNames', {'participant','period','level','hits','fas','miss','cr','rt_hit_sum','rt_fa_sum'});

% SDT rates (Hautus), z, d′
nSignal = metrics.hits + metrics.miss;
nNoise  = metrics.fas  + metrics.cr;
metrics.hit_rate_adj = (metrics.hits + 0.5) ./ (nSignal + 1);
metrics.fa_rate_adj  = (metrics.fas  + 0.5) ./ (nNoise  + 1);
metrics.zH  = safeNormInv(metrics.hit_rate_adj);
metrics.zFA = safeNormInv(metrics.fa_rate_adj);
metrics.dprime = metrics.zH - metrics.zFA;

% Mean RTs (per condition)
metrics.mean_rt_hit = metrics.rt_hit_sum ./ metrics.hits;
metrics.mean_rt_hit(metrics.hits==0) = NaN;
metrics.mean_rt_fa  = metrics.rt_fa_sum ./ metrics.fas;
metrics.mean_rt_fa(metrics.fas==0) = NaN;

% Order periods robustly
metrics.period = categorical(string(metrics.period));
present = categories(metrics.period);
desired = cellstr(PERIOD_ORDER(:)');
keep_order = intersect(desired, present, 'stable');
metrics.period = reordercats(metrics.period, keep_order);
metrics = sortrows(metrics, {'participant','period','level'});

% Save CSV outputs
writetable(file_tbl, csv_out_filelevel);
writetable(metrics, csv_out_metrics);
writetable(metrics, csv_out_all);
fprintf('Metrics saved: %s\n', csv_out_metrics);

%% ---------- STATISTICS ----------
measure_vars = {'hit_rate_adj','fa_rate_adj','dprime','mean_rt_hit','mean_rt_fa'};
group_summary = makeGroupSummaryLong(metrics, {'period','level'}, measure_vars);
paired_stats = makePairedStatsVsBaseline(metrics, {'level'}, measure_vars, "Baseline TMSEEG");
try
    writetable(group_summary, csv_out_group_summary);
    writetable(paired_stats, csv_out_paired_stats);
catch
    warning('Could not write N-back statistics CSV outputs.');
end

%% ---------- PLOTTING ----------
parts  = unique(metrics.participant,'stable');
levels_all = unique(metrics.level,'stable');

for ip = 1:numel(parts)
    sub = metrics(metrics.participant==parts(ip), :);
    sub = sortrows(sub, {'period','level'});

    participantMultiLine(sub, keep_order, levels_all, 'dprime',       "d'",            participant_fig_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'hit_rate_adj', 'Hit rate (adj)',participant_fig_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'fa_rate_adj',  'FA rate (adj)', participant_fig_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'mean_rt_hit',  'Mean RT (Hits) (ms)',participant_fig_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'mean_rt_fa',   'Mean RT (FAs) (ms)', participant_fig_dir, parts(ip));
end

groupTrendByLevel(metrics, keep_order, levels_all, 'dprime',       "d'",            group_fig_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'hit_rate_adj', 'Hit rate (adj)',group_fig_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'fa_rate_adj',  'FA rate (adj)', group_fig_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'mean_rt_hit',  'Mean RT (Hits)',group_fig_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'mean_rt_fa',   'Mean RT (FAs)', group_fig_dir);

fprintf('Done.\nOrganized N-back results -> %s\n', nback_dir);
end


%% STATS HELPERS
function stats_tbl = makeGroupSummaryLong(T, group_vars, measure_vars)
    rows = [];
    if isempty(T)
        stats_tbl = table();
        return;
    end
    G = unique(T(:, group_vars), 'rows', 'stable');
    for ig = 1:height(G)
        mask = true(height(T),1);
        for gv = 1:numel(group_vars)
            v = group_vars{gv};
            mask = mask & (string(T.(v)) == string(G.(v)(ig)));
        end
        for im = 1:numel(measure_vars)
            mv = measure_vars{im};
            vals = T.(mv)(mask);
            vals = vals(~isnan(vals));
            row = struct();
            for gv = 1:numel(group_vars)
                v = group_vars{gv};
                row.(v) = string(G.(v)(ig));
            end
            row.measure = string(mv);
            row.n = double(numel(vals));
            if isempty(vals)
                row.mean = NaN; row.sd = NaN; row.sem = NaN;
                row.median = NaN; row.min = NaN; row.max = NaN;
            else
                row.mean = mean(vals, 'omitnan');
                row.sd = std(vals, 0, 'omitnan');
                row.sem = row.sd ./ sqrt(row.n);
                row.median = median(vals, 'omitnan');
                row.min = min(vals);
                row.max = max(vals);
            end
            rows = [rows; row]; %#ok<AGROW>
        end
    end
    if isempty(rows), stats_tbl = table(); else, stats_tbl = struct2table(rows); end
end

function paired_tbl = makePairedStatsVsBaseline(T, group_vars, measure_vars, baseline_period)
    rows = [];
    if isempty(T) || ~ismember('period', T.Properties.VariableNames) || ~ismember('participant', T.Properties.VariableNames)
        paired_tbl = table();
        return;
    end
    if isempty(group_vars)
        G = table("all", 'VariableNames', {'dummy_group'});
        group_vars_local = {'dummy_group'};
        T.dummy_group = repmat("all", height(T), 1);
    else
        G = unique(T(:, group_vars), 'rows', 'stable');
        group_vars_local = group_vars;
    end
    all_periods = unique(string(T.period), 'stable');
    comp_periods = all_periods(all_periods ~= string(baseline_period));
    for ig = 1:height(G)
        gmask = true(height(T),1);
        for gv = 1:numel(group_vars_local)
            v = group_vars_local{gv};
            gmask = gmask & (string(T.(v)) == string(G.(v)(ig)));
        end
        for cp = 1:numel(comp_periods)
            comp_period = comp_periods(cp);
            for im = 1:numel(measure_vars)
                mv = measure_vars{im};
                B = T(gmask & string(T.period)==string(baseline_period), {'participant', mv});
                C = T(gmask & string(T.period)==string(comp_period), {'participant', mv});
                B.Properties.VariableNames{2} = 'baseline_value';
                C.Properties.VariableNames{2} = 'comparison_value';
                J = innerjoin(B, C, 'Keys', 'participant');
                diff_vals = J.comparison_value - J.baseline_value;
                diff_vals = diff_vals(~isnan(diff_vals));
                row = struct();
                for gv = 1:numel(group_vars_local)
                    v = group_vars_local{gv};
                    if strcmp(v, 'dummy_group'), continue; end
                    row.(v) = string(G.(v)(ig));
                end
                row.baseline_period = string(baseline_period);
                row.comparison_period = string(comp_period);
                row.measure = string(mv);
                row.n_pairs = double(numel(diff_vals));
                row.mean_diff = mean(diff_vals, 'omitnan');
                row.sd_diff = std(diff_vals, 0, 'omitnan');
                row.sem_diff = row.sd_diff ./ sqrt(row.n_pairs);
                row.tstat = NaN;
                row.df = row.n_pairs - 1;
                row.p = NaN;
                if row.n_pairs >= 2
                    try
                        [~, p, ~, st] = ttest(diff_vals, 0);
                        row.tstat = st.tstat;
                        row.df = st.df;
                        row.p = p;
                    catch
                        row.tstat = NaN;
                        row.p = NaN;
                    end
                end
                rows = [rows; row]; %#ok<AGROW>
            end
        end
    end
    if isempty(rows), paired_tbl = table(); else, paired_tbl = struct2table(rows); end
end

%% ======================= HELPERS ===========================
function pid_fmt = normalizeSubject(pid_raw, width)
    pid_fmt = string(pid_raw);
    num = str2double(pid_raw);
    if ~isnan(num)
        pid_fmt = string(sprintf(['%0',num2str(width),'d'], num));
    end
end

function [pid, sesnum] = parseFilename(name)
% Robust parser for:
%   1) "...ses-<num>..." anywhere
%   2) "^nBack_<pid>_<ses>_" at start (case-insensitive)
% Fallback: first two numeric tokens near the start.
    pid = "unknown";
    sesnum = NaN;
    s = lower(name);

    % Case 1: explicit ses-#
    m = regexp(s, 'ses[-_]?(\d+)', 'tokens', 'once');
    if ~isempty(m)
        sesnum = str2double(m{1});
        sub = regexp(s, 'subject[-_]?(\d+)', 'tokens', 'once');
        if ~isempty(sub)
            pid = string(sub{1});
        else
            pre = regexp(s, '(\d+).*?ses[-_]?(\d+)', 'tokens', 'once');
            if ~isempty(pre), pid = string(pre{1}); end
        end
        return
    end

    % Case 2: nBack_<pid>_<ses>_
    m = regexp(s, '^nback[-_]?(\d+)[-_](\d+)[-_]', 'tokens', 'once');
    if ~isempty(m)
        pid    = string(m{1});
        sesnum = str2double(m{2});
        return
    end

    % Fallback: first two numeric tokens near the start
    toks = regexp(s, '^.*?([0-9]+)[-_]([0-9]+)[-_]', 'tokens', 'once');
    if ~isempty(toks)
        pid    = string(toks{1});
        sesnum = str2double(toks{2});
    end
end

function M = coerceResultsToNumeric(S, minCols)
% Return a numeric matrix with >= minCols columns from S.results.
% Accepts:
%   - numeric matrix already
%   - cell array (string/char/numeric) → str2double
% Returns [] if not found/usable.
    if ~isfield(S, 'results'); M = []; return; end
    R = S.results;

    % Numeric matrix
    if isnumeric(R) && ismatrix(R) && size(R,2) >= minCols
        M = R; return;
    end

    % Cell array → numeric
    if iscell(R) && ismatrix(R) && size(R,2) >= minCols
        % Convert everything via str2double; numeric stays numeric
        M = nan(size(R));
        for c = 1:size(R,2)
            col = R(:,c);
            % Convert each cell to string then str2double
            try
                % faster vectorized conversion where possible
                M(:,c) = cellfun(@(x) str2double(string(x)), col);
            catch
                % fallback per element
                for i = 1:numel(col)
                    M(i,c) = str2double(string(col{i}));
                end
            end
        end
        return;
    end

    % Some files store a cell array of char arrays (e.g., as a char matrix)
    if ischar(R)
        % unlikely in your case; not handled
        M = []; return;
    end

    % Last chance: search any field that looks like results table
    fns = fieldnames(S);
    for i = 1:numel(fns)
        val = S.(fns{i});
        if isnumeric(val) && ismatrix(val) && size(val,2) >= minCols
            M = val; return;
        end
        if iscell(val) && ismatrix(val) && size(val,2) >= minCols
            M = nan(size(val));
            for c = 1:size(val,2)
                M(:,c) = cellfun(@(x) str2double(string(x)), val(:,c));
            end
            return;
        end
    end

    M = [];
end

function z = safeNormInv(p)
    epsv = 1e-12;
    p = max(min(p, 1-epsv), epsv);
    try
        z = norminv(p,0,1);
    catch
        z = -sqrt(2) * erfcinv(2*p);
    end
end

function row = makeFileRow(pid, period, level, hits, fas, miss, cr, rt_hit_sum, rt_fa_sum)
    row.participant = string(pid);
    row.period      = string(period);
    row.level       = double(level);
    row.hits        = double(hits);
    row.fas         = double(fas);
    row.miss        = double(miss);
    row.cr          = double(cr);
    row.rt_hit_sum  = double(rt_hit_sum);
    row.rt_fa_sum   = double(rt_fa_sum);
end

function [a,b,c] = firstTriple(A,B,C)
    a = A(1); b = B(1); c = C(1);
end

function participantMultiLine(sub, keep_order, levels_all, field, ylab, out_dir, pID)
% One participant, lines per N-back level, x=periods (ordered)
    P = categorical(string(sub.period));
    present = categories(P);
    local_order = intersect(keep_order, present, 'stable');
    if isempty(local_order)
        local_order = present;
    end
    P = reordercats(P, local_order);
    sub.period = P;
    sessions = categories(P); x = 1:numel(sessions);

    f = figure('Color','k','Position',[220 80 1000 450]); hold on; box on; grid on;
    title(sprintf('Participant %s — %s', pID, ylab), 'Color', 'w');

    for L = levels_all'
        Lmask = (sub.level==L);
        Y = NaN(numel(sessions),1);
        if any(Lmask)
            tmp = sub(Lmask, {'period', field});
            [~,loc] = ismember(string(tmp.period), sessions);
            Y(loc) = tmp.(field);
        end
        plot(x, Y, '-o', 'LineWidth', 1.8, 'DisplayName', sprintf('Level %g', L));
    end

    set(gca,'XTick',x,'XTickLabel',sessions,'XTickLabelRotation',45,'Color','k','XColor','w','YColor','w');
    ylabel(ylab,'Color','w'); xlabel('Period','Color','w');
    legend('Location','best','TextColor','w','Color',[0.1 0.1 0.1]);

    if ~exist(out_dir,'dir'); mkdir(out_dir); end
    exportgraphics(f, fullfile(out_dir, sprintf('participant_%s_%s.png', pID, matlab.lang.makeValidName(ylab))), 'Resolution', 300);
    close(f);
end

function groupTrendByLevel(M, keep_order, levels_all, field, ylab, out_dir)
% Mean±SEM across participants; lines per level
    P = categorical(string(M.period));
    present = categories(P);
    local_order = intersect(keep_order, present, 'stable');
    if isempty(local_order)
        local_order = present;
    end
    P = reordercats(P, local_order);
    M.period = P;
    sessions = categories(P); x = 1:numel(sessions);

    f = figure('Color','k','Position',[220 100 1000 450]); hold on; box on; grid on;
    title(sprintf('Group trend — %s', ylab), 'Color', 'w');

    for L = levels_all'
        sub = M(M.level==L, {'participant','period',field});
        parts = unique(sub.participant,'stable');
        Y = NaN(numel(sessions), numel(parts));
        for ip = 1:numel(parts)
            sp = sub(sub.participant==parts(ip), :);
            [~,loc] = ismember(string(sp.period), sessions);
            Y(loc, ip) = sp.(field);
        end
        mu  = mean(Y, 2, 'omitnan');
        sem = std(Y, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(Y),2));

        if numel(x) >= 2
            px = [x, fliplr(x)];
            py = [mu-sem; flipud(mu+sem)]';
            patch(px, py, [0.7 0.7 0.7], 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility','off');
        end
        plot(x, mu, '-o', 'LineWidth', 2.2, 'DisplayName', sprintf('Level %g', L));
    end

    set(gca,'XTick',x,'XTickLabel',sessions,'XTickLabelRotation',45,'Color','k','XColor','w','YColor','w');
    ylabel(ylab,'Color','w'); xlabel('Period','Color','w');
    legend('Location','best','TextColor','w','Color',[0.1 0.1 0.1]);

    if ~exist(out_dir,'dir'); mkdir(out_dir); end
    exportgraphics(f, fullfile(out_dir, sprintf('group_%s.png', matlab.lang.makeValidName(ylab))), 'Resolution', 300);
    close(f);
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
