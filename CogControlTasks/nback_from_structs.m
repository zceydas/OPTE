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

%% ---------- PATHS ----------
root_dir = 'C:\Users\nasak\Documents\Notes\JHU\CPCR\OPTE';

data_dir = fullfile(root_dir, 'Data');
out_dir  = fullfile(root_dir, 'Results');
if ~exist(data_dir,'dir'); error('Data folder not found: %s', data_dir); end
if ~exist(out_dir,'dir');  mkdir(out_dir); end

csv_out_metrics = fullfile(out_dir, 'nback_metrics_by_period_level.csv');
csv_out_filemap = fullfile(out_dir, 'nback_filemap.csv');

% Subject ID formatting
SUBJECT_PAD_WIDTH = 4;  % e.g., 5 -> "0005"

% Session → period mapping (edit if needed)
ses_to_period = containers.Map( ...
    {1,                   2,               3,                  4,                   5}, ...
    {'Baseline TMSEEG',   'Dosing Session','1-week follow-up', '2-week follow-up',  '1-month follow-up'} ...
);
PERIOD_ORDER = unique(string(values(ses_to_period)),'stable');

%% ---------- FIND FILES ----------
files = dir(fullfile(data_dir, '**', '*.mat'));
if isempty(files)
    error('No .mat files found under: %s', data_dir);
end

%% ---------- PREFLIGHT: file → subject/session/period map ----------
frows = [];
for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);
    [pid_raw, sesnum] = parseFilename(files(k).name);
    pid_fmt = normalizeSubject(pid_raw, SUBJECT_PAD_WIDTH);

    mapped = true; reason = "";
    period = "";
    if ~isnan(sesnum) && isKey(ses_to_period, sesnum)
        period = string(ses_to_period(sesnum));
    else
        mapped = false;
        if isnan(sesnum)
            reason = "No session number found in filename";
        else
            reason = "Session number not in ses_to_period map";
        end
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

% Save metrics CSV
writetable(metrics, csv_out_metrics);
fprintf('Metrics saved: %s\n', csv_out_metrics);

%% ---------- PLOTTING ----------
parts  = unique(metrics.participant,'stable');
levels_all = unique(metrics.level,'stable');

for ip = 1:numel(parts)
    sub = metrics(metrics.participant==parts(ip), :);
    sub = sortrows(sub, {'period','level'});

    participantMultiLine(sub, keep_order, levels_all, 'dprime',       "d'",            out_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'hit_rate_adj', 'Hit rate (adj)',out_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'fa_rate_adj',  'FA rate (adj)', out_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'mean_rt_hit',  'Mean RT (Hits) (ms)',out_dir, parts(ip));
    participantMultiLine(sub, keep_order, levels_all, 'mean_rt_fa',   'Mean RT (FAs) (ms)', out_dir, parts(ip));
end

groupTrendByLevel(metrics, keep_order, levels_all, 'dprime',       "d'",            out_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'hit_rate_adj', 'Hit rate (adj)',out_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'fa_rate_adj',  'FA rate (adj)', out_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'mean_rt_hit',  'Mean RT (Hits)',out_dir);
groupTrendByLevel(metrics, keep_order, levels_all, 'mean_rt_fa',   'Mean RT (FAs)', out_dir);

fprintf('Done.\nCSV & figures -> %s\n', out_dir);
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
    P = reordercats(P, keep_order);
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
    saveas(f, fullfile(out_dir, sprintf('participant_%s_%s.png', pID, matlab.lang.makeValidName(ylab))));
    close(f);
end

function groupTrendByLevel(M, keep_order, levels_all, field, ylab, out_dir)
% Mean±SEM across participants; lines per level
    P = categorical(string(M.period));
    P = reordercats(P, keep_order);
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
    saveas(f, fullfile(out_dir, sprintf('group_%s.png', matlab.lang.makeValidName(ylab))));
    close(f);
end