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

%% PATHS
root_dir = 'C:\Users\nasak\Documents\Notes\JHU\CPCR\OPTE';


data_dir = fullfile(root_dir, 'Data');
out_dir  = fullfile(root_dir, 'Results');
if ~exist(data_dir,'dir'); error('Data folder not found: %s', data_dir); end
if ~exist(out_dir,'dir');  mkdir(out_dir); end

csv_out = fullfile(out_dir, 'gonogo_metrics_by_period.csv');

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
files = dir(fullfile(data_dir, '**', '*.mat'));
if isempty(files)
    error('No .mat files found under: %s', data_dir);
end

%% PER-FILE COUNTS
rows = [];
skipped_map = 0;
skipped_nores = 0;

for k = 1:numel(files)
    fpath = fullfile(files(k).folder, files(k).name);

    % Parse subject ID and ses-# from filename (e.g., Subject_0005_ses-5_CCPT_*.mat)
    [pid, sesnum] = parseFilename(files(k).name);

    % Map ses-# to period label
    if ~isnan(sesnum) && isKey(ses_to_period, sesnum)
        period = string(ses_to_period(sesnum));
    else
        skipped_map = skipped_map + 1;
        fprintf('Skipping (cannot map ses-#): %s\n', fpath);
        continue;
    end

    % Load .mat; find results struct
    S = load(fpath);
    if isfield(S, fld.results)
        R = S.(fld.results);
    else
        R = pickFirstStructArray(S);
    end
    if ~isstruct(R) || isempty(R)
        skipped_nores = skipped_nores + 1;
        fprintf('Skipping (no results struct): %s\n', fpath);
        continue;
    end

    % Subject inside file overrides filename
    if isfield(S, fld.subject) && ~isempty(S.(fld.subject))
        pid = string(S.(fld.subject));
    end

    % Convert to table
    Rt = struct2table(R);
    mustHave(Rt, fld.isGo,     fpath);
    mustHave(Rt, fld.response, fpath);
    mustHave(Rt, fld.correct,  fpath);
    mustHave(Rt, fld.RT,       fpath);

    isGo      = logical(Rt.(fld.isGo));
    responded = logical(Rt.(fld.response));
    correct   = logical(Rt.(fld.correct));
    RT        = Rt.(fld.RT);

    % Counts per file
    nGo    = sum(isGo == 1);
    nNoGo  = sum(isGo == 0);
    hits   = sum(isGo & correct);                % correct Go
    fas    = sum(~isGo & responded);             % responses on No-Go
    rt_sum = sum(RT(isGo & correct), 'omitnan'); % for weighted mean

    rows = [rows; makeFileRow(pid, period, nGo, nNoGo, hits, fas, rt_sum, fpath, sesnum)];
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
    writetable(metrics_by_period, csv_out);
catch
    warning('Could not write CSV to %s (is the file open?)', csv_out);
end

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

    saveas(f, fullfile(out_dir, sprintf('participant_%s_timeseries.png', participants(ip))));
    close(f);
end

% Group trend plots (use only periods actually present)
plotGroupTrend(metrics_by_period, 'H',         'Hit Rate',      keep_order, out_dir);
plotGroupTrend(metrics_by_period, 'FA',        'False Alarm',   keep_order, out_dir);
plotGroupTrend(metrics_by_period, 'dprime',    "d'",            keep_order, out_dir);
plotGroupTrend(metrics_by_period, 'criterion', 'Criterion (c)', keep_order, out_dir);
plotGroupTrend(metrics_by_period, 'mean_rt',   'Mean RT',       keep_order, out_dir);

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

    saveas(f, fullfile(out_dir, sprintf('group_trend_%s.png', varname)));
    close(f);
end
