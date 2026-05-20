function metrics_effort = effort_selection_from_structs
% effort_selection_from_structs.m
% Reads ALL Effort Selection .mat files under <root_dir>\Data (recursively).
% Expects filenames of the form:
%   Learning_subject_<ID>_Session<1-3>.mat
%   Decision_subject_<ID>_Session<1-3>.mat
%
% For each participant × period × EffortLevel (from 1–4), computes:
%
%   From Learning phase (Responses.Learning, 14 columns):
%       BlockNo, EffortLevel, Accuracy, RT, SwitchNumber, TrialType,
%       Response, Number, CorrectAnswer, Order, StartTrial, Trial,
%       TrialStartTime, TrialStartSince
%
%       → Using block-wise Order changes to classify trials as:
%             - Switch trial: Order(i) ~= Order(i-1) within same BlockNo
%             - Repeat trial: Order(i) == Order(i-1) within same BlockNo
%           (first trial in each block is treated as Repeat)
%
%       Metrics per EffortLevel:
%           n_switch_learning      (# valid switch trials, Accuracy < 9)
%           n_repeat_learning      (# valid repeat trials, Accuracy < 9)
%           switch_acc_learning    (prop correct on switch, Accuracy == 1)
%           repeat_acc_learning    (prop correct on repeat, Accuracy == 1)
%           switch_cost_acc_learning  = switch_acc - repeat_acc
%           switch_rt_learning     (mean RT on correct switch trials)
%           repeat_rt_learning     (mean RT on correct repeat trials)
%           switch_cost_rt_learning  = switch_rt - repeat_rt
%
%   From Decision phase (Responses.DecisionPhase.Selection, 12 columns):
%       BlockNo, LeftSide_Effort, RightSide_Effort, ChosenLevel,
%       ChosenSide, RT, MagParKind, OfferStartTime, OfferStartTimeSince,
%       ExecutionBeginTime, ExecutionEndTime, Run
%
%       Per EffortLevel:
%           n_offered              (# trials where that effort appeared on L or R)
%           n_chosen               (# trials where ChosenLevel == EffortLevel)
%           p_select               = n_chosen / n_offered (NaN if never offered)
%           decision_rt            mean RT on trials where ChosenLevel == EffortLevel
%
%   From Decision phase Execution (Responses.DecisionPhase.Execution, same 14
%   columns as Learning):
%       Same switch/repeat metrics as Learning, but labeled "..._execution".
%
% Outputs:
%   metrics_effort : table with one row per (participant × period × EffortLevel)
%   CSV            : <root_dir>\Results\effort_selection_metrics.csv
%   Figures        :
%       - Per-participant selection probabilities over periods for each EffortLevel
%       - Group-level selection probabilities (mean ± SEM) over periods
%       - Per-participant switch vs repeat accuracy/RT over periods (learning + execution)
%       - Group-level switch vs repeat accuracy/RT over periods (learning + execution)

%% PATHS
root_dir = 'C:\Users\nasak\Documents\Notes\JHU\CPCR\OPTE';
data_dir = fullfile(root_dir, 'Data');
out_dir  = fullfile(root_dir, 'Results');
if ~exist(data_dir,'dir'); error('Data folder not found: %s', data_dir); end
if ~exist(out_dir,'dir');  mkdir(out_dir); end

csv_out = fullfile(out_dir, 'effort_selection_metrics.csv');

%% Session → period mapping
ses_to_period = containers.Map( ...
    {1, 2, 3}, ...
    {'Baseline', '1-week', '1-month'} ...
);

PERIOD_ORDER = [ "Baseline", "1-week", "1-month" ];

%% FIND LEARNING & DECISION FILES
learn_files = dir(fullfile(data_dir, '**', 'Learning_*.mat'));
dec_files   = dir(fullfile(data_dir, '**', 'Decision_*.mat'));

if isempty(learn_files) && isempty(dec_files)
    error('No Learning_*.mat or Decision_*.mat files found under: %s', data_dir);
end

learn_map = containers.Map('KeyType','char','ValueType','char');
dec_map   = containers.Map('KeyType','char','ValueType','char');

for k = 1:numel(learn_files)
    fname = learn_files(k).name;
    fpath = fullfile(learn_files(k).folder, fname);
    [pid, sesnum] = parseEffortFilename(fname);
    if isnan(sesnum) || isempty(pid)
        fprintf('Skipping Learning (cannot parse subject/session): %s\n', fpath);
        continue;
    end
    key = makeKey(pid, sesnum);
    learn_map(key) = fpath;
end

for k = 1:numel(dec_files)
    fname = dec_files(k).name;
    fpath = fullfile(dec_files(k).folder, fname);
    [pid, sesnum] = parseEffortFilename(fname);
    if isnan(sesnum) || isempty(pid)
        fprintf('Skipping Decision (cannot parse subject/session): %s\n', fpath);
        continue;
    end
    key = makeKey(pid, sesnum);
    dec_map(key) = fpath;
end

all_keys = unique([learn_map.keys, dec_map.keys]);

if isempty(all_keys)
    error('No participant-session mapping could be built from filenames.');
end

%% PER (participant × session) PROCESSING
rows = [];
for ik = 1:numel(all_keys)
    key = all_keys{ik};
    [pid, sesnum] = parseKey(key);
    if isnan(sesnum) || ~isKey(ses_to_period, sesnum)
        fprintf('Skipping key %s (no valid session mapping)\n', key);
        continue;
    end
    period = string(ses_to_period(sesnum));

    learning_metrics = [];
    if isKey(learn_map, key)
        fpathL = learn_map(key);
        S = load(fpathL);
        if ~isfield(S, 'Results')
            warning('File %s has no Results struct. Skipping Learning.', fpathL);
        else
            try
                M_learn = extractLearningMatrix(S.Results);
                learning_metrics = computeSwitchRepeatMetrics(M_learn);
            catch ME
                warning('Error processing Learning in %s: %s', fpathL, ME.message);
            end
        end
    end

    selection_metrics = [];
    exec_metrics      = [];
    if isKey(dec_map, key)
        fpathD = dec_map(key);
        S = load(fpathD);
        if ~isfield(S, 'Results')
            warning('File %s has no Results struct. Skipping Decision.', fpathD);
        else
            try
                [M_sel, M_exec] = extractDecisionMatrices(S.Results);
                selection_metrics = computeSelectionMetrics(M_sel);
                exec_metrics      = computeSwitchRepeatMetrics(M_exec);
            catch ME
                warning('Error processing Decision in %s: %s', fpathD, ME.message);
            end
        end
    end

    if isempty(learning_metrics) && isempty(selection_metrics) && isempty(exec_metrics)
        fprintf('No usable data for %s (participant %s, session %d).\n', key, pid, sesnum);
        continue;
    end

    eff_all = [];
    if ~isempty(learning_metrics)
        eff_all = [eff_all; [learning_metrics.effort_level]'];
    end
    if ~isempty(selection_metrics)
        eff_all = [eff_all; [selection_metrics.effort_level]'];
    end
    if ~isempty(exec_metrics)
        eff_all = [eff_all; [exec_metrics.effort_level]'];
    end
    eff_levels = unique(eff_all);

    for ie = 1:numel(eff_levels)
        E = eff_levels(ie);

        L    = findRowByEffort(learning_metrics, E);
        Ssel = findRowByEffort(selection_metrics, E);
        Lexe = findRowByEffort(exec_metrics, E);

        row.participant   = string(pid);
        row.sesnum        = double(sesnum);
        row.period        = string(period);
        row.effort_level  = double(E);

        row.n_switch_learning           = getFieldOrNaN(L, 'n_switch');
        row.n_repeat_learning           = getFieldOrNaN(L, 'n_repeat');
        row.switch_acc_learning         = getFieldOrNaN(L, 'switch_acc');
        row.repeat_acc_learning         = getFieldOrNaN(L, 'repeat_acc');
        row.switch_cost_acc_learning    = getFieldOrNaN(L, 'switch_cost_acc');
        row.switch_rt_learning          = getFieldOrNaN(L, 'switch_rt');
        row.repeat_rt_learning          = getFieldOrNaN(L, 'repeat_rt');
        row.switch_cost_rt_learning     = getFieldOrNaN(L, 'switch_cost_rt');

        row.n_switch_execution          = getFieldOrNaN(Lexe, 'n_switch');
        row.n_repeat_execution          = getFieldOrNaN(Lexe, 'n_repeat');
        row.switch_acc_execution        = getFieldOrNaN(Lexe, 'switch_acc');
        row.repeat_acc_execution        = getFieldOrNaN(Lexe, 'repeat_acc');
        row.switch_cost_acc_execution   = getFieldOrNaN(Lexe, 'switch_cost_acc');
        row.switch_rt_execution         = getFieldOrNaN(Lexe, 'switch_rt');
        row.repeat_rt_execution         = getFieldOrNaN(Lexe, 'repeat_rt');
        row.switch_cost_rt_execution    = getFieldOrNaN(Lexe, 'switch_cost_rt');

        row.n_offered   = getFieldOrNaN(Ssel, 'n_offered');
        row.n_chosen    = getFieldOrNaN(Ssel, 'n_chosen');
        row.p_select    = getFieldOrNaN(Ssel, 'p_select');
        row.decision_rt = getFieldOrNaN(Ssel, 'decision_rt');

        rows = [rows; row]; %#ok<AGROW>
    end
end

if isempty(rows)
    error('No valid Effort Selection metrics could be computed from available files.');
end

metrics_effort = struct2table(rows);

%% ORDERING & CSV
metrics_effort.period = categorical(string(metrics_effort.period));
present = categories(metrics_effort.period);
desired = cellstr(PERIOD_ORDER(:)');
keep_order = intersect(desired, present, 'stable');
metrics_effort.period = reordercats(metrics_effort.period, keep_order);

metrics_effort = sortrows(metrics_effort, {'participant','period','effort_level'});

try
    writetable(metrics_effort, csv_out);
catch
    warning('Could not write CSV to %s (is the file open?)', csv_out);
end

%% PLOTTING: SELECTION PROBABILITY
participants = unique(metrics_effort.participant, 'stable');
eff_levels   = unique(metrics_effort.effort_level);
eff_levels(eff_levels == 999) = [];

for ip = 1:numel(participants)
    pid = participants(ip);
    sub = metrics_effort(metrics_effort.participant == pid, :);
    sub = sortrows(sub, {'period','effort_level'});

    periods = categories(sub.period);
    x = 1:numel(periods);

    f = figure('Color','k','Position',[200 200 900 600]); hold on; grid on;
    title(sprintf('Participant %s — P(select EffortLevel) over periods', pid), 'Interpreter','none');
    xlabel('Period');
    ylabel('P(select)');

    cmap = lines(numel(eff_levels));
    for ie = 1:numel(eff_levels)
        E = eff_levels(ie);
        mask = sub.effort_level == E;

        y = NaN(numel(periods),1);
        for iper = 1:numel(periods)
            maskP = mask & (sub.period == periods{iper});
            y(iper) = mean(sub.p_select(maskP), 'omitnan');
        end

        plot(x, y, '-o', 'LineWidth', 1.5, 'Color', cmap(ie,:), ...
            'DisplayName', sprintf('Effort %d', E));
    end
    set(gca,'XTick',x,'XTickLabel',periods,'XTickLabelRotation',30);
    legend('Location','best');
    saveas(f, fullfile(out_dir, sprintf('participant_%s_selection_probabilities.png', pid)));
    close(f);
end

plotEffortGroupTrend(metrics_effort, eff_levels, PERIOD_ORDER, out_dir);
plotSwitchRepeatAcrossPeriods(metrics_effort, eff_levels, PERIOD_ORDER, out_dir);

fprintf('Done.\nEffort Selection metrics (CSV & figures) -> %s\n', out_dir);

end

%% HELPERS
function [pid, sesnum] = parseEffortFilename(name)
    pid = "";
    sesnum = NaN;
    try
        tok = regexp(name, 'subject[_-]?(\d+).*Session[_-]?(\d+)', 'tokens', 'once');
        if ~isempty(tok)
            pid = string(tok{1});
            sesnum = str2double(tok{2});
        end
    catch
    end
end

function key = makeKey(pid, sesnum)
    key = sprintf('%s_ses%d', pid, sesnum);
end

function [pid, sesnum] = parseKey(key)
    tok = regexp(key, '^(.+)_ses(\d+)$', 'tokens', 'once');
    if isempty(tok)
        pid = "";
        sesnum = NaN;
    else
        pid = string(tok{1});
        sesnum = str2double(tok{2});
    end
end

function M = extractLearningMatrix(Results)
    M = [];
    if ~isfield(Results, 'Subject')
        error('Results struct has no field "Subject".');
    end
    subs = Results.Subject;
    for i = 1:numel(subs)
        sub = subs(i);
        if isfield(sub, 'Session') && ~isempty(sub.Session)
            sessArr = sub.Session;
            for j = 1:numel(sessArr)
                sess = sessArr(j);
                if isfield(sess, 'Responses') && isfield(sess.Responses, 'Learning') ...
                        && ~isempty(sess.Responses.Learning)
                    M = sess.Responses.Learning;
                    return;
                end
            end
        end
    end
    if isempty(M)
        error('Could not find Responses.Learning in Results.Subject.Session.');
    end
end

function [Sel, Exec] = extractDecisionMatrices(Results)
    Sel  = [];
    Exec = [];
    if ~isfield(Results, 'Subject')
        error('Results struct has no field "Subject".');
    end
    subs = Results.Subject;
    for i = 1:numel(subs)
        sub = subs(i);
        if isfield(sub, 'Session') && ~isempty(sub.Session)
            sessArr = sub.Session;
            for j = 1:numel(sessArr)
                sess = sessArr(j);
                if isfield(sess, 'Responses') && isfield(sess.Responses, 'DecisionPhase')
                    DP = sess.Responses.DecisionPhase;
                    if isfield(DP, 'Selection') && ~isempty(DP.Selection)
                        Sel = DP.Selection;
                    end
                    if isfield(DP, 'Execution') && ~isempty(DP.Execution)
                        Exec = DP.Execution;
                    end
                    if ~isempty(Sel) || ~isempty(Exec)
                        return;
                    end
                end
            end
        end
    end
    if isempty(Sel) && isempty(Exec)
        error('Could not find DecisionPhase.Selection/Execution in Results.Subject.Session.');
    end
end

function metrics = computeSwitchRepeatMetrics(M)
    if isempty(M)
        metrics = [];
        return;
    end

    varNames = { ...
        'BlockNo','EffortLevel','Accuracy','RT','SwitchNumber','TrialType', ...
        'Response','Number','CorrectAnswer','Order','StartTrial','Trial', ...
        'TrialStartTime','TrialStartSince'};

    if size(M,2) ~= numel(varNames)
        error('Expected 14 columns for Learning/Execution matrix, got %d.', size(M,2));
    end

    T = array2table(M, 'VariableNames', varNames);

    n = height(T);
    if n == 0
        metrics = [];
        return;
    end

    valid   = T.Accuracy < 9;
    correct = T.Accuracy == 1;

    isSwitch = false(n,1);
    blocks = unique(T.BlockNo(~isnan(T.BlockNo)));
    for ib = 1:numel(blocks)
        b = blocks(ib);
        idx = find(T.BlockNo == b);
        if numel(idx) <= 1
            continue;
        end
        ord = T.Order(idx);
        swLocal = [false; diff(ord) ~= 0];
        isSwitch(idx) = swLocal;
    end
    isRepeat = ~isSwitch;

    eff_levels = unique(T.EffortLevel(~isnan(T.EffortLevel)));

    metrics = struct( ...
        'effort_level', {}, ...
        'n_switch', {}, ...
        'n_repeat', {}, ...
        'switch_acc', {}, ...
        'repeat_acc', {}, ...
        'switch_cost_acc', {}, ...
        'switch_rt', {}, ...
        'repeat_rt', {}, ...
        'switch_cost_rt', {} );

    for ie = 1:numel(eff_levels)
        E = eff_levels(ie);
        idxEff = (T.EffortLevel == E) & valid;

        swIdx  = idxEff & isSwitch;
        repIdx = idxEff & isRepeat;

        n_sw  = sum(swIdx);
        n_rep = sum(repIdx);

        if n_sw > 0
            sw_acc = sum(swIdx & correct) / n_sw;
            sw_rt  = mean(T.RT(swIdx & correct), 'omitnan');
        else
            sw_acc = NaN;
            sw_rt  = NaN;
        end

        if n_rep > 0
            rep_acc = sum(repIdx & correct) / n_rep;
            rep_rt  = mean(T.RT(repIdx & correct), 'omitnan');
        else
            rep_acc = NaN;
            rep_rt  = NaN;
        end

        row.effort_level    = double(E);
        row.n_switch        = double(n_sw);
        row.n_repeat        = double(n_rep);
        row.switch_acc      = sw_acc;
        row.repeat_acc      = rep_acc;
        row.switch_cost_acc = sw_acc - rep_acc;
        row.switch_rt       = sw_rt;
        row.repeat_rt       = rep_rt;
        row.switch_cost_rt  = sw_rt - rep_rt;

        metrics(end+1,1) = row; %#ok<AGROW>
    end
end

function metrics = computeSelectionMetrics(Msel)
    if isempty(Msel)
        metrics = [];
        return;
    end

    varNames = { ...
        'BlockNo','LeftSide_Effort','RightSide_Effort','ChosenLevel', ...
        'ChosenSide','RT','MagParKind','OfferStartTime','OfferStartTimeSince', ...
        'ExecutionBeginTime','ExecutionEndTime','Run'};

    if size(Msel,2) ~= numel(varNames)
        error('Expected 12 columns for Selection matrix, got %d.', size(Msel,2));
    end

    T = array2table(Msel, 'VariableNames', varNames);

    eff_all = [T.LeftSide_Effort; T.RightSide_Effort; T.ChosenLevel];
    eff_levels = unique(eff_all(~isnan(eff_all)));

    metrics = struct( ...
        'effort_level', {}, ...
        'n_offered', {}, ...
        'n_chosen', {}, ...
        'p_select', {}, ...
        'decision_rt', {} );

    for ie = 1:numel(eff_levels)
        E = eff_levels(ie);

        offered = (T.LeftSide_Effort == E) | (T.RightSide_Effort == E);
        chosen  = (T.ChosenLevel == E);

        n_offered = sum(offered);
        n_chosen  = sum(chosen);

        if n_offered > 0
            p_select = n_chosen / n_offered;
        else
            p_select = NaN;
        end

        if n_chosen > 0
            decision_rt = mean(T.RT(chosen), 'omitnan');
        else
            decision_rt = NaN;
        end

        row.effort_level = double(E);
        row.n_offered    = double(n_offered);
        row.n_chosen     = double(n_chosen);
        row.p_select     = p_select;
        row.decision_rt  = decision_rt;

        metrics(end+1,1) = row; %#ok<AGROW>
    end
end

function row = findRowByEffort(metrics, E)
    row = [];
    if isempty(metrics)
        return;
    end
    for i = 1:numel(metrics)
        if isfield(metrics(i), 'effort_level') && metrics(i).effort_level == E
            row = metrics(i);
            return;
        end
    end
end

function v = getFieldOrNaN(s, fld)
    if isempty(s) || ~isfield(s, fld) || isempty(s.(fld))
        v = NaN;
    else
        v = s.(fld);
    end
end

function plotEffortGroupTrend(M, eff_levels, PERIOD_ORDER, out_dir)
    if isempty(M)
        return;
    end

    eff_levels(eff_levels == 999) = [];

    P = categorical(string(M.period));
    present = categories(P);
    desired = cellstr(PERIOD_ORDER(:)');
    keep_order = intersect(desired, present, 'stable');
    P = reordercats(P, keep_order);
    M.period = P;

    periods = categories(M.period);
    x = 1:numel(periods);

    f = figure('Color','k','Position',[200 200 900 600]); hold on; grid on;
    title('Group selection probability — Effort levels over periods');
    xlabel('Period');
    ylabel('P(select)');

    cmap = lines(numel(eff_levels));

    for ie = 1:numel(eff_levels)
        E = eff_levels(ie);
        y   = NaN(numel(periods),1);
        sem = NaN(numel(periods),1);

        for ip = 1:numel(periods)
            mask = (M.effort_level == E) & (M.period == periods{ip});
            vals = M.p_select(mask);
            y(ip)   = mean(vals, 'omitnan');
            sem(ip) = std(vals, 0, 'omitnan') ./ sqrt(sum(~isnan(vals)));
        end

        errorbar(x, y, sem, '-o', 'LineWidth', 1.5, 'Color', cmap(ie,:), ...
                 'DisplayName', sprintf('Effort %d', E));
    end

    set(gca, 'XTick', x, 'XTickLabel', periods, 'XTickLabelRotation', 30);
    legend('Location','best');

    saveas(f, fullfile(out_dir, 'group_selection_probabilities.png'));
    close(f);
end

function plotSwitchRepeatAcrossPeriods(M, eff_levels, PERIOD_ORDER, out_dir)

    if isempty(M)
        return;
    end

    eff_levels(eff_levels == 999) = [];

    P = categorical(string(M.period));
    present = categories(P);
    desired = cellstr(PERIOD_ORDER(:)');
    keep_order = intersect(desired, present, 'stable');
    P = reordercats(P, keep_order);
    M.period = P;

    periods = categories(M.period);
    x = 1:numel(periods);

    participants = unique(M.participant, 'stable');

    phases = { ...
        'learning', ...
        'execution' ...
    };

    for ph = 1:numel(phases)
        phase = phases{ph};

        sw_acc_col = sprintf('switch_acc_%s', phase);
        rp_acc_col = sprintf('repeat_acc_%s', phase);
        sw_rt_col  = sprintf('switch_rt_%s',  phase);
        rp_rt_col  = sprintf('repeat_rt_%s',  phase);

        if ~ismember(sw_acc_col, M.Properties.VariableNames)
            continue;
        end

        for ip = 1:numel(participants)
            pid = participants(ip);
            subP = M(M.participant == pid, :);

            for ie = 1:numel(eff_levels)
                E = eff_levels(ie);
                sub = subP(subP.effort_level == E, :);
                if isempty(sub); continue; end

                y_sw_acc = NaN(numel(periods),1);
                y_rp_acc = NaN(numel(periods),1);
                y_sw_rt  = NaN(numel(periods),1);
                y_rp_rt  = NaN(numel(periods),1);

                for k = 1:numel(periods)
                    mask = sub.period == periods{k};
                    y_sw_acc(k) = mean(sub.(sw_acc_col)(mask), 'omitnan');
                    y_rp_acc(k) = mean(sub.(rp_acc_col)(mask), 'omitnan');
                    y_sw_rt(k)  = mean(sub.(sw_rt_col)(mask),  'omitnan');
                    y_rp_rt(k)  = mean(sub.(rp_rt_col)(mask),  'omitnan');
                end

                f1 = figure('Color','k','Position',[200 200 950 450]); hold on; grid on;
                title(sprintf('Participant %s — %s accuracy (Effort %d)', pid, phase, E), 'Interpreter','none');
                plot(x, y_sw_acc, '-o', 'LineWidth', 1.8, 'DisplayName','Switch');
                plot(x, y_rp_acc, '-o', 'LineWidth', 1.8, 'DisplayName','Repeat');
                set(gca,'XTick',x,'XTickLabel',periods,'XTickLabelRotation',30);
                xlabel('Period'); ylabel('Accuracy'); ylim([0 1]);
                legend('Location','best');
                saveas(f1, fullfile(out_dir, sprintf('participant_%s_%s_accuracy_effort_%d.png', pid, phase, E)));
                close(f1);

                f2 = figure('Color','k','Position',[200 200 950 450]); hold on; grid on;
                title(sprintf('Participant %s — %s RT (Effort %d)', pid, phase, E), 'Interpreter','none');
                plot(x, y_sw_rt, '-o', 'LineWidth', 1.8, 'DisplayName','Switch');
                plot(x, y_rp_rt, '-o', 'LineWidth', 1.8, 'DisplayName','Repeat');
                set(gca,'XTick',x,'XTickLabel',periods,'XTickLabelRotation',30);
                xlabel('Period'); ylabel('RT');
                legend('Location','best');
                saveas(f2, fullfile(out_dir, sprintf('participant_%s_%s_rt_effort_%d.png', pid, phase, E)));
                close(f2);
            end
        end

        for ie = 1:numel(eff_levels)
            E = eff_levels(ie);

            y_sw_acc = NaN(numel(periods),1);
            y_rp_acc = NaN(numel(periods),1);
            se_sw_acc = NaN(numel(periods),1);
            se_rp_acc = NaN(numel(periods),1);

            y_sw_rt = NaN(numel(periods),1);
            y_rp_rt = NaN(numel(periods),1);
            se_sw_rt = NaN(numel(periods),1);
            se_rp_rt = NaN(numel(periods),1);

            for k = 1:numel(periods)
                mask = (M.effort_level == E) & (M.period == periods{k});

                v_sw_acc = M.(sw_acc_col)(mask);
                v_rp_acc = M.(rp_acc_col)(mask);
                v_sw_rt  = M.(sw_rt_col)(mask);
                v_rp_rt  = M.(rp_rt_col)(mask);

                y_sw_acc(k) = mean(v_sw_acc, 'omitnan');
                y_rp_acc(k) = mean(v_rp_acc, 'omitnan');
                se_sw_acc(k) = std(v_sw_acc, 0, 'omitnan') ./ sqrt(sum(~isnan(v_sw_acc)));
                se_rp_acc(k) = std(v_rp_acc, 0, 'omitnan') ./ sqrt(sum(~isnan(v_rp_acc)));

                y_sw_rt(k) = mean(v_sw_rt, 'omitnan');
                y_rp_rt(k) = mean(v_rp_rt, 'omitnan');
                se_sw_rt(k) = std(v_sw_rt, 0, 'omitnan') ./ sqrt(sum(~isnan(v_sw_rt)));
                se_rp_rt(k) = std(v_rp_rt, 0, 'omitnan') ./ sqrt(sum(~isnan(v_rp_rt)));
            end

            f3 = figure('Color','k','Position',[200 200 950 450]); hold on; grid on;
            title(sprintf('Group — %s accuracy (Effort %d)', phase, E));
            errorbar(x, y_sw_acc, se_sw_acc, '-o', 'LineWidth', 1.8, 'DisplayName','Switch');
            errorbar(x, y_rp_acc, se_rp_acc, '-o', 'LineWidth', 1.8, 'DisplayName','Repeat');
            set(gca,'XTick',x,'XTickLabel',periods,'XTickLabelRotation',30);
            xlabel('Period'); ylabel('Accuracy'); ylim([0 1]);
            legend('Location','best');
            saveas(f3, fullfile(out_dir, sprintf('group_%s_accuracy_effort_%d.png', phase, E)));
            close(f3);

            f4 = figure('Color','k','Position',[200 200 950 450]); hold on; grid on;
            title(sprintf('Group — %s RT (Effort %d)', phase, E));
            errorbar(x, y_sw_rt, se_sw_rt, '-o', 'LineWidth', 1.8, 'DisplayName','Switch');
            errorbar(x, y_rp_rt, se_rp_rt, '-o', 'LineWidth', 1.8, 'DisplayName','Repeat');
            set(gca,'XTick',x,'XTickLabel',periods,'XTickLabelRotation',30);
            xlabel('Period'); ylabel('RT');
            legend('Location','best');
            saveas(f4, fullfile(out_dir, sprintf('group_%s_rt_effort_%d.png', phase, E)));
            close(f4);
        end
    end
end
