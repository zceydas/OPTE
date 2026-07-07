function multitaper_group_spectrogram_figures(spectrogram_index, save_dir)

if isempty(spectrogram_index)
    warning('No spectrogram index available for group spectrogram figures.');
    return;
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

T = normalize_psd_table_ids(spectrogram_index);
conditions = psd_condition_order();
eyes_list = unique(T.Eyes, 'stable');

for e = 1:numel(eyes_list)
    eyes_name = string(eyes_list(e));
    fig_path = fullfile(save_dir, ['GROUP_multitaper_spectrogram_' char(eyes_name) '.png']);
    mat_path = fullfile(save_dir, ['GROUP_multitaper_spectrogram_' char(eyes_name) '.mat']);

    if exist(fig_path, 'file') && exist(mat_path, 'file')
        fprintf('Group spectrogram figure already exists. Skipping:\n%s\n', fig_path);
        continue;
    end

    group_specs = cell(height(conditions), 1);
    group_times = cell(height(conditions), 1);
    group_freqs = cell(height(conditions), 1);
    n_files = nan(height(conditions), 1);

    for c = 1:height(conditions)
        subset = T(T.Eyes == eyes_name & ...
                   T.Session == conditions.Session(c) & ...
                   T.Epoch == conditions.Epoch(c), :);

        if isempty(subset)
            continue;
        end

        loaded = cell(height(subset), 1);
        minF = inf;
        minT = inf;

        for i = 1:height(subset)
            temp = load(subset.SpectrogramMAT(i), 'spec_power', 'spec_times', 'spec_freqs');
            loaded{i} = temp;
            minF = min(minF, size(temp.spec_power, 1));
            minT = min(minT, size(temp.spec_power, 2));
        end

        all_log_specs = nan(minF, minT, height(subset));

        for i = 1:height(subset)
            all_log_specs(:,:,i) = log10(loaded{i}.spec_power(1:minF, 1:minT) + eps);
        end

        group_specs{c} = mean(all_log_specs, 3, 'omitnan');
        group_times{c} = loaded{1}.spec_times(1:minT);
        group_freqs{c} = loaded{1}.spec_freqs(1:minF);
        n_files(c) = height(subset);
    end

    if all(cellfun(@isempty, group_specs))
        continue;
    end

    save(mat_path, 'group_specs', 'group_times', 'group_freqs', 'conditions', 'eyes_name', 'n_files', '-v7.3');

    fig = figure('Color','w','Position',[50 50 1800 600]);
    set(fig, 'InvertHardcopy','off');
    tiledlayout(1, height(conditions), 'TileSpacing','compact', 'Padding','compact');

    for c = 1:height(conditions)
        nexttile;

        if isempty(group_specs{c})
            axis off;
            title([char(conditions.Condition(c)) newline 'No data'], 'Interpreter','none', 'Color','k');
            continue;
        end

        imagesc(group_times{c}, group_freqs{c}, group_specs{c});
        axis xy;
        xlabel('Time (s)', 'Color','k');
        ylabel('Frequency (Hz)', 'Color','k');
        title(sprintf('%s\nN=%d', char(conditions.Condition(c)), n_files(c)), 'Interpreter','none', 'Color','k');
        set(gca, 'Color','w', 'XColor','k', 'YColor','k');
    end

    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Color = 'k';
    ylabel(cb, 'Mean log10 power');

    sgtitle(['Group multitaper spectrogram - ' char(eyes_name)], 'Interpreter','none', 'Color','k');
    saveas(fig, fig_path);
    close(fig);
end

end
