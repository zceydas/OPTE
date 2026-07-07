%% Path selection fix for build_lzc_v7_heatmap_outputs.m
%
% In build_lzc_v7_heatmap_outputs.m, replace the block that starts with:
%
%   input_root = uigetdir(...)
%   ...
%   v7_root = fullfile(input_root, 'LZc_Schartner_ProductionV7_Results');
%   csv_root = fullfile(v7_root, 'CSV');
%
% with this block. Then add the local function at the bottom of the file.

%% ---------------- REPLACEMENT BLOCK ----------------

selected_root = uigetdir(pwd, ['Select folder containing PostICA EEGLAB .set files, ' ...
    'or LZc_Schartner_ProductionV7_Results, or its CSV folder']);

if isequal(selected_root, 0)
    disp('Folder selection cancelled.');
    return;
end

[input_root, v7_root, csv_root] = resolve_v7_paths(selected_root);

if ~exist(csv_root, 'dir')
    error('V7 CSV folder not found: %s', csv_root);
end

%% ---------------- ADD THIS LOCAL FUNCTION AT THE BOTTOM ----------------

function [input_root, v7_root, csv_root] = resolve_v7_paths(selected_root)

selected_root = char(selected_root);
[parent_dir, folder_name] = fileparts(selected_root);

if strcmp(folder_name, 'CSV')
    csv_root = selected_root;
    v7_root = parent_dir;
    input_root = fileparts(v7_root);

elseif strcmp(folder_name, 'LZc_Schartner_ProductionV7_Results')
    v7_root = selected_root;
    csv_root = fullfile(v7_root, 'CSV');
    input_root = fileparts(v7_root);

else
    input_root = selected_root;
    v7_root = fullfile(input_root, 'LZc_Schartner_ProductionV7_Results');
    csv_root = fullfile(v7_root, 'CSV');
end

end
