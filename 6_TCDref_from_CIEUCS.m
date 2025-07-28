% Script: plot_DeltaE_from_MultiFiles.m
% Description: Process updated L*a*b* % files and plot ΔE*ab with spacing fix

clear; clc;

%% File Setup
base_path = uigetdir(pwd, 'Select the folder with CAM02UCS_sRGB_D65_XXX.txt files');
if base_path == 0
    error('No folder selected.');
end

file_suffixes = [300, 600, 900, 1200, 1500];
all_thickness = [];
all_deltaE = [];

ref_file = fullfile(base_path, 'CAM02UCS_sRGB_D65_300.txt');
if ~isfile(ref_file)
    error('Reference file Lab_D65_300.txt not found.');
end

%% Extract reference from 0% of Lab_D65_300.txt
[L0, a0, b0] = get_reference_values(ref_file);

%% Process each file
for i = 1:length(file_suffixes)
    suffix = file_suffixes(i);
    filename = sprintf('CAM02UCS_sRGB_D65_%d.txt', suffix);
    filepath = fullfile(base_path, filename);
    
    if ~isfile(filepath)
        warning('File not found: %s', filename);
        continue;
    end

    % Parse percentage-based file
    [percentages, L_vals, a_vals, b_vals] = parse_nm_file(filepath);

    % Remove 0% (reference) for all files except the first
    if i > 1
        mask = percentages > 0;
        percentages = percentages(mask);
        L_vals = L_vals(mask);
        a_vals = a_vals(mask);
        b_vals = b_vals(mask);
    end

    % Convert % to nm: (percent/100)*300 + offset
    offset = suffix - 300;
    thickness_nm = (percentages / 100) * 300 + offset;

    % Compute ΔE*ab
    delta_E = sqrt((L_vals - L0).^2 + (a_vals - a0).^2 + (b_vals - b0).^2);

    % Accumulate
    all_thickness = [all_thickness; thickness_nm];
    all_deltaE = [all_deltaE; delta_E];
end

%% Sort by thickness
[all_thickness, sort_idx] = sort(all_thickness);
all_deltaE = all_deltaE(sort_idx);

%% Ask user if normalization is required
normalize_choice = input('Normalize colorbar to [0, 1]? (y/n): ', 's');
normalize_flag = strcmpi(normalize_choice, 'y');

if normalize_flag
    delta_E_plot = (all_deltaE - min(all_deltaE)) / (max(all_deltaE) - min(all_deltaE));
else
    delta_E_plot = all_deltaE;
end


%% Plot heatmap
figure;
imagesc(all_thickness', [min(delta_E_plot) max(delta_E_plot)], delta_E_plot');  % 1D heatbar
colormap("turbo");
cb = colorbar;
xlabel('Effective Packing', 'Interpreter', 'latex');
yticks([]);
title('$\Delta E^*_{ab}$ vs. Thickness ($PS$)', 'Interpreter','latex','FontWeight','bold');
set(gca, 'FontSize', 12);

% Set X-ticks dynamically at 300 nm intervals
xticks(0:300:max(all_thickness));

cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';

if normalize_flag
    cb.Ticks = 0:0.05:1;
    cb.Label.String = '$\Delta E^*_{ab}$ (Normalized)';
else
    cb.Label.String = '$\Delta E^*_{ab}$';
end

%% Save full normalized ΔE data
output_data = [all_thickness, delta_E_plot];
output_data = sortrows(output_data, 1);
output_filename = fullfile(base_path, 'DeltaEab_vs_Thickness_Normalized.txt');

fid_out = fopen(output_filename, 'w');
fprintf(fid_out, 'Thickness_nm\tDeltaEab_Normalized\n');
fprintf(fid_out, '%.2f\t%.8f\n', output_data');
fclose(fid_out);
fprintf('Combined ΔE values saved to:\n%s\n', output_filename);

%% Save per-file % vs ΔE data
fprintf('\nSaving per-file ΔE values (normalized, percentage scale without %% sign):\n');

start_idx = 1;

for i = 1:length(file_suffixes)
    suffix = file_suffixes(i);
    filepath = fullfile(base_path, sprintf('Lab_D65_%d.txt', suffix));
    if ~isfile(filepath)
        continue;
    end

    [percentages, ~, ~, ~] = parse_nm_file(filepath);
    if i > 1
        percentages = percentages(percentages > 0);
    end
    n_points = length(percentages);

    idx_end = start_idx + n_points - 1;
    deltaE_segment = delta_E_plot(start_idx:idx_end);

    % Adjust percentages: multiply by 3, then add (suffix - 300)
    percentage_offset = suffix - 300;
    adjusted_percentages = percentages * 3 + percentage_offset;

    % Save
    outname = fullfile(base_path, sprintf('DeltaEab_Lab_D65_%d_Percent.txt', suffix));
    fid_pct = fopen(outname, 'w');
    fprintf(fid_pct, 'Percent\tDeltaEab_Normalized\n');
    for k = 1:n_points
        fprintf(fid_pct, '%.0f\t%.8f\n', adjusted_percentages(k), deltaE_segment(k));
    end
    fclose(fid_pct);

    fprintf('Saved: %s\n', outname);
    start_idx = idx_end + 1;
end

%% --- Helper Functions ---
function [L0, a0, b0] = get_reference_values(filepath)
    fid = fopen(filepath, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if contains(line, '0nm:')
            parts = split(line, ':');
            values = sscanf(strtrim(parts{2}), '%f, %f, %f');
            L0 = values(1); a0 = values(2); b0 = values(3);
            fclose(fid);
            return;
        end
    end
    fclose(fid);
    error('0nm reference not found in %s', filepath);
end
function [percents, L, a, b] = parse_nm_file(filepath)
    fid = fopen(filepath, 'r');
    percents = []; L = []; a = []; b = [];
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line) || ~contains(line, 'nm:')
            continue;
        end
        try
            parts = split(line, ':');
            nm_val = str2double(strrep(parts{1}, 'nm', ''));  % e.g., 5nm -> 5
            percent = (nm_val / 100) * 100;  % Normalize to 100%, still returns 5, 10...
            vals = sscanf(strtrim(parts{2}), '%f, %f, %f');
            if numel(vals) == 3
                percents(end+1,1) = percent;
                L(end+1,1) = vals(1);
                a(end+1,1) = vals(2);
                b(end+1,1) = vals(3);
            end
        catch
            warning('Skipping malformed line: %s', line);
        end
    end
    fclose(fid);
end