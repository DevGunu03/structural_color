% === GUI FILE SELECTION ===
[filename, folder_path] = uigetfile('*.txt', 'Select an sRGB .txt file to convert');
if isequal(filename, 0)
    error('No file selected. Operation cancelled.');
end

file_pattern = filename;

% === LIST SELECTED FILE ONLY (for future compatibility with wildcards) ===
files = dir(fullfile(folder_path, file_pattern));

% === PROCESS EACH FILE ===
for k = 1:length(files)
    filename = files(k).name;
    filepath = fullfile(folder_path, filename);

    % === READ LINES ===
    try
        lines = readlines(filepath);
    catch
        warning('Failed to read file: %s', filename);
        continue;
    end

    % === SKIP HEADER ===
    lines = lines(2:end);

    % === PREP OUTPUT FILE ===
    output_filename = fullfile(folder_path, ['CAM02UCS_' filename]);
    fid = fopen(output_filename, 'w');
    fprintf(fid, 'Thickness (nm): J′, a′, b′\n');

    % === PARSE AND CONVERT EACH LINE ===
    for i = 1:numel(lines)
        line = strtrim(lines(i));
        tokens = regexp(line, '(\d+)nm:\s*([\d\.Ee+-]+),\s*([\d\.Ee+-]+),\s*([\d\.Ee+-]+)', 'tokens');

        if isempty(tokens)
            warning('Skipping malformed line: %s', line);
            continue;
        end

        % Extract values
        thickness_str = tokens{1}{1};
        rgb = str2double(tokens{1}(2:4));

        % Convert to CAM02-UCS
        isd = true;                % normalize J' for deltaE
        K_L = 1;                   % K_L parameter
        c1  = 0.007;               % c1 parameter
        c2  = 0.0228;              % c2 parameter
        
        cam02 = sRGB_to_CAM02UCS(rgb, isd, K_L, c1, c2);

        % Write to output file
        fprintf(fid, '%snm: %.8f, %.8f, %.8f\n', thickness_str, cam02(1), cam02(2), cam02(3));
    end

    fclose(fid);
    fprintf('Processed: %s → %s\n', filename, output_filename);
end

disp('All files processed.');
