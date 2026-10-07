% Define file paths
txt_file = 'PS_beads.txt';       % Replace with your .txt file
xls_file = 'Index_of_Refraction_library.xls';       % Replace with your .xls file
output_file = 'updated_file.xls'; % Output Excel file

% --- Read the .txt file (wavelength, n, k) ---
txt_data = readmatrix(txt_file);
wavelength_txt = txt_data(:, 1)*1000;
n_txt = txt_data(:, 2);
k_txt = txt_data(:, 3);

% --- Interpolate over 300–800 nm ---
interp_range = 300:1:800;
n_interp = interp1(wavelength_txt, n_txt, interp_range, 'linear', 'extrap');
k_interp = interp1(wavelength_txt, k_txt, interp_range, 'linear', 'extrap');

% --- Read the Excel file ---
raw_data = readcell(xls_file);              % Returns cell array with headers
wavelength_xls = cell2mat(raw_data(2:end, 1));  % Column 1, skip header

% --- Interpolate onto Excel's wavelength grid ---
n_new = interp1(interp_range, n_interp, wavelength_xls, 'linear', 'extrap');
k_new = interp1(interp_range, k_interp, wavelength_xls, 'linear', 'extrap');

% --- Append column headers at AH1 and AI1 ---
raw_data{1, 36} = 'PS_beads_n';   % Column AJ - 36
raw_data{1, 37} = 'PS_beads_k';   % Column AK - 37

% --- Append interpolated values ---
for i = 1:length(wavelength_xls)
    raw_data{i+1, 36} = n_new(i);  
    raw_data{i+1, 37} = k_new(i);  
end

% --- Write back to Excel ---
writecell(raw_data, output_file);
disp(['Done! File saved as: ', output_file]);
