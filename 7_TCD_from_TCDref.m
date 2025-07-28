clear; clc;

%% Step 1: Load normalized ΔE vs. thickness table
data = readmatrix('DeltaEab_CIECAM_normalized.txt');
thickness_nm = data(:,1);
norm_deltaE_table = round(data(:,2), 6);  % Round to decimal places here

%% Step 2: Load RGB microscope image
[imgfile, imgpath] = uigetfile({'*.png;*.jpg;*.tif'}, 'Select RGB OM image');
if isequal(imgfile, 0), error('No image selected.'); end
img = im2double(imread(fullfile(imgpath, imgfile)));

figure, imshow(img);
title('Select rectangular region for 0 nm (bare substrate)');
ref_rect = round(getrect);  % [x, y, width, height]
ref_crop = img(ref_rect(2):(ref_rect(2)+ref_rect(4)), ...
               ref_rect(1):(ref_rect(1)+ref_rect(3)), :);
ref_RGB = squeeze(mean(reshape(ref_crop, [], 3), 1));

%% Step 3: Convert image and reference to CAM02-UCS and compute ΔE
isd = true;     % whether to normalize J'
K_L = 1;        % K_L parameter
c1  = 0.007;    % c1 parameter
c2  = 0.0228;   % c2 parameter

% Reshape image to Nx3 for batch conversion
[rows, cols, ~] = size(img);
rgb_reshaped = reshape(img, [], 3);

% Convert image and reference to CAM02-UCS
cam02_img = sRGB_to_CAM02UCS(rgb_reshaped, isd, K_L, c1, c2);  % Nx3
cam02_img = reshape(cam02_img, rows, cols, 3);

cam02_ref = sRGB_to_CAM02UCS(ref_RGB, isd, K_L, c1, c2);  % 1x3

% Compute ΔE in CAM02-UCS
delta_E = sqrt((cam02_img(:,:,1) - cam02_ref(1)).^2 + ...
               (cam02_img(:,:,2) - cam02_ref(2)).^2 + ...
               (cam02_img(:,:,3) - cam02_ref(3)).^2);

deltaE_map = imgaussfilt(delta_E, 1.2);

%% Step 5: Normalize and Display ΔE*ab map
max_deltaE = max(deltaE_map(:));
deltaE_map_norm = deltaE_map / max_deltaE;

%% Step 5.1: Ask user for a specific ΔE range to display (up to 4 decimals)
prompt = {'Enter lower ΔE bound [0–1]:', 'Enter upper ΔE bound [0–1]:'};
dlg_title = 'Select ΔE Region to Display';
dims = [1 40];
definput = {'0.2000','0.8000'};
answer = inputdlg(prompt, dlg_title, dims, definput);

if isempty(answer)
    error('No input provided for ΔE range.');
end

lower_bound = round(str2double(answer{1}), 4);
upper_bound = round(str2double(answer{2}), 4);

if isnan(lower_bound) || isnan(upper_bound) || ...
   lower_bound < 0 || upper_bound > 1 || lower_bound >= upper_bound
    error('Invalid bounds. Please enter values in [0,1] with lower < upper.');
end

% Create mask for pixels within the selected ΔE range
deltaE_mask = (deltaE_map_norm >= lower_bound) & (deltaE_map_norm <= upper_bound);

% Apply mask: set out-of-range pixels to black (0)
deltaE_filtered = deltaE_map_norm;
deltaE_filtered(~deltaE_mask) = 0;

%% Display filtered ΔE map
figure, imagesc(deltaE_filtered);
axis image off;
colormap(turbo(16));
cb = colorbar;
cb.Ticks = 0:0.0625:1;
cb.TickLabelInterpreter = 'latex';
clim([0 1]);
title(sprintf('ΔE in Range [%.4f, %.4f]', lower_bound, upper_bound), ...
      'Interpreter', 'latex');