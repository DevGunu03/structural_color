clear; clc;

% Ask user for final thickness of PS
start_thickness = input('Enter starting thickness for PS (in nm): ');
gap = input('Enter gap for PS (in nm): ');
max_thickness = input('Enter maximum thickness for PS (in nm): ');

% Define wavelength range
lambda_range = 360:830;

% Base layers and fixed thicknesses
layers = {'Air', 'PS-beads', 'SiO2-Franta', 'Si-Franta'};
thicknesses_base = [0, 0, 100, 0.43e6];  % PS starts at 0 nm
incl = 1; % For using Maxwell-Garnett equation or the EMA approximation

% Loop over thickness from 0 to max_thickness (inclusive)
for t_val = start_thickness:gap:max_thickness

    % Update thickness of aMoO3
    thicknesses = thicknesses_base;
    thicknesses(2) = t_val;
    
    % Call TransferMatrix_Updated
    Reflectance = TransferMatrix_Updated_multiple(layers, thicknesses, lambda_range, incl); 

    % Prepare data to save
    output_data = [lambda_range(:), Reflectance(:)];
    filename = sprintf('PS_%dnm.txt', t_val);
    writematrix(output_data, filename);
    fprintf('Saved reflectance to %s\n', filename);
end
