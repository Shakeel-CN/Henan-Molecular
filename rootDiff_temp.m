% Define the domain and initial conditions
L = 0.002; % Length of domain (cm)
N = 100; % Number of grid points
dx = L/N; % Grid spacing
x = linspace(0, L, N+1); % Grid points

% Varying parameters
velocity = 20e-6; % Velocity of E.coli in m/s
pore_size_range = linspace(3.3e-9, 6.2e-9, 100); % Pore size range in m
temp_range_C = linspace(30, 50, 100); % Temperature range in Celsius
temp_range_K = temp_range_C + 273.15; % Convert temperature to Kelvin

% Preallocate matrices for total time and concentration
total_time_matrix = zeros(length(pore_size_range), length(temp_range_K));
concentration_matrix = zeros(N+1, length(pore_size_range), length(temp_range_K));

% Free water diffusion coefficient for E. coli at 25 degrees Celsius
D_free = 200e-12; % m^2/s

% Define the porosity, tortuosity and constrictivity
porosity = 0.35; % Assumed constant
tortuosity = 2.11; % Assumed constant
constrictivity = 1; % Assumed constant

% Loop through all parameter combinations
for p_idx = 1:length(pore_size_range)
    for temp_idx = 1:length(temp_range_K)
        pore_size = pore_size_range(p_idx);
        T = temp_range_K(temp_idx);
        temp_effect = exp(-(T - min(temp_range_K))/(max(temp_range_K) - min(temp_range_K))); % Temperature effect

        % Calculate the specific growth rate of E. coli
        growth_rate = 0.75/3600; % Set Avg growth rate to 0.75 h^-1 and change from hr to sec

        % Calculate the adjusted velocity based on pore size and temperature
        adjusted_velocity = velocity * (pore_size / max(pore_size_range)) * temp_effect;

        % Compute the effective diffusion coefficient using given equation
        D_eff = D_free * (porosity / tortuosity) * constrictivity;

        % Calculate time taken for bacteria to pass through porous medium
        time_diffusion = (L^2) / (2 * D_eff);
        time_advection = L / adjusted_velocity;
        
        % Include running and tumbling times
        running_time = 1.25; % sec
        tumbling_time = 0.17; % sec
        
        time_total = time_diffusion + time_advection + running_time - tumbling_time;

        % Store the results in the matrices
        total_time_matrix(p_idx, temp_idx) = time_total;
    end
end

% Calculate the average total time
average_total_time = mean(total_time_matrix(:));

% Display the average total time
fprintf('The average total time taken by E. coli to pass through the porous medium is %.2f seconds.\n', average_total_time);

% Create a 2D surface plot for time taken by E. coli w.r.t. temperature and pore size
[X, Y] = meshgrid(pore_size_range * 1e9, temp_range_C);
figure('Color', 'white'); % Setting the background to white
surf(X, Y, squeeze(total_time_matrix(:, :)), 'EdgeColor', 'none');
colormap(parula);
cbar = colorbar;
ylabel(cbar, 'Total Time (s)', 'FontSize', 12, 'FontWeight', 'bold');

% Add labels and title
xlabel('Pore Size (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (\circC)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Average Total Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Impact of Temperature on E. coli Passage Through Porous Medium', 'FontSize', 14, 'FontWeight', 'bold');

% Customize axis and grid
axis tight;
grid on;
set(gca, 'GridAlpha', 0.3, 'Box', 'off', 'FontSize', 12, 'LineWidth', 1.2);

% Create a 2D surface plot for time taken by E. coli w.r.t. temperature and pore size
[X, Y] = meshgrid(pore_size_range * 1e9, temp_range_C);

% Save the data to an Excel file
excel_filename = 'total_time_matrix.xlsx';

% Prepare data for Excel
excel_data = [pore_size_range' * 1e9, total_time_matrix];

% Create headers for the data
header_pore_size = 'Pore Size (nm)';
header_time = arrayfun(@(x) sprintf('Time (s) at Temp %d C', x), temp_range_C, 'UniformOutput', false);
header = [header_pore_size, header_time];

% Write data to Excel
xlswrite(excel_filename, header, 'Sheet1', 'A1');
xlswrite(excel_filename, excel_data, 'Sheet1', 'A2');

fprintf('Excel file "%s" has been created.\n', excel_filename);
% Create a 2D surface plot for time taken by E. coli w.r.t. temperature and pore size
[X, Y] = meshgrid(pore_size_range * 1e9, temp_range_C);

% Save the data to an Excel file
excel_filename = 'total_time_pore.xlsx';

% Prepare data for Excel
excel_data = [pore_size_range' * 1e9, total_time_matrix];

% Create a header for the data
header = ['Pore Size (nm)', arrayfun(@(x) sprintf('Temp %d C', x), temp_range_C, 'UniformOutput', false)];

% Write data to Excel
xlswrite(excel_filename, header, 'Sheet1', 'A1');
xlswrite(excel_filename, excel_data, 'Sheet1', 'A2');

fprintf('Excel file "%s" has been created.\n', excel_filename);
