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

% Define the porosity, tortuosity, and constrictivity
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

        % Compute the effective diffusion coefficient using the given equation
        D_eff = D_free * (porosity / tortuosity) * constrictivity;

        % Calculate time taken for bacteria to pass through porous medium
        time_diffusion = (L^2) / (2 * D_eff);
        time_advection = L / adjusted_velocity;
            
        % Include running and tumbling times
        running_time = 1.25; % sec
        tumbling_time = 0.17; % sec
            
        time_total = time_diffusion + time_advection + running_time - tumbling_time;

        % Define the number of time steps for simulation
        number_of_time_steps = 1000; % Adjust as needed
            
        % Initialize the concentration profile
        concentration = zeros(N+1, 1);
        % Set the initial condition (for example, a Dirac delta function at the center)
        concentration(round(N/2)) = 1;

        % Time step and number of time steps
        dt = time_total / number_of_time_steps;
        num_steps = round(time_total / dt);

        % Calculate the concentration profile using the finite difference scheme
        for step = 1:num_steps
            new_concentration = zeros(N+1, 1);  % Create a new array to store updated concentrations

            % Update concentrations using the finite difference formula
            for i = 2:N
                dCdx2 = (concentration(i+1) - 2 * concentration(i) + concentration(i-1)) / (dx^2);
                new_concentration(i) = concentration(i) + dt * D_eff * dCdx2;
            end

            % Update concentration due to growth
            growth_factor = exp(growth_rate * dt);
            new_concentration = new_concentration * growth_factor;

            % Set updated concentrations as the new concentration profile
            concentration = new_concentration;
        end

        % Store the results in the matrices
        total_time_matrix(p_idx, temp_idx) = time_total;
        concentration_matrix(:, p_idx, temp_idx) = concentration;
    end
end

% Calculate the average total time
average_total_time = mean(total_time_matrix(:));

% Display the average total time
fprintf('The average total time taken by E. coli to pass through the porous medium is %.2f seconds.\n', average_total_time);

% Create a 3D surface plot for time taken by E. coli
[X, Y] = meshgrid(pore_size_range * 1e9, temp_range_C);
figure('Color', 'white'); % Setting the background to white
h = surf(X, Y, squeeze(total_time_matrix(:, :, 1)), 'EdgeColor', 'none', 'FaceAlpha', 0.8);

% Improve lighting
light('Position', [-1 0 1], 'Style', 'infinite');
shading interp;
material dull;

% Add labels and title
xlabel('Pore Size (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Temperature (C)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Average Total Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Average Total Time taken by E. coli to Pass Through 20\mu\rm{m} Plant Root', 'FontSize', 14, 'FontWeight', 'bold');

% Add a color bar
cbar = colorbar;
ylabel(cbar, 'Total Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
% Set colormap
colormap(parula(256));

% Adjust view angle for better perspective
view(3);
% Create a 3D surface plot for time taken by E. coli
[X, Y] = meshgrid(pore_size_range * 1e9, pore_length_range * 1e9);
figure;
surf(X, Y, squeeze(total_time_matrix(:, :, 1))); % Choose a specific temperature slice

% Add labels and title
xlabel('Pore Size (nm)');
ylabel('Temperature (C)');
zlabel('Average Total Time (s)');
title('Average Total Time taken by E. coli to Pass Through 20\mu\rm{m} Plant Root');

% Customize axis and grid
axis tight;
grid on;
set(gca, 'GridAlpha', 0.3, 'Box', 'off', 'FontSize', 12, 'LineWidth', 1.2);

% Add a color bar
cbar = colorbar;
ylabel(cbar, 'Total Time (s)', 'FontSize', 12, 'FontWeight', 'bold');


