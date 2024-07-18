% Load the X data set
Ltracers = 10;
sigmax = 0.1;
obsdatafile = dir(fullfile(sprintf('obs_Ltracers%i_sigmax%.2g.mat', Ltracers, sigmax)));
if isempty(obsdatafile)
    error('Data file not found.');
end
obsdata = load(fullfile(obsdatafile.folder, obsdatafile.name));
X_data = cell2mat(obsdata.X);

% Define the number of elements for x and y velocities
num_elements = size(X_data, 1) / 2;  % This assumes the data is structured as 2*100001 rows

% Define spatial coordinates assuming uniform spacing (adjust as needed)
x = linspace(0, 1, num_elements);  % Spatial coordinate for x
y = linspace(0, 1, num_elements);  % Spatial coordinate for y

% Expand x and y to full grid matrices
[X, Y] = meshgrid(x, y);

% Initialize matrices for velocity and stream function
velocity_field_x = zeros(num_elements, size(X_data, 2));
velocity_field_y = zeros(num_elements, size(X_data, 2));
psi_field = zeros(num_elements, size(X_data, 2));  % Stream function matrix

% Loop through each column (assuming each column is a separate set of observations)
for i = 1:size(X_data, 2)
    % Extract the x and y velocity components
    vx = X_data(1:num_elements, i);
    vy = X_data(num_elements+1:end, i);

    % Integration to find stream function ψ from vy assuming ψ(x=0, y) = 0
    psi = cumsum(vy) * (y(2) - y(1));  % Discrete integration of vy w.r.t y

    % Adjust stream function to zero at one boundary if necessary (here, at y=0)
    psi = psi - psi(1);  % Making sure ψ starts from zero at the boundary

    % Calculate vx and vy from stream function
    diff_psi_dy = diff(psi) / diff(y(1:2));
    vx_new = [0; diff_psi_dy];  % Prepending zero to restore dimension
    diff_psi_dx = diff(psi) / diff(x(1:2));
    vy_new = -[0; diff_psi_dx];  % Prepending zero and applying negative sign

    % Store results
    velocity_field_x(:, i) = vx_new;
    velocity_field_y(:, i) = vy_new;
    psi_field(:, i) = psi;
end

% Visualization of the first velocity field
figure;
quiver(X, Y, velocity_field_x(:, 1), velocity_field_y(:, 1));
title('Simulated Velocity Field from ψ');
xlabel('X');
ylabel('Y');
