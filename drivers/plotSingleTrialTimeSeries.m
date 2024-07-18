function plotSingleTrialTimeSeries(XT, m, trial_index, time_point)
 % Assume your data is stored in 'XT' and 'm' matrices with dimensions (20, 100000)
% Choose any index that you would like to compare (e.g., idx = 5)
idx = 5;

% Extract selected rows
selected_row_XT = XT(idx, :);
selected_row_m = m(idx, :);

% Define the start date and time for the x-axis - adjust accordingly
start_datetime = datetime('2023-03-28 00:00:00'); % YYYY-MM-DD hh:mm:ss

% Calculate the end date based on the desired length of the plot
end_datetime = start_datetime + hours(0) + minutes(0) + seconds(2000);

% Generate a datetime vector with 1 second intervals covering 2000 seconds
time_vector = seconds(linspace(0, 2000, numel(selected_row_XT)));
time_array = datenum(start_datetime) + time_vector;

% Convert timenums to datetime objects
time_array = datetime(time_array, 'ConvertFrom', 'datenum');

% Initialize the figure
figure;
title('Comparison Between Truth XT and Predictive Lago (Filter Recovery)');
xlabel('Time (Seconds)');
ylabel('Value');
grid on;

% Plot the data
hold on;
plot(time_array, selected_row_XT, '-b', 'LineWidth', 1.5);
hold on;
stairs(time_array, selected_row_m, 'r', 'LineStyle', '--', 'LineWidth', 1.5);

% Legend
legend('Truth XT', 'Predictive Lago', 'Location', 'best');
hold off;
end
