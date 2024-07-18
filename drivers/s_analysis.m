function T = s_analysis(pin, m, Xt)

Tend = pin.Tend;
damping_type = pin.damping_type;
d = pin.d;
mu = pin.mu;
spectrum_type = pin.spectrum_type;
E0 = pin.E0;
num_modes_flat = pin.num_modes_flat;
sigmax = pin.sigmax;
N = pin.computed.N;
L = 2*pi;

% Basic Statistical Analysis
mean_m = mean(m, 2);  % Mean along rows
median_m = median(m, 2);  % Median along rows
std_m = std(m, 0, 2);  % Standard Deviation along rows
var_m = var(m, 0, 2);  % Variance along rows
skew_m = skewness(m, 0, 2);  % Skewness along rows
kurt_m = kurtosis(m, 0, 2);  % Kurtosis along rows

mean_Xt = mean(Xt, 2);  % Mean along rows

% Create a table with the statistics
T = table(mean_m, median_m, std_m, var_m, skew_m, kurt_m, mean_Xt);

% Display the table
disp(T); 

% Plotting
figure;
subplot(2,2,1);
histogram(m);
title('Histogram of m');

subplot(2,2,2);
boxplot(m);
title('Box Plot of m');

subplot(2,2,3);
plot(mean_m);
title('Time Series of Mean Values');

subplot(2,2,4);
scatter(mean_m, var_m);
title('Mean vs. Variance');

% Hypothesis Testing (Example: Test if mean of m is zero)
[h, p] = ttest(m);
disp(['t-test p-value: ', num2str(p)]);

end
