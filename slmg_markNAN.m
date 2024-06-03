function  [xnum_nans, ynum_nans] = slmg_markNAN(x, y)
% Create a figure
figure;

% Plot the numbers, ignoring NaNs for the main plot
subplot(2,1,1);
plot(x, 'b'); % 'b-o' specifies a blue line with circle markers
hold on;
xnan_indices = find(isnan(x)); % Identify NaN values
% Plot 'x' at y=0 for NaN values
plot(xnan_indices, zeros(size(xnan_indices)), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % 'rx' specifies red crosses
title('X Coordinates');
legend('Original',  'NaN Values Marked at y=0');

subplot(2,1,2);
plot(y, 'b'); % 'b-o' specifies a blue line with circle markers
hold on;
ynan_indices = find(isnan(y));% Identify NaN values
% Plot 'x' at y=0 for NaN values
plot(ynan_indices, zeros(size(ynan_indices)), 'rx', 'MarkerSize', 10, 'LineWidth', 2); % 'rx' specifies red crosses
title('Y Coordinates');
legend('Original',  'NAN');

% Turn on grid
grid on;

% Display the plot
hold off;

% Count the number of NaN values
xnum_nans = sum(isnan(x));
ynum_nans = sum(isnan(y));

% Display the result
fprintf('The number of NaN values in x_clean is: %d\n', xnum_nans);
fprintf('The number of NaN values in y_clean is: %d\n', ynum_nans);

end