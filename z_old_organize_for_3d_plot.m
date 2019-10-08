% Function to organize generated data for plotting

% Input descriptions
% x_in | pH values calculated during calc_speciation
% y_in | phantom network plateau moduli calculated in calc_Gp_linear_Miller_Macosko
% x_min | minimum pH we want to consider (must have data down to this pH range for each metal concentration)
% x_max | maximum pH we want to consider (must have data up to this pH range for each metal concentration)
% x_inc | the increment or step size we want to discretize our pH data

function [x_out, y_out] = organize_for_3d_plot(x_in, y_in, x_min, x_max, x_inc)

% number of inputed pH values
x_num = (x_max - x_min) ./ x_inc + 1;

% discretized pH values
target_x = linspace(x_min,x_max,x_num)';

% initialize vector for interpolated phantom network moduli
target_y = zeros(1,length(target_x))';

for target_number = 1:length(target_x)
    
    % pH we want to interpolate
    x_target = target_x(target_number);
    
    % difference between inputed pH values and the our target
    x_subtraction = x_in - x_target;
    
    % inputed pH indices that are right above and below the target
    index_low = find(x_subtraction > 0, 1, 'first') - 1;
    index_high = index_low + 1;
    
    % get the calculated pH and moduli values that are right above and below the target pH value
    x_low = x_in(index_low);
    x_high = x_in(index_high);
    y_low = y_in(index_low);
    y_high = y_in(index_high);
    
    % how far from x_low is the traget (0 is at x_target = x_low, 1 is at x_target = x_high)
    x_error = (x_target - x_low) ./ (x_high - x_low);
    
    % interpolate moduli value
    y_target = abs(y_high - y_low) .* x_error + y_low;
    
    % if moduli is less than 10 Pa (0.01kPa), set it to 0
    if y_target < 10^1
        y_target = 0;
    end
    
    % record corresponding discretized moduli
    target_y(target_number) = y_target;
    
end

% output discretized pH
x_out = target_x;
% output discretized moduli
y_out = target_y;

end
