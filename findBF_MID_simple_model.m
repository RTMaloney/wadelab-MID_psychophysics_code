
function [sumsquarederror] = findBF_MID_simple_model(params, maxparams, data, inv_alpha)

% Find Mean Absolute Deviation between data and predictions
%Generate predictions
w_CD = params(1);
w_IOVD = params(2);

predictions = w_CD .* inv_alpha(:,1) + w_IOVD .* inv_alpha(:,2);

%meanabsdev = mean(abs(data - predictions)); %Dave's method (produces more or less same result)
sumsquarederror = sum((data - predictions).^2); %minimise the sum-squared error
% Artificially adjust to account for parameters above and below specified
% min and max (assumes min allowable parameter value is 0 for now)
if any([min(params) < 0, min(maxparams - params) < 0])
disp('Out of Bounds')
    sumsquarederror = realmax;
end