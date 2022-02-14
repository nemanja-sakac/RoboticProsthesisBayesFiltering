function p_signal_x = laplace_model(signal_val, x)
% LAPLACE_MODEL Models signal as a Laplacian probability density function
%
% p_signal_x = LAPLACE_MODEL(signal_val, x) returns the estimated 
% probability density p(signal|x) for a current time sample t, approximated
% by a Laplace probability density function. 'signal_val' is the signal 
% value at time sample t and 'x' is the underlying "driving process" that 
% generated the signal.

p_signal_x = zeros(length(x), length(signal_val));
for i = 1:length(signal_val)
    p_signal_x(:, i) = exp(-signal_val(i) ./ x) ./ x;
end