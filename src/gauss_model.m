function p_signal_x = gauss_model(signal_val, x)
% GAUSS_MODEL Models signal as amplitude-modulated zero-mean Gaussian
% noise.
%
% p_signal_x = GAUSS_MODEL(signal_val, x) returns the estimated probability
% density p(signal|x) for a current time sample t, based on a Gaussian 
% process model. 'signal_val' is the signal value at time sample t and 'x' 
% is the underlying "driving process" that generated the signal.

p_signal_x = 2 * exp(-signal_val ^ 2 ./ (2 * x .^ 2))...
    ./ sqrt(2 * pi * x .^ 2);