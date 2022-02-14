function map_x = bayes_filt(signal, fs, model_type, varargin)
% BAYES_FILT Nonlinear Bayesian filtering developed for the estimation of
% an underlying "driving process" of a signal. Based on paper and algorithm 
% demonstration by Sanger [1].
%
%   map_x = BAYES_FILT(signal, fs, model_type, ...) - returns the 
%   envelope of the input signal, given by parameter 'signal'. 'fs' is the 
%   sampling frequency of the input signal. 'model_type' is the specified 
%   model of the signal based on an underlying "driving process". Two types
%   of models have been implemented:
%
%       'gauss' - Models the signal as amplitude-modulated zero-mean
%       Gaussian noise.
%
%       'laplace' - Models the signal with a Laplace distribution.
%       
%   Parameters subsequently described are optional. In the case of not 
%   being defined by the user, optional parameters will take on default
%   values. The first optional parameter is a vector of three 'free 
%   parameters', their meanings defined as follows:
%
%       alpha - Expected rate of gradual drift in the signal (signal units 
%       / unit time).
%
%       beta - The expected rate of sudden shifts in signal value (number 
%       of shifts / unit time).
%
%       gamma - Probability of error per single measurement (measurement 
%       uncertainty).
% 
%   Current implementation requires all three parameters to be defined. 
%   The parameters in the vector should be specified in the following order: 
%       free_params(1) = alpha, 
%       free_params(2) = beta, 
%       free_params(3) = gamma.
%   Alternatively, an empty vector can be defined, in which case default 
%   values for the free parameters are used. 
%       
%   'n_outputs' is the number of "quantization levels" (possible values) 
%   for the estimation of posterior probability functions. 'scaling_fact' 
%   is the scaling factor by which the rectified signal is amplified before
%   its amplitude is cut off to a maximum value of 1. This ensures that the
%   variability of the signal doesn't affect the overall estimation.
%
%
%   Author: Nemanja Sakac, Faculty of Technical Sciences, Novi Sad, 2021.
%   References:
%       [1] Sanger, T.D., "Bayesian Filtering of Myoelectric Signals", 
%           Journal of Neurophysiology, vol. 97, pp. 1839-1845, 2007.


% Default parameter definitions
% Maximum (cut-off) value for the rectified signal 
rectified_max = 1;
% Sampling period
delta_t = 1 / fs;
% Mapping default values to variable names
keySet = {'alpha', 'beta', 'gamma', 'n_outputs', 'scaling_fact'};
valueSet = [(0.001 * delta_t), (10E-24 * delta_t), 0, 50, 4];
defaults = containers.Map(keySet, valueSet);

% Setting the optional parameters (free_params, n_outputs, scaling_fact)
switch(max(size(varargin)))
    case 0
        alpha = defaults('alpha');
        beta = defaults('beta');
        gamma = defaults('gamma');
        n_outputs = defaults('n_outputs');
        scaling_fact = defaults('scaling_fact');
    case 1
        if isempty(varargin{1})
            alpha = defaults('alpha');
            beta = defaults('beta');
            gamma = defaults('gamma');
        else     
            alpha = varargin{1}(1);
            beta = varargin{1}(2);
            gamma = varargin{1}(3);
        end
        n_outputs = defaults('n_outputs');
        scaling_fact = defaults('scaling_fact');
    case 2
        if isempty(varargin{1})
            alpha = defaults('alpha');
            beta = defaults('beta');
            gamma = defaults('gamma');
        else     
            alpha = varargin{1}(1);
            beta = varargin{1}(2);
            gamma = varargin{1}(3);
        end
        n_outputs = varargin{2};
        scaling_fact = defaults('scaling_fact');
    case 3
        if isempty(varargin{1})
            alpha = defaults('alpha');
            beta = defaults('beta');
            gamma = defaults('gamma');
        else     
            alpha = varargin{1}(1);
            beta = varargin{1}(2);
            gamma = varargin{1}(3);
        end
        n_outputs = varargin{2};
        scaling_fact = varargin{3};
end

% Obtaining the handle for the specified stochastic process model given by
% parameter 'model_type'
switch(model_type)
    case 'gauss'
        spec_model = @gauss_model;
    case 'laplace'
        spec_model = @laplace_model;
    otherwise
        spec_model = @laplace_model;
end


% Signal preprocessing
% Mean of the signal is removed and the signal is rectified
mean_vals = mean(signal);
signal = abs(signal);
for i = 1:size(signal, 2)
    signal(:, i) = signal(:, i) - mean_vals(i);
    % Prescaling the rectified signal 
    signal(:, i) = scaling_fact * rectified_max * signal(:, i)...
        / max(signal(:, i));
end
% Cutting off amplitudes of the scaled rectified signal if they are above
% the preset maximum value
signal(signal > rectified_max) = rectified_max;

% Initialization of the underlying driving process for the generation of
% the signal
x = linspace(rectified_max / n_outputs, rectified_max, n_outputs)';
% Initialization of the estimated MAP(x(t))
map_x = zeros(size(signal));
% Approximation of the second derivative with second differences
dp_x_t = [alpha, (1 - 2 * alpha), alpha];


% Implementation of the algorithm according to the steps described in
% Sanger [1].

% 1. Initialization of the prior probability p(x,t) as a 
% uniform distribution. Each column corresponds to a different channel.
prior_x = ones(n_outputs, size(signal, 2)) / n_outputs;

% The envelope is estimated for each sample of the signal
for t = 1:size(signal, 1)
    
    % 2. Forward propagate p(x,t-) using the convolution function
    % NOTE: The original implementation by the author of the paper uses
    % filtfilt, while conv is used here for an overall faster running time.
    % Set probability of a sudden jump (step) value
    beta_prior = beta + (1 - beta) * prior_x;
    % Calculates the second derivative
    for i = 1:size(prior_x, 2)
        prior_x(:, i) = conv(prior_x(:, i), dp_x_t, 'same');
    end
    prior_x = prior_x + beta_prior;
    
    % 3. Evaluation of the signal value at current sample
    signal_val = signal(t, :);
    
    % 4. Calculate the posterior likelihood function
    % Calculate the probability density p(signal|x) according to the
    % specified model.
    p_signal_x = spec_model(signal_val, x);
    % Calculate posterior density using Bayes rule
    posterior_x = (1 - gamma) * p_signal_x .* prior_x + gamma; 
    
    % 5. Output the signal estimate MAP(x(t)) = argmax p(x,t)
    % Finding the value for which probability is maximum.
    % min ensures that only the first value is selected in case of
    % multiple maximum value samples
    peak_x = zeros(size(posterior_x, 2), 1);
    for i = 1:length(peak_x)
        peak_x(i) = min(find(posterior_x(:, i) == max(posterior_x(:, i))));
    end
    
    % Interpolating the peak for a better estimation
%     peak_index = peak_x;
    peak_index = zeros(length(peak_x), 1); 
    for i = 1:length(peak_x)
        if (peak_x(i) > 1 && peak_x(i) < n_outputs)
            % Left difference
            dL = posterior_x(peak_x(i) - 1, i) - posterior_x(peak_x(i), i);  
            % Right difference
            dR = posterior_x(peak_x(i), i) - posterior_x(peak_x(i) + 1, i);
            peak_index(i) = (peak_x(i) - 0.5 - (dL / (dR - dL)));
        % If maximum occurs at an endpoint do not interpolate
        else
            peak_index = peak_x(i);    
        end
    end
    % Convert peak index to scaled signal
    map_x(t, :) = (rectified_max / (n_outputs - 1)) * peak_index;
    
    % 6. Rescaling the posterior density so sum(posterior_x) = 1
    sum_post_x = sum(posterior_x);
    for i = 1:size(posterior_x, 2)
        posterior_x(:, i) = posterior_x(:, i) / sum_post_x(i);
    end

    % 7. Repeat from step 2
    % Prior for next iteration is posterior from this iteration
    prior_x = posterior_x;
    
end
