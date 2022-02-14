%%
% Load EMG signals
file_name = {'emg_data.mat', 'EMGdata.txt'};
x = load(file_name {1});
fs = 1000;
% fs = x.fs;

%%
% EMGdata.txt is organized into three columns: 
% torque, bicepsEMG and tricepsEMG
% torque = x(:, 1);
% biceps = x(:, 2);
% triceps = x(:, 3);

% signal = x.emg_signali{4}';
% signal = biceps;
signal = round(20 * rand(200, 10));

%%
free_params = [];
n_quant = 10;
scaling_fact = 4;
model_type = 'laplace';

%%
% Obtain the MAP estimation of the EMG signal
map_x = bayes_filt(signal, fs, model_type, free_params, n_quant, scaling_fact);

%%
% Plot normalized MAP estimation
Ts = 1 / fs;
t = 0:Ts:length(signal) / fs - Ts; 
plot(t, [map_x / max(map_x), abs(signal) / max(abs(signal))])
