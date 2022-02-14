%%
clear
clc
addpath('../src');

%%
% Load EMG signals
filename = '../bin/emg_data.mat';
x = load(filename);
fs = x.fs;

%%
% Bayes filter parameter settings
free_params = [];
n_quant = 50;
scaling_fact = 4;
model_type = 'laplace';

%%
% min_length = length(x.emg_signali{1});
% for i = 1:length(x.emg_signali)
%     if (length(x.emg_signali{i}) < min_length)
%         min_length = length(x.emg_signali{i});
%     end
% end
% 
% multichannel_sig = zeros(min_length, length(x.emg_signali));
% for i = 1:length(x.emg_signali)
%     multichannel_sig(:, i) = (x.emg_signali{i}(1:min_length))';
% end

%%
% Test performance for each signal
% Sampling period and frequency
Ts = 1 / fs;
% Time vector
meas_time = zeros(length(x.emg_signali), 1);

for i = 1:length(x.emg_signali)
    signal = x.emg_signali{i}';
    
    % Signal filtering, envelope calculation, performance measurement
    tic
    map_x = bayes_filt(signal, fs, model_type, free_params, n_quant,...
        scaling_fact);
    meas_time(i) = toc;
    
    % Normalize values for comparison
    signal = abs(signal);
    for j = 1:size(signal, 2)
        signal(:, j) = signal(:, j) / max(signal(:, j));
        map_x(:, j) = (map_x(:, j) - min(map_x(:, j))) / max(map_x(:, j));
    end
    
    t = 0:Ts:size(signal, 1) / fs - Ts;
    
    % Display signal and envelope
    figure
    plot(t, [map_x, signal])
    xlabel('t(s)')
    legend('MAP estimate', 'EMG')
    annotation('textbox', [0.15, 0.8, 0.25, 0.1],...
        'String', ['Execution time: ', num2str(meas_time(i)), ' s']);
    
    % Save figure
    filename = ['../img/emg', num2str(i), '.png'];
    saveas(gcf, filename);
end

%%
% Test performance for "multichannel" signal
% tic
% map_x = bayes_filt(multichannel_sig, fs, model_type, free_params, n_quant,...
%     scaling_fact);
% meas_time = toc;
% 
% multichannel_sig = abs(multichannel_sig);
% for j = 1:size(multichannel_sig, 2)
%     multichannel_sig(:, j) = multichannel_sig(:, j)...
%         / max(multichannel_sig(:, j));
%     map_x(:, j) = (map_x(:, j)) / max(map_x(:, j));
% end
% 
% % - min(map_x(:, j))
% for i = 1:size(multichannel_sig, 2)
%     
%     t = 0:Ts:size(multichannel_sig, 1) / fs - Ts;
%     
%     % Display signal and envelope
%     figure
%     plot(t, [map_x(:, i), multichannel_sig(:, i)])
%     xlabel('t(s)')
%     legend('MAP estimate', 'EMG')
% %     annotation('textbox', [0.15, 0.8, 0.25, 0.1],...
% %         'String', ['Execution time: ', num2str(meas_time(i)), ' s']);
% 
%     % Save figures
%     filename = ['../img/emg', num2str(i), '.png'];
% %     saveas(gcf, filename);
% end
%%
close all
