%% RADAR MODEL
% This code illustrates a method for finding the delay between two signals.

close all;
clear all;
clc;
%% 
% The signal is constructed with given parameters.

Fs = 1e5;  % Sampling frequency
dt = 1/Fs; % Sampling period
Range = 1000; % Range
signal_length = 4e-3; % Signal length
antenna_gain = 20; % Antenna gain in dB
rcs = 0.1; % Radar cross section
Time = 0.05; % Half time axis lengh
light_speed = physconst("LightSpeed"); % Light speed
dt_div = 1000; % Division of dt for cross correlation
fc = 1e4; % Operating frequency
t = 0:dt:Time-dt; % Time scale
transmitter_power = 100; % Transmitted signal power
delay = -2*Range/light_speed; % Delay of received signal
%% 
% Transmitted sinusoidal pulse signal

window_transmitted = 1*(t>=0 & t<=signal_length); % Transmitted signal width
signal_model_t = sin(2*pi*fc*t).*window_transmitted; % Sinusoidal pulse signal for transmitter
transmitter_value_square = sum(signal_model_t.^2)/(length(nonzeros(window_transmitted))*transmitter_power);
transmitted_signal = signal_model_t/sqrt(transmitter_value_square); % Transmitted signal
%% 
% Delay, loss and noise are added and the received signal is constructed.

free_spl = (4*pi*Range*fc/light_speed)^2; % Free space path loss
antenna_gain_p = db2pow(antenna_gain); % Antenna gain in power
loss = antenna_gain_p^2*rcs/(free_spl*4*pi*Range^2); % Total loss from radar equation
tpower = sum(transmitted_signal.^2)/length(nonzeros(window_transmitted)); % Transmitted signal power
rpower = tpower*loss; % Received signal power
window_received = 1*((t+delay)>=0 & (t+delay)<=signal_length); % Received signal width
signal_model_r = sin(2*pi*fc*(t+delay)).*window_received; % Sinusoidal pulse signal for receiver
receiver_value_square = sum(signal_model_r.^2)/(length(nonzeros(window_received))*rpower); 
received_signal = signal_model_r/sqrt(receiver_value_square); % Received signal
receiver_power = sum(received_signal.^2)/length(nonzeros(window_received)); % Received signal power
received_signal_w_noise = awgn(received_signal,0,receiver_power); % Received signal with noise 
power_of_noise = sum((received_signal_w_noise-received_signal).^2)/length(t); % Power of noise 
snr = 10*log10(receiver_power/power_of_noise) % SNR value
figure;
plot(t,transmitted_signal)
figure;
plot(t,received_signal)
figure;
plot(t,received_signal_w_noise)
matched_filter = transmitted_signal(end:-1:1)';
filtered_signal = conv(received_signal_w_noise1,matched_filter);
filtered_signal = filtered_signal(round(length(filtered_signal)/4) : 3*round(length(filtered_signal)/4));
figure;
plot(t,filtered_signal,"Color","b");
title("Received Sinusoidal Pulse Signal after Matched Filtering");
xlabel("Time(s)");
ylabel("Voltage(V)");
%% 
% The signals are plotted.

figure;
plot(t,transmitted_signal,"Color","b"); 
hold on;
plot(t,received_signal_w_noise1,"Color","g");
title("Transmitted and Received Sinusoidal Pulse Signals Noise 1");
legend("Transmitted Signal","Received Signal");
xlabel("Time(s)");
ylabel("Voltage(V)");

received_signal_w_noise2 = received_signal + 0.5*randn(1,length(t)); % Received signal with noise 2
figure;
plot(t,transmitted_signal,"Color","b"); 
hold on;
plot(t,received_signal_w_noise2,"Color","g");
title("Transmitted and Received Sinusoidal Pulse Signals Noise 2");
legend("Transmitted Signal","Received Signal");
xlabel("Time(s)");
ylabel("Voltage(V)");
power_of_noise2 = sum((received_signal_w_noise2-received_signal).^2)/length(t); % Power of noise 1
%% 
% To find the fractional delay, cross correlation is used. Both auto and cross 
% correlations are computed and the highest value of the cross correlation is 
% found. Then, the index of this value is obtained, which is the delay time.

[cross_cor,t_cor] = xcorr(transmitted_signal,received_signal_w_noise1); % Cross correlation of the original and delayed signals
[auto_cor,t_acor] = xcorr(transmitted_signal); % Auto correlation of the original signal
%% 
% Both auto and cross correlations are plotted.

% Plots of the auto correlation and cross correlation functions
figure;
plot(t_cor/Fs,cross_cor,"LineWidth",1,"Color","b");
hold on;
plot(t_acor/Fs,auto_cor,"LineWidth",1,"Color","g");
title("Auto Correlation of Original Signal and Cross Correlation of Original and Delayed Signals");
legend("Cross Correlation","Auto Correlation");
xlabel("Time(s)");
ylabel("Voltage(V)");
%% 
% Peak value of the cross correlation is found with the index.

[pk,i] = max(cross_cor); % Peak value and index value of the peak for cross correlation function 
x_value = t_cor/Fs;
peak_value = pk;
unit_delay = x_value(i);
%% 
% As can be seen, this value is only the unit part of the delay. To compute 
% the fractional part, first, the original signal is delayed by the unit part. 
% Then, the cross correlation is computed in a single sampling period by dividing 
% this range to small parts and the delay is obtained in a similar way as before.

% unit = finddelay(signal,signal_del1) % Unit part of the delay
% window = 1*((t+del)>=0 & (t+del)<=8e-3); 
% signal = sin(2*pi*fc*(t+del)).*window; % Original signal is delayed according to the unit found 

% Finding the cross correlations in small interval for fractional delay
% computation
correlation_values = [];
for mini_delays = -dt:dt/dt_div:dt
    window_mini_del = 1*((t+unit_delay+mini_delays)>=0 & (t+unit_delay+mini_delays)<=signal_length);
    signal_mini_del = sin(2*pi*fc*(t+unit_delay+mini_delays)).*window_mini_del;
    cor_value = sum(signal_mini_del.*received_signal_w_noise1);
    correlation_values = [correlation_values;cor_value];
end
mini_delays = -dt:dt/dt_div:dt;

% Plot of the cross correlation values
figure;
plot(mini_delays,correlation_values,"LineWidth",1,"Color","r");
[pk,i] = max(correlation_values); % Peak value of the cross correlation and its index
fractional_delay = mini_delays(i); % Fractional delay
total_delay = unit_delay + fractional_delay;
range_of_target = light_speed*abs(total_delay)/2
