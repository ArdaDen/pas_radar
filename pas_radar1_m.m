%% RADAR MODEL
% This code illustrates a method for finding the delay between two signals.

close all;
clear all;
clc;
%% 
% The signal is constructed with given parameters.

Fs = 1e8;  % Sampling frequency
dt = 1/Fs; % Sampling period
Range = 5837.463; % Range
min_range = 200; % Minimum range
max_range = 10000; % Maximum range
antenna_gain = 20; % Antenna gain in dB
rcs = 0.1; % Radar cross section
light_speed = physconst("LightSpeed"); % Light speed
signal_length = 2*min_range/light_speed; % Signal length
Time = 2*max_range/light_speed; % PRI
prf = 1/Time;
snr = 3;
dt_div = 100; % Division of dt for cross correlation
fc = 1e7; % Operating frequency
transmitter_power_per_pulse = 100; % Transmitted signal power
pulse_number = 10; % Number of pulses 
t = 0:dt:Time;
rad_velocity_in = 200; % Velocity of the target m/s
rad_acc_in = 10; % Accelaration of the target m/s^2
antenna_gain_p = 10*log10(antenna_gain); % Antenna gain in linear
%% 
% Transmitted sinusoidal pulse signaln

transmitter_data = zeros(length(t),pulse_number);
for i = 0:pulse_number-1
    window_transmitted_i = 1*(t>=0 & t<=signal_length);
    signal_model_t_i = sin(2*pi*fc*t).*window_transmitted_i;
    transmitter_value_square_i = sum(signal_model_t_i.^2)/(length(nonzeros(window_transmitted_i))*transmitter_power_per_pulse);
    transmitted_signal_i = signal_model_t_i/sqrt(transmitter_value_square_i); 
    tpower = sum(transmitted_signal_i.^2)/length(nonzeros(window_transmitted_i)); % Transmitted signal power
    transmitter_data(:,i+1) = transmitted_signal_i.';
end

%% 
% Delay, loss and noise are added and the received signal is constructed.

receiver_power_c = zeros(1,pulse_number);
receiver_power_t = zeros(1,pulse_number);
snr_values_1 = zeros(1,pulse_number);
snr_values_2 = zeros(1,pulse_number);
noise1_powers = zeros(1,pulse_number);
noise2_powers = zeros(1,pulse_number);
radar_data_w_noise1 = zeros(length(t),pulse_number);
radar_data_w_noise2 = zeros(length(t),pulse_number);
for j = 0:pulse_number-1
    free_spl_j = (4*pi*(Range-j*rad_velocity_in*Time)*fc/light_speed)^2; % Free space path loss
    loss_j = antenna_gain_p^2*rcs/(free_spl_j*4*pi*(Range-j*rad_velocity_in*Time)^2); % Total loss from radar equation
    rpower_j = tpower*loss_j; % Received signal power
    delay_j = -2*(Range-rad_velocity_in*j*Time)/light_speed; % Delay of received signal
    delay_j_signal = -4*pi*fc*(Range-j*rad_velocity_in*Time)/light_speed;
    window_received_j = 1*((t+delay_j)>=0 & (t+delay_j)<=signal_length); % Received signal width
    signal_model_r_j = sin(2*pi*fc*t + delay_j_signal).*window_received_j; % Sinusoidal pulse signal for receiver
    receiver_value_square_j = sum(signal_model_r_j.^2)/(length(nonzeros(window_received_j))*rpower_j); 
    received_signal_j = signal_model_r_j/sqrt(receiver_value_square_j); % Received signal
    receiver_power_j = sum(received_signal_j.^2)/length(nonzeros(window_received_j)); % Received signal power
    received_signal_j_w_noise1 = awgn(received_signal_j,snr,10*log10(receiver_power_j));
    radar_data_w_noise1(:,j+1) = received_signal_j_w_noise1.';
    power_of_noise1_j = sum((received_signal_j_w_noise1-received_signal_j).^2)/length(t); % Power of noise 1
    received_signal_j_w_noise2 = received_signal_j + sqrt(power_of_noise1_j)*randn(1,length(t)); % Received signal with noise 2
    radar_data_w_noise2(:,j+1) = received_signal_j_w_noise2.';
    power_of_noise2_j = sum((received_signal_j_w_noise2-received_signal_j).^2)/length(t); % Power of noise 1
    snr1_j = 10*log10(sum(receiver_power_j)/power_of_noise1_j); 
    snr2_j = 10*log10(sum(receiver_power_j)/power_of_noise2_j); 
    receiver_power_c(j+1) = receiver_power_j;
    snr_values_1(j+1) = snr1_j;
    snr_values_2(j+1) = snr2_j;
    noise1_powers(j+1) = power_of_noise1_j;
    noise2_powers(j+1) = power_of_noise2_j;
    receiver_power_t(j+1) = rpower_j;
end

first_signal_transmitted = transmitter_data(:,1);
matched_filter = first_signal_transmitted(end:-1:1).';

radar_data_matched_filtered_w_noise1 = zeros(2*length(t)-1,pulse_number);
for q = 1:pulse_number
    matched_filtered = conv(matched_filter,radar_data_w_noise1(:,q));
    radar_data_matched_filtered_w_noise1(:,q) = matched_filtered;
end

radar_data_matched_filtered_w_noise2 = zeros(2*length(t)-1,pulse_number);
for q = 1:pulse_number
    matched_filtered = conv(matched_filter,radar_data_w_noise2(:,q));
    radar_data_matched_filtered_w_noise2(:,q) = matched_filtered;
end

first_matched_filtered_signal_with_noise1 = radar_data_matched_filtered_w_noise1(:,1);
first_matched_filtered_signal_with_noise2 = radar_data_matched_filtered_w_noise2(:,1);

figure;
tiledlayout(2,1)
ax1 = nexttile;
plot(t,radar_data_w_noise1(:,1),"Color","b")
title("Received Sinusoidal Pulse Signal before Matched Filtering Noise 1");
xlabel("Time(s)");
ylabel("Voltage(V)");

ax2 = nexttile;
plot(t,first_matched_filtered_signal_with_noise1(round(length(matched_filtered)/2):end),"Color","b");
title("Received Sinusoidal Pulse Signal after Matched Filtering Noise 1");
xlabel("Time(s)");
ylabel("Voltage(V)");
%%
% The signals are plotted.

figure;
tiledlayout(2,1);
ax3 = nexttile;
plot(t,20*log10(abs(transmitter_data(:,1))),"Color","b");
hold on;
plot(t,20*log10(abs(radar_data_w_noise1(:,1))),"Color","g");
title("Transmitted and Received Sinusoidal Pulse Signals Noise 1 without Matched Filtering");
legend("Transmitted Signal","Received Signal");
xlabel("Time(s)");
ylabel("dB Volt");

ax4 = nexttile;
plot(t,20*log10(abs(transmitter_data(:,1))),"Color","b");
hold on;
plot(t,20*log10(abs(first_matched_filtered_signal_with_noise1(round(length(matched_filtered)/2):end))),"Color","g");
title("Transmitted and Received Sinusoidal Pulse Signals Noise 1 with Matched Filtering");
legend("Transmitted Signal","Received Signal");
xlabel("Time(s)");
ylabel("dB Volt");

figure;
tiledlayout(2,1)
ax5 = nexttile;
plot(t,radar_data_w_noise2(:,1),"Color","b")
title("Received Sinusoidal Pulse Signal before Matched Filtering Noise 2");
xlabel("Time(s)");
ylabel("Voltage(V)");

ax6 = nexttile;
plot(t,first_matched_filtered_signal_with_noise2(round(length(matched_filtered)/2):end),"Color","b");
title("Received Sinusoidal Pulse Signal after Matched Filtering Noise 2");
xlabel("Time(s)");
ylabel("Voltage(V)");
%%
% The signals are plotted.

figure;
tiledlayout(2,1);
ax7 = nexttile;
plot(t,20*log10(abs(transmitter_data(:,1))),"Color","b");
hold on;
plot(t,20*log10(abs(radar_data_w_noise2(:,1))),"Color","g");
title("Transmitted and Received Sinusoidal Pulse Signals Noise 2 without Matched Filtering");
legend("Transmitted Signal","Received Signal");
xlabel("Time(s)");
ylabel("dB Volt");

ax8 = nexttile;
plot(t,20*log10(abs(transmitter_data(:,1))),"Color","b");
hold on;
plot(t,20*log10(abs(first_matched_filtered_signal_with_noise2(round(length(matched_filtered)/2):end))),"Color","g");
title("Transmitted and Received Sinusoidal Pulse Signals Noise 2 with Matched Filtering");
legend("Transmitted Signal","Received Signal");
xlabel("Time(s)");
ylabel("dB Volt");


%%
% To find the fractional delay, cross correlation is used. Both auto and cross
% correlations are computed and the highest value of the cross correlation is
% found. Then, the index of this value is obtained, which is the delay time.
time_val = zeros(1,pulse_number);
for m = 1:pulse_number
    [cross_cor,t_cor] = xcorr(transmitter_data(:,m).', radar_data_w_noise1(:,m).'); % Cross correlation of the original and delayed signals
    [auto_cor,t_acor] = xcorr(transmitter_data(:,m).'); % Auto correlation of the original signal
    %%
    %Both auto and cross correlations are plotted.

    %Plots of the auto correlation and cross correlation functions
    figure;
    tiledlayout(2,1)
    ax1 = nexttile;
    plot(t_acor/Fs,auto_cor,"LineWidth",1,"Color","b");
    title("Auto Correlation of Transmitted Signal"); 
    legend("Auto Correlation");
    xlabel("Time(s)");
    ylabel("Voltage(V)");

    ax2 = nexttile;
    plot(t_cor/Fs,cross_cor,"LineWidth",1,"Color","g");
    title("Cross Correlation of Transmitted and Received Signals");
    legend("Cross Correlation");
    xlabel("Time(s)");
    ylabel("Voltage(V)");

    guard_cell_number = 20000;
    reference_cell_number = 20000;
    pfa = 1e-6;
    N = 2*reference_cell_number;
    alpha = N*(pfa^(-1/N)-1);
    cfar_vec = [ones(1,reference_cell_number),zeros(1,guard_cell_number),ones(1,reference_cell_number)];
    cfar_threshold = conv(cross_cor,cfar_vec,"same");
    figure;
    plot(cfar_threshold)
    hold on;
    plot(cross_cor(end:-1:1))
    %%
    % Peak value of the cross correlation is found with the index.
    [peak,in] = max(cross_cor); % Peak value and index value of the peak for cross correlation function
    x_value = t_cor/Fs;
    peak_value = peak;
    unit_delay = x_value(in);

    %%
    % As can be seen, this value is only the unit part of the delay. To compute
    % the fractional part, first, the original signal is delayed by the unit part.
    % Then, the cross correlation is computed in a single sampling period by dividing
    % this range to small parts and the delay is obtained in a similar way as before.

    % Finding the cross correlations in small interval for fractional delay
    % computation
    correlation_values = [];
    for mini_delays = -dt:dt/dt_div:dt
        window_mini_del = 1*((t+unit_delay+mini_delays)>=0 & (t+unit_delay+mini_delays)<=signal_length);
        signal_mini_del = sin(2*pi*fc*(t+unit_delay+mini_delays)).*window_mini_del;
        cor_value = sum(signal_mini_del.*radar_data_w_noise1(:,m)');
        correlation_values = [correlation_values cor_value];
    end
    mini_delays = -dt:dt/dt_div:dt;

    % Plot of the cross correlation values
    figure;
    plot(mini_delays,correlation_values,"LineWidth",1,"Color","r");
    [pk,i] = max(correlation_values); % Peak value of the cross correlation and its index
    fractional_delay = mini_delays(i); % Fractional delay
    total_delay = unit_delay + fractional_delay;
    time_val(m) = total_delay;
    range_of_target = light_speed*abs(total_delay)/2;
    waitforbuttonpress
end
