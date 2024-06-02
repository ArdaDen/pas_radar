%% IQ DATA
% I data and Q data values are found.
for l = 1:pulse_number
    q_data = 2*radar_data_w_noise1(:,l).*cos(2*pi*fc*t); % Q Data
    i_data = 2*radar_data_w_noise1(:,l).*sin(2*pi*fc*t); % I Data
    fft_q_data = fft(q_data); %FFT of Q Data
    fft_i_data = fft(i_data); %FFT of I Data
    freq_scale = linspace(0,Fs,length(t)); % Frequency scale of FFT
    %%
    % Their plots

    % Plot of Q Data
    figure;
    subplot(3,2,1)
    plot(t,q_data,"LineWidth",1,"Color","b");
    title("Q Data");
    xlabel("Time(s)")
    ylabel("Voltage(V)")

    % Plot of I Data
    subplot(3,2,2)
    plot(t,i_data,"LineWidth",1,"Color","b");
    title("I Data");
    xlabel("Time(s)")
    ylabel("Voltage(V)")
    % Plot of absolute
    subplot(3,2,3)
    plot(freq_scale,abs(fft_q_data)/length(t),"LineWidth",1,"Color","b");
    title("FFT of Q Data")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")

    % Plot of absolute
    subplot(3,2,4)
    plot(freq_scale,abs(fft_i_data)/length(t),"LineWidth",1,"Color","b");
    title("FFT of I Data")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")

    % Plot of shifted FFT of Q Data
    shifted_freq_scale = linspace(-Fs/2,Fs/2,length(t)); % Shifted frequency scale
    subplot(3,2,5)
    plot(shifted_freq_scale,fftshift(abs(fft_q_data)/length(t)),"LineWidth",1,"Color","b");
    title("Shifted FFT of Q Data")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")

    % Plot of shifted FFT of I Data
    subplot(3,2,6)
    plot(shifted_freq_scale,fftshift(abs(fft_i_data)/length(t)),"LineWidth",1,"Color","b");
    title("Shifted FFT of I Data")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")
    %%
    % The signals are passed from a lowpass filter to obtain the angle. Finding
    % the valus of the Q and I Data and then calculating the angle. From there, the
    % delay of the signal can be calculated.

    low_pass_filter = designfilt("lowpassfir","FilterOrder",128,"DesignMethod","maxflat","HalfPowerFrequency",fc*0.02,"SampleRate",Fs);
    % low_pass_filter = fir1(32,200/Fs);
    % fvtool(low_pass_filter,"MagnitudeDisplay","magnitude")
    q_data_filtered = filter(low_pass_filter,q_data); % Sine theta
    i_data_filtered = filter(low_pass_filter,i_data); % Cosine theta
    % q_data_filtered = lowpass(q_data,0.0001*fc,Fs); % Sine theta
    % i_data_filtered = lowpass(i_data,0.0001*fc,Fs); % Cosine theta

    % Plot of Q Data after filtering
    figure;
    subplot(3,2,1)
    plot(t,q_data_filtered,"LineWidth",1,"Color","b");
    title("Q Data after Lowpass Filtering")
    xlabel("Time(s)")
    ylabel("Voltage(V)")

    % Plot of I Data after filtering
    subplot(3,2,2)
    plot(t,i_data_filtered,"LineWidth",1,"Color","b");
    title("I Data after Lowpass Filtering")
    xlabel("Time(s)")
    ylabel("Voltage(V)")

    % FFT of Q Data after filtering
    fft_q_data_filtered = fft(q_data_filtered); % FFT

    % FFT of I Data after filtering
    fft_i_data_filtered = fft(i_data_filtered); % FFT

    % Plot of absolute after filtering Q Data
    subplot(3,2,3)
    plot(freq_scale,abs(fft_q_data_filtered)/length(t),"LineWidth",1,"Color","b");
    title("FFT of Q Data after Filtering")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")

    % Plot of absolute after filtering I Data
    subplot(3,2,4)
    plot(freq_scale,abs(fft_i_data_filtered)/length(t),"LineWidth",1,"Color","b");
    title("FFT of I Data after Filtering")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")

    % Plot of shifted FFT after filtering Q Data
    subplot(3,2,5)
    plot(shifted_freq_scale,fftshift(abs(fft_q_data_filtered)/length(t)),"LineWidth",1,"Color","b");
    title("Shifted FFT of Q Data after Filtering")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")

    % Plot of shifted FFT after filtering I Data
    subplot(3,2,6)
    plot(shifted_freq_scale,fftshift(abs(fft_i_data_filtered)/length(t)),"LineWidth",1,"Color","b");
    title("Shifted FFT of I Data after Filtering")
    xlabel("Frequency(Hz)")
    ylabel("Amplitude")
    waitforbuttonpress
end