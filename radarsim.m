close all;
clear all;
clc;
fc = 1e10;                              % Operating frequency
pd = 0.9;                               % Probability of detection
pfa = 1e-6;                             % Probability of false alarm
coor_of_target = [-5000;-5000;5000];    % Target coordinates
speed_of_target = [0;0;0];              % Target speed
tgt_rcs = 100;                          % Target radar cross section
pulse_width = 200e-6;                   % Pulse width
pulse_bandwidth = 1/pulse_width;        % Pulse bandwidth
tx_power = 25e6;                        % Transmitter power
antenna_gain = 33;                      % Antenna gain
system_noise_T = 290;                   % System noise temperature
num_of_antenna_elements = 6;            % Number of antenna elements
num_of_pulses = 60;                     % Number of pulses

prop_speed = physconst('LightSpeed');   % Propagation speed
lambda = prop_speed/fc;                 % Wavelength
range_res = prop_speed*pulse_width/2;   % Range resolution            
noise_bandwidth = pulse_bandwidth;      % Noise Bandwidth

figure;
rocsnr([4 8 12 15],"MaxPfa",pfa)    
snr_min = albersheim(pd, pfa);          % Snr_min calculation
max_range = (tx_power*(db2pow(antenna_gain))^2*lambda^2*tgt_rcs/(4^3*pi^3*db2pow(snr_min) ...
    *physconst("Boltzmann")*system_noise_T*noise_bandwidth))^(1/4);  % Maximum range calculation
prf = prop_speed/(2*max_range);         % Pulse repetition frequency
fs = 2*pulse_bandwidth*prf;             % Sampling rate

antenna = phased.GaussianAntennaElement('FrequencyRange',[5e9 15e15],"Beamwidth",[10 10]);
figure;
pattern(antenna,fc,-180:180,CoordinateSystem="polar");
array = phased.ULA("Element",antenna,"ElementSpacing",lambda/2,"NumElements",num_of_antenna_elements,"ArrayAxis","x");
figure;
pattern(array,fc,-180:180,CoordinateSystem="polar");

collector = phased.Collector("Sensor",array,"PropagationSpeed",prop_speed,"OperatingFrequency",fc);
radiator = phased.Radiator("Sensor",array,"PropagationSpeed",prop_speed,"OperatingFrequency",fc);

transmitter = phased.Transmitter("PeakPower",tx_power,"Gain",antenna_gain,"InUseOutputPort",true);
receiver = phased.ReceiverPreamp("Gain",antenna_gain,"SampleRate",fs,"EnableInputPort",true);

doa = phased.MVDREstimator("SensorArray",array,"OperatingFrequency",fc,"DOAOutputPort",true);

target = phased.RadarTarget("MeanRCS",tgt_rcs,"OperatingFrequency",fc);
env = phased.FreeSpace("OperatingFrequency",fc,"SampleRate",fs,"TwoWayPropagation",true);
waveform = phased.LinearFMWaveform('PulseWidth',pulse_width,'PRF',prf,'SampleRate',fs);
pulse = waveform();
t = (0:size(pulse)-1)/fs;
figure;
plot(t,real(pulse))

transmitterplatform = phased.Platform("InitialPosition",[0;0;0]);
receiverplatform1 = phased.Platform("InitialPosition",[0;0;0]);
receiverplatform2 = phased.Platform("InitialPosition",[10000;0;0]);
targetplatform = phased.Platform("InitialPosition",coor_of_target,"Velocity",speed_of_target);

sv = phased.ScenarioViewer("Name","Radar 1","ShowBeam","All","BeamWidth",antenna.Beamwidth,"BeamRange",max_range,"UpdateRate",0.1, ...
    "ShowPosition",true,"ShowPosition",true,"ShowSpeed",true,"ShowAltitude",true,"ShowLegend",true);

fast_time_grid = unigrid(0,1/fs,1/prf,'[)');
slow_time_grid = (0:num_of_pulses-1)/prf;

rxpulses = zeros(numel(fast_time_grid),num_of_pulses,num_of_antenna_elements);

for m = 1:num_of_pulses

    [transmitterpos,transmittervel] = step(transmitterplatform,1);
    [receiverpos,receivervel] = step(receiverplatform2,1);
    [targetpos,targetvel] = step(targetplatform,1);

    [targetranget,targetanglet] = rangeangle(targetpos,transmitterpos);
    [targetranger,targetangler] = rangeangle(targetpos,receiverpos);

    sv.BeamSteering = [targetanglet targetangler];
    step(sv,[transmitterpos,receiverpos],[transmittervel,receivervel],targetpos,targetvel);
    
    pulse = waveform();
    [tsignal,status] = transmitter(pulse);
    tprop = radiator(tsignal,targetanglet);
    signal_at_target = env(tprop,transmitterpos,targetpos,transmittervel,targetvel);

    signal_reflected = target(signal_at_target);

    signal_at_receiver = env(signal_reflected,targetpos,receiverpos,targetvel,receivervel);
    received_signal = collector(signal_at_receiver,targetangler);
    final_signal = receiver(received_signal,~(status>0));
    rxpulses(:,m,:) = final_signal;
    doas = doa(final_signal);
    doas = broadside2az(doas,10);
    plotSpectrum(doa)
    waitforbuttonpress
    pause(0.1)
end
npower = noisepow(noise_bandwidth,receiver.NoiseFigure,receiver.ReferenceTemperature);
threshold = npower * db2pow(npwgnthresh(pfa,num_of_pulses,'noncoherent'));

num_pulse_plot = 2;
figure;
plot(real(rxpulses(:,1,1)));
figure;
helperRadarPulsePlot(rxpulses,threshold,fast_time_grid,slow_time_grid,num_pulse_plot);
matchingcoeff = getMatchedFilter(waveform);
matchedfilter = phased.MatchedFilter('Coefficients',matchingcoeff,'GainOutputPort',true);
[rxpulses, mfgain] = matchedfilter(rxpulses);

matchingdelay = size(matchingcoeff,1)-1;
rxpulses = buffer(rxpulses(matchingdelay+1:end),size(rxpulses,1));

threshold = threshold * db2pow(mfgain);

figure;
helperRadarPulsePlot(rxpulses,threshold,fast_time_grid,slow_time_grid,num_pulse_plot);