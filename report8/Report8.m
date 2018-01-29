close all; clear all; clc
 
load('spike_neural.mat') % Load the neural_sig signal
L = length(neural_sig); % Duration of the signal in samples
fs = 10240; % Sampling frequency in Hz
WinSize = [0.2:0.125:2]; % Window size in seconds 
WinSize = round(WinSize.*fs); % Window size in samples
 
f_ax = (-pi:2*pi/fs:pi-2*pi/fs)./(2*pi).*fs; % Frequency axis in Hz
for uu = 1 : length(WinSize)
    window = rectwin(WinSize(uu))'; % Window type
    for n = 1:35 %Why is n=35?
        wind_signal=neural_sig((n-1)*WinSize(uu)+(1:WinSize(uu))).*window;
        Segm_spect{uu}(n,:) = fftshift(abs(fft(wind_signal,fs)).^2)./WinSize(uu); % Periodogram
    end    
    variabil(uu) = max(var(Segm_spect{uu})); % Variance of estimate (max variance over the frequency axis)
end;
%% 

%window size 0.2
Mean_periodogram1 = mean(Segm_spect{1}); %Found by doing 0.2*fs and looking at Win_Size to see the 

%number of samples this corresponds to
figure; plot(f_ax, Mean_periodogram1);
title('Mean of the periodogram for rectangular window size = 0.2s');
xlabel('Frequency (Hz)'); ylabel('Amplitude (AU)');

peak = find(f_ax == 0)
Bias02 = sum(Mean_periodogram1(peak-5:peak+5))/Mean_periodogram1(peak)

%% 

%window size 0.5
Mean_periodogram2 = mean(Segm_spect{4}); 
figure; plot(f_ax, Mean_periodogram2);
title('Mean of the periodogram for rectangular window size = 0.5s');
xlabel('Frequency (Hz)'); ylabel('Amplitude (AU)');
Bias05 = sum(Mean_periodogram2(peak-5:peak+5))/Mean_periodogram2(peak)

%% 

%window size 2
Mean_periodogram3 = mean(Segm_spect{15}); 
figure; plot(f_ax, Mean_periodogram3);
title('Mean of the periodogram for rectangular window size = 2s');
xlabel('Frequency (Hz)'); ylabel('Amplitude (AU)');
Bias2 = sum(Mean_periodogram3(peak-5:peak+5))/Mean_periodogram3(peak)


%% 

%variance
plot(WinSize/10240, variabil);
title('Variance of estimation of the periodogram for rect window');
xlabel('Window size (s)'); ylabel('Variance');
