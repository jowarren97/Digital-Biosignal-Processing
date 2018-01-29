% Fifth Tutorial, Friday 10.11.2017.
clear all; close all; clc;
 
load 'emg'; % Load the EMG signal
fs = 1600; % Sampling frequency
L = length(emg); % Duration of the signal in samples
time_ax=[0:1/fs:(L-1)/fs]; % Time axis of the signal in seconds

%figure;
%plot(time_ax, abs(emg));

fc = 2; % Filter cut-off frequency in (Hz)
wc = 2*fc/fs*pi; % Normalized cut-off frequency
 
% Generate Sinc function with a given sampling rate
t = -floor(L):floor(L); % Form the time axis of the Sinc function (theoretically it can be infinitely long)
sinc_func=wc*sinc(wc*t);
 
figure; plot(t, sinc_func);
title('Sinc function');
xlabel('n');
ylabel('AU');
 
f_snc=fftshift(fft(sinc_func)); % Find the DFT of the sinc function 
f_ax =(-pi+pi/length(sinc_func):2*pi/length(sinc_func):pi-pi/length(sinc_func))./pi; % Frequency axis for the DFT of sinc function
 
figure; plot(f_ax,abs(f_snc));
title('Magnitude of discrete time Fourier transform of the Sinc');
xlabel('Frequency (rad)')
ylabel('AU');
 
% Create the FIR filter by truncating the very long sinc function and by using
% one of the following windows: hanning; rectwin;
FIR_duration=100;
truncation_section=floor(length(sinc_func)/2-FIR_duration/2):floor(length(sinc_func)/2+FIR_duration/2);
fir_filt = sinc_func(truncation_section).*hanning(length(truncation_section))';
 
figure;freqz(fir_filt);
 
% Create moving avarege filter with the lenght of MA_coef_num
MA_coef_num=100;
MA = ones(1,MA_coef_num)/MA_coef_num;

figure;freqz(MA);

rect_emg = abs(emg);
env_FIR=conv(rect_emg,fir_filt, 'same');
figure;
plot(time_ax, env_FIR);
title('Hanning filtered EMG, filter length = 100')
 
% Filter the rectified EMG using moving avarage filter
env_MA =conv(rect_emg,MA, 'same');
figure;
plot(time_ax, env_MA);
title('MA filtered EMG, filter length = 10')
