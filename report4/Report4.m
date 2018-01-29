% Forth Tutorial, Friday 03.11.2017.
clear all; close all; clc;
 
load('ECG.mat') % Load the signal
ECG=ECG-mean(ECG); % Remove mean
fs = 1000; % Sample frequency in Hz
t_ax = (1:length(ECG))/fs; % Time axis of the signal

% Plot the signal
figure(1), plot(t_ax, ECG );
title('ECG signal');
xlabel('Time (s)'),ylabel('AU');
xlim([0 5]);
 
ECG_duration=size(ECG,2); % Duration of the ECG signal in samples
f_ax=fftshift([-pi+pi/ECG_duration:2*pi/ECG_duration:pi-pi/ECG_duration]); % Frequency axis for DFT
F_ECG = fft(ECG); % Calculate DFT 
 
% Plot the DFT of the signal
figure(2), plot(f_ax,abs(F_ECG));
title('Magnitude of discrete time Fourier transform of the ECG signal');
xlabel('Frequency (rad)')
ylabel('AU');
 
% Create moving average filter
MA_coef_num = 50;
MA = ones(1,MA_coef_num)/MA_coef_num;

ECG_filt = conv(ECG,MA,'same');
 
% Plot the filtered ECG signal 
figure(3); plot(t_ax, ECG_filt, 'r');
title('Filtered ECG signal, filter length = 50');
xlabel('Time (s)'),ylabel('AU');
xlim([0 5]);
 
F_ECG_filt = fft(ECG_filt);% Calculate DFT of the filtered ECG
 
% Plot the Fourier transform of the filtered signal
figure(4), plot(f_ax,abs(F_ECG_filt));
title('Magnitude of discrete time Fourier transform of the filtered ECG signal, filter length = 50');
xlabel('Frequency (rad)')
ylabel('AU');
 
% Create system function of the z-transform of the MA filter
H = tf(MA,1,1/fs,'variable','z^-1');
 
% Evaluate the magnitude and the phase of the filter with respect to the normalized frequency 
figure(5), freqz(MA,1);
% Go from coefficients to zeroes and poles of the moving average filters
[MA_zeros,MA_poles] = tf2zp(H.Numerator{1,1},(H.Denominator{1,1})); 
% Z-plane of the moving average filter
figure(6), zplane(MA_zeros,MA_poles)
title('Filter coefficients represented in the z-plane, filter length = 50')
