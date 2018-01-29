% Ninth tutorial, Friday 15.12.2017.
clear all; close all; clc;
 
load('spike_neural(1).mat') %load the signal
nerual_signal = neural_sig - mean(neural_sig); %eliminate the mean
fs = 10240; % sampling frequency in Hz
 
win_size_sec=0.250; %in seconds
win_size=fix(win_size_sec*fs); % window size in samples
win_type=hanning(win_size);
 
overlap_size_percent=0; % size of overlap as fraction of window size
overlap_size=fix(overlap_size_percent *win_size); % size of overlap in samples
 
%Obtain periodogram using Welch approach with a given window type, size and
%overlap
[neural_periodogram,f_ax]=pwelch(nerual_signal,win_type,overlap_size,2*fs,fs);
figure; plot(f_ax,neural_periodogram);
title('Periodogram of the neural signal, hanning win size = 250ms');
xlabel('Frequency (Hz)')
ylabel('AU');
xlim([0 50]);

[Y,I] = max(neural_periodogram);
firing_rate = f_ax(I);

%% 
spike_count = 0;

for i = 1:length(nerual_signal)
    if (nerual_signal(i) > -0)
       spike_count = spike_count + 1;
    end
end

real_firing_rate = spike_count/(length(nerual_signal)/fs)
