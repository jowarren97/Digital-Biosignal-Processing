% Third tutorial, Friday 27.10.2017
clear all
close all
clc
 
% Signal loading
load('EEG.mat');
 
% Sampling frequency
fsamp = 512;
 
% Select duration to analyse
Duration = 1; % Duration in seconds (max 15 seconds)
Duration = round(Duration*fsamp);
EEG = EEG(1:Duration); 
 
% Signal duration in samples
L = length(EEG);
 
% Plot the EEG signal
time_ax = [1:L]./fsamp;
figure(1)
plot(time_ax, EEG);
xlabel('Time (s)')
title(['EEG signal'])
ylabel('Amplitude EEG (Arbitrary Units)')
 
% Compute DFT of the channel with DC offset removal
X1 = fft( EEG - mean(EEG) );
 
% Compute PSD (power spectral density) of the channel
PSD1 = fftshift(abs(X1).^2);
 
% Build the frequency axis in radians
freq_a_rad = [-pi+pi/L:2*pi/L:pi-pi/L];
% Convert the frequency axis in Hz
freq_a_Hz = freq_a_rad./(2*pi).*fsamp;
 
% Plot the PSD of the channel with frequencies in radians and Hz
figure(2), subplot(2,1,1), plot(freq_a_rad,PSD1);
xlabel('Frequency (radiants)')
title(['Power Spectral Density of EEG; Duration = 1s'])
ylabel('PSD (Arbitrary Units)')
xlim([-1 1])
 
figure(2), subplot(2,1,2), plot(freq_a_Hz,PSD1);
xlabel('Frequency (Hz)')
title(['Power Spectral Density of EEG; Duration = 1s'])
ylabel('PSD (Arbitrary Units)')
xlim([(-512/(2*pi)) (512/(2*pi))])

 
%% Compute the percentage of power in different subbands
%[Here please complete with instructions for computing the relative power in the frequency bands.]
N = length(PSD1);

totsum = 0;
partialsum = 0;
k1 = find(freq_a_Hz>0.5);
k1 = k1(1)
k2 = find(freq_a_Hz<4);
k2 = k2(end)

for k = 1:N/2
    totsum = totsum + PSD1(k);
end

for k = k1:k2
    partialsum = partialsum + PSD1(k);
end

delta = partialsum*100/totsum

%% Compute the percentage of power in different subbands
%[Here please complete with instructions for computing the relative power in the frequency bands.]
N = length(PSD1);

totsum = 0;
partialsum = 0;
k1 = find(freq_a_Hz>4);
k1 = k1(1)
k2 = find(freq_a_Hz<8);
k2 = k2(end)

for k = 1:N/2
    totsum = totsum + PSD1(k);
end

for k = k1:k2
    partialsum = partialsum + PSD1(k);
end

theta = partialsum*100/totsum

%% Compute the percentage of power in different subbands
%[Here please complete with instructions for computing the relative power in the frequency bands.]
N = length(PSD1);

totsum = 0;
partialsum = 0;
k1 = find(freq_a_Hz>8);
k1 = k1(1)
k2 = find(freq_a_Hz<13);
k2 = k2(end)

for k = 1:N/2
    totsum = totsum + PSD1(k);
end

for k = k1:k2
    partialsum = partialsum + PSD1(k);
end

alpha = partialsum*100/totsum

%% Compute the percentage of power in different subbands
%[Here please complete with instructions for computing the relative power in the frequency bands.]
N = length(PSD1);

totsum = 0;
partialsum = 0;
k1 = find(freq_a_Hz>13);
k1 = k1(1)
k2 = find(freq_a_Hz<30);
k2 = k2(end)

for k = 1:N/2
    totsum = totsum + PSD1(k);
end

for k = k1:k2
    partialsum = partialsum + PSD1(k);
end

beta = partialsum*100/totsum

%% Compute the percentage of power in different subbands
%[Here please complete with instructions for computing the relative power in the frequency bands.]
N = length(PSD1);

totsum = 0;
partialsum = 0;
k1 = find(freq_a_Hz>30);
k1 = k1(1)
k2 = find(freq_a_Hz<42);
k2 = k2(end)

for k = 1:N/2
    totsum = totsum + PSD1(k);
end

for k = k1:k2
    partialsum = partialsum + PSD1(k);
end

gamma = partialsum*100/totsum




