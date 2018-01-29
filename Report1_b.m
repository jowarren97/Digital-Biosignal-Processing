close all
clear all
fclose('all')

M = 1;

while M<=8
    load('Signals_EMG.mat'); % Loading the recorded EMGs (two channels)
    n_step = 100; % Number of steps in the loop for optimal alignment
    stepCount = 0.2; % Size of the step of non-integer delay applied to one of the signals for alignment
    ied=24; % Interelectrode distance of the recordings in mm
    Fs = 2048; % Recording sampling frequency
    Ts = 1/Fs; % Sampling interval

    % Downsampling by an integer
    
    channel1 = channel1(1:M:end); % Downsampling first channel
    channel2 = channel2(1:M:end); % Downsampling second channel
    Fs=Fs/M;
    Ts=Ts*M;
    % End downsampling

    timeAxis=[1:length(channel1)].*Ts.*1000; % Definition of time axis in ms
    freqAxis=fftshift([-0.5:1/(length(channel1)):0.5-1/(length(channel1))]); % Definition of discrete frequency axis

    channel1_ft = fft(channel1); % Fourier transform of the first channel

    for uu = 1 : n_step
        channel1_dt = (channel1_ft).*exp(-i*2*pi*stepCount*uu*freqAxis); % complex exponential multiplication (delay in frequency domain)
        channel1_dt = real(ifft((channel1_dt))); % inverse transform to go back to the time domain 
        MSE_vect(uu)= sum((channel1_dt - channel2).^ 2)./sum(channel2.^ 2).*100; % normalized mean square error between aligned signals
        delay(uu) = stepCount*uu; % Imposed delay in samples
    end;
    xlabel('Time (ms)')
    ylabel('Signal amplitude (AU)')

    % Identification of the optimal delay (minimum mean sqaure error)
    [MSEopt, optDelay] = min(MSE_vect);

    figure(1)
    plot(delay*Ts*1000, MSE_vect); 

    hold on;
    
    M = M*2; % Downsampling factor
end;

legend(['M=1'],['M=2'],['M=4'],['M=8']);
xlabel('Delay/ ms')
ylabel('Mean Squared Error')
