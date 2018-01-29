% Seventh tutorial, Friday 01.12.2017.
close all; clear all; clc
 
load('EEG.mat') % Load the EEG signal
fs = 200; % Sample frequency in Hz
L = length(EEG); % Duration of the signal in samples
 
winStep = 0.025; % window size increase in ms
%winSizeVector = floor([1:(L)/(winStep*fs)]*winStep*fs); % Vector containing window sizes ranging from 25ms to 15s in 25ms step increases 
winSizeVector = floor([1:(L)/(winStep*fs)]*winStep*fs); % Vector containing window sizes ranging from 25ms to 15s in 25ms step increases 

% values for different window lengths
for iter_winSize = 1:numel(winSizeVector) %iterate thru elements of winSizeVector
    %calculate the first window values
    V_fixWin(1) = var(EEG(1:winSizeVector(iter_winSize))); % variance; V_fixWin(i,j) contains variance for winSizeVector(i) at jth interval section
    M_fixWin(1) = mean(EEG(1:winSizeVector(iter_winSize))); % mean
 
    % calculate the remaining values for the same window size
    for iterSig=2:floor((length(EEG))/winSizeVector(iter_winSize))-1
        start_idx=(iterSig-1)*winSizeVector(iter_winSize);
        %V_fixWin(iterSig) = var(EEG(next interval));
        V_fixWin(iterSig) = var(EEG(iterSig*winSizeVector(iter_winSize):iterSig*winSizeVector(iter_winSize)+winSizeVector(iter_winSize)));% Calculate variance for the given window length
        M_fixWin(iterSig) = mean(EEG(iterSig*winSizeVector(iter_winSize):iterSig*winSizeVector(iter_winSize)+winSizeVector(iter_winSize)));% Calculate mean for the given window length
    end
    V{iter_winSize,:} = V_fixWin;% store variance for the given window length; V(i,:) stores variance for window size 25*i ms (?)
    M{iter_winSize,:} = M_fixWin;% store mean for the given window length
    
    V_avg(iter_winSize) = mean(V_fixWin);% calculate the mean value of the variance for the given window length
    M_avg(iter_winSize) = mean(M_fixWin);% calculate the mean value of the mean for the given window length
        
    V_fixWin=[]; %empty contents for next iteration
    M_fixWin=[];
end
%%
 winSizeIdx=[1.5/0.025]; %?
 x_axis=[1:length(V{winSizeIdx,:})].*(winSizeVector(winSizeIdx)/fs);
 figure, plot(x_axis,V{winSizeIdx,:},'k','LineWidth',2);
 xlabel('Time [s]'), ylabel('[AU]')
 title('Estimates of the variance for window length = 1.5s')

 figure, plot(x_axis,M{winSizeIdx,:},'k','LineWidth',2);
 xlabel('Time [s]'), ylabel('[AU]')
 title('Estimates of the mean for window length = 1.5s')
 %% 
 
 x_ax_length = [1:(L)/(winStep*fs)]*winStep; %build the x axis for the plot
%Plot the variance with respect to the window length
 figure, plot(x_ax_length,V_avg,'k','LineWidth',2),hold on
 xlabel('Duration of the windows [s]'), ylabel('[AU]')
 title('Values of variance depending on the window size')

 
%Plot the variance with respect to the window length
 figure, plot(x_ax_length,M_avg,'k','LineWidth',2),hold on
 xlabel('Duration of the windows [s]'), ylabel('???')
 title('Values of mean depending on the window size')
 
figure
histogram(EEG,'Normalization','probability');
xlabel('bins'), ylabel('[AU]')
title('PDF of the random process')

