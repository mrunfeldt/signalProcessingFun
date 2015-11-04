

% Example of plotting a spectrogram from a long timeseries
% Instructions: 
% (1) download "callwaveform.mat" and "spikeTimes" and set "sigPath" to
% % % location of files
% (2) add "spectro_chunkAndSetTime" function to MATLAB path

clear all

sigPath = 'C:\Users\mel\Desktop\gitProj\signalProcessingFun\'; % Path to signal
call=load([sigPath,'callwaveform.mat']);  
load([sigPath,'spikes.mat']); % load spike times (

signal = call.callwaveform; % timeseries INPUT
sRate = 48000; % Temporal Sampling Rate for Call : 48 kHz
tres = 1e-3 ; % Desired temporal resolution (sec)

% % % Option (1): Set Input Parameters % % %
Nfft = 512; % Fourier resolution
faxis = logspace(log10(50),log10(4e4)) ;
chunkSize = 0.75; % in seconds
olap = 0.1; % in seconds

[psd, time, f, fH] = plot_pwrSpec_psth_overlay...
    (spikes, signal, sRate, tres, chunkSize,olap,Nfft,faxis) ;


%% % % Option (2): Allow default parameters % % 
[psd, time, f, fH] = plot_pwrSpec_psth_overlay...
    (spikes, signal, sRate, tres) ;