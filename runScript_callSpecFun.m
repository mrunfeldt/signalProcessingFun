

% Example of plotting a spectrogram from a long timeseries
% Instructions: 
% (1) download "callwaveform.mat" and set "sigPath" to location
% (2) add "spectro_chunkAndSetTime" function to MATLAB path

clear all

sigPath = 'C:\Users\mel\Desktop\gitProj\signalProcessingFun\callwaveform.mat'; % Path to signal
call=load(sigPath);  

signal = call.callwaveform; % timeseries INPUT
% signal = signal(1:round(length(signal/4))); % option to truncate signal

sRate = 48000; % Temporal Sampling Rate for Call : 48 kHz
tres = 1e-3 ; % Desired temporal resolution (sec)
Nfft = 512; % Fourier resolution
faxis = logspace(log10(50),log10(4e4)) ;
chunkSize = 0.75; % in seconds
olap = 0.1; % in seconds

[psd, time, f, fH] = spectro_chunkAndSetTime...
    (signal, sRate, tres, chunkSize,olap,Nfft,faxis) ;