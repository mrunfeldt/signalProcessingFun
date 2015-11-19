
% Inputs: default is seconds and Hz (count/seconds)
% (1) "signal" - timeseries sigal in discretized time.
% (2) "sRate" - sampling rate of "signal"
% (3) "tres" - desired temporal resolution of output spectrogram
% (4: optional)[v1] "chunkSize" - (in seconds) process spectrogram in temporal
% % % chunks of this size. Default is to set max at 1e5 samples.
% (5: opt.)[v2] "olap" - amount of time (sec) to overlap in temporal binning for FFT
% (6: opt.)[v3] "Nfft" - Fourier resoution. Number of bins for FFT.
% (7: opt)[v4] "faxis" - vector of frequency bins for FFT
% (8: opt)[v5] "noPlot" - STRING - set to "noPlot" if output plot is NOT desired.
% % % % % default is to produce plot

% Outputs:
% (1) "psd" - binned spectral power density
% (2) "time" - time corresponding to "psd" xaxis
% (3) "f" - frequency bins corresp. to "psd" yaxis
% (4, opt.) "fH" = figure handle (ONLY IFF Input(8) ~= "noPlot")

function [psd, time, f, varargout] = spectro_chunkAndSetTime(signal, sRate, tres, varargin)
% Full Input: (signal, sRate, tres, chunkSize, olap, Nfft, faxis, noPlot)
%             (signal, sRate, tres, chunkSize, olap, Nfft, faxis, plotMe);

dur = length(signal) ; % count of signal samples
binSize = round(tres * sRate) ; % # of samples (@sRate res) for each time bin in spectrogram

if nargin  == 3 % No optional inputs provided
% % Generate INPUT (4): "chunkSize" and "chunks" % % %
% if signal is less than "cMax," process as one chunk. Else chunk @("cMax")
    cMax = 1e5; % maximum number of samples per chunk
    if cMax < dur; chunk_size = dur ; chunks = 1 ;  
    else chunk_size = cMax ;
        chunks =  [ 1:chunk_size: dur-1 ] ; % 
        if chunks(end) ~= dur-1; chunks=[chunks dur-1];end % augment if end is rounded off
    end 
    olap = binSize/2 ; % # of time samples of overlap
else % "chunkSize" was provided
    chunk_size = (varargin{1} *sRate)+1 ; % seconds, converted to samples
    % chunk_size = (chunkSize *sRate)+1 ; % for DEVEL
    chunks =  [ 1:chunk_size: dur-1 ] ; % 
    if chunks(end) ~= dur-1; chunks=[chunks dur-1];end % augment if end is rounded off
end

if  nargin < 5 % "chunkSize" was provided. "olap" was not
    olap = binSize/2 ; % # of time samples of overlap
else olap = ceil(varargin{2}*sRate) ;
end

if nargin  < 6 ; Nfft = 512; else Nfft = varargin{3};
end

% % VIP for controlling output temporal resoltion % %
window =  round( (dur - olap) / (dur/binSize) + olap ) ; 

if window > chunk_size; 
    disp(['Increase chunkSize and/or decrease olap: window = ',num2str(window), ' chunk_size = ',num2str(chunk_size)]);end

if nargin < 7 
% NO "faxis" provided - let "spectrogram" provide frequency axis
disp('yo')
    psd = []; time = []; tCount = 0;
    for k = 1:length(chunks)-1 % % % Calculate Spectrogram in Chunks % % % 
        tSeries = signal(chunks(k):chunks(k+1)+1) ; 
        [~,f,t,pwrD] = spectrogram(tSeries,window,olap);
        psd = [psd pwrD] ; time = [time t+tCount] ; % grow spectrogram PSD and Time
        tCount = tCount + t(end-1); % internal count for chunks
        disp(strcat(num2str(k/length(chunks)*100) , ' % processed')) % disp progress
        fLabel = 'Normalized frequency' ;
    end
else
% Use "faxis" for frequency range % % 
faxis = varargin{4} ; 
    psd = []; time = []; tCount = 0;
    for k = 1:length(chunks)-1% % % Calculate Spectrogram in Chunks % % % 
        tSeries = signal(chunks(k):chunks(k+1)+1) ; 
        [~,f,t,pwrD] = spectrogram(tSeries,window,olap,faxis,sRate); 
        psd = [psd pwrD] ; time = [time t+tCount] ; % grow spectrogram PSD and Time
        tCount = tCount + t(end-1); % internal count for chunks
        disp(strcat(num2str(k/length(chunks)*100) , ' % processed')) % disp progress
        fLabel = 'Frequency (kHz)' ;
    end % END (k) produce spectrogram with provided "faxis"
end

% % Determine if "noPlot" was requested % % 
if nargin == 8; pltComm = varargin{5}; else pltComm = 'yesPlot' ; end

if strcmp(pltComm,'yesPlot')
    spec=10*log10(psd);yVals = log10(f);
    
    fH = figure;
    surf(time,yVals,spec,'edgecolor','none');view(0,90); 
    colormap(jet); colorbar
    xlabel('Time (sec)');ylabel(fLabel)
    ylim([min(yVals) max(yVals)]);xlim([time(1) time(end)])
    varargout{1} = fH ; drawnow
end
