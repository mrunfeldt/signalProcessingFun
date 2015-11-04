
% Inputs: default is seconds and Hz (count/seconds)
% (1) "spikeTimes" - vector - time of spikes in seconds
% (2) "signal" - timeseries sigal in discretized time.
% (3) "sRate" - sampling rate of "signal"
% (4) "tres" - desired temporal resolution of output spectrogram
% (5: optional)[v1] "chunkSize" - (in seconds) process spectrogram in temporal
% % % chunks of this size. Default is to set max at 1e5 samples.
% (6: opt.)[v2] "olap" - amount of time (sec) to overlap in temporal binning for FFT
% (7: opt.)[v3] "Nfft" - Fourier resoution. Number of bins for FFT.
% (8: opt)[v4] "faxis" - vector of frequency bins for FFT
% (9: opt)[v5] "noPlot" - STRING - set to "noPlot" if output plot is NOT desired.
% % % % % default is to produce plot

% Outputs:
% (1) "psd" - binned spectral power density
% (2) "time" - time corresponding to "psd" xaxis
% (3) "f" - frequency bins corresp. to "psd" yaxis
% (4, opt.) "fH" = figure handle (ONLY IFF Input(8) ~= "noPlot")

function [psd, time, f, varargout] = plot_pwrSpec_psth_overlay...
    (spikes, signal, sRate, tres, varargin)
% Full Input: (spikeTimes, signal, sRate, tres, chunkSize, olap, Nfft, faxis, noPlot)

dur = length(signal) ; % count of signal samples
binSize = round(tres * sRate) ; % # of samples (@sRate res) for each time bin in spectrogram

if nargin  == 4 % No optional inputs provided
% % Generate INPUT (4): "chunkSize" and "chunks" % % %
% if signal is less than "cMax," process as one chunk. Else chunk @("cMax")
    cMax = 1e5; % maximum number of samples per chunk
    if cMax < dur; chunk_size = dur ; chunks = [1 dur-1];  
    else chunk_size = cMax ;
        chunks =  [ 1:chunk_size: dur-1 ] ; % 
        if chunks(end) ~= dur-1; chunks=[chunks dur-1];end % augment if end is rounded off
    end 
    olap = binSize*2 ; % # of time samples of overlap
else % "chunkSize" was provided
    chunk_size = (varargin{1} *sRate)+1 ; % seconds, converted to samples
    % chunk_size = (chunkSize *sRate)+1 ; % for DEVEL
    chunks =  [ 1:chunk_size: dur-1 ] ; % 
    if chunks(end) ~= dur-1; chunks=[chunks dur-1];end % augment if end is rounded off
end

if  nargin < 6 % "chunkSize" was provided. "olap" was not
    olap = binSize*2 ; % # of time samples of overlap
else olap = ceil(varargin{2}*sRate) ;
end

if nargin  < 7 ; Nfft = 512; else Nfft = varargin{3};
end

% % VIP for controlling output temporal resoltion % %
window =  round( (dur - olap) / (dur/binSize) + olap ) ; 

if nargin < 8 
% NO "faxis" provided - let "spectrogram" provide frequency axis
    psd = []; time = []; tCount = 0;
    for k = 1:length(chunks)-1 % % % Calculate Spectrogram in Chunks % % % 
        tSeries = signal(chunks(k):chunks(k+1)+1) ; 
        [~,f,t,pwrD] = spectrogram(tSeries,window,olap); 
        t = t - t(1) ; % start with zero
        psd = [psd pwrD] ; time = [time t+tCount] ; % grow spectrogram PSD and Time
        tCount = tCount + t(end) + mode(diff(t)); % internal count for chunks
        disp(strcat(num2str(k/length(chunks)*100) , ' % processed')) % disp progress
        fLabel = 'Normalized frequency' ;
    end
    time = [1:length(time)].*tres; % convert to seconds
else
% Use "faxis" for frequency range % % 
faxis = varargin{4} ; 
    psd = []; time = []; tCount = 0;
    for k = 1:length(chunks)-1% % % Calculate Spectrogram in Chunks % % % 
        tSeries = signal(chunks(k):chunks(k+1)+1) ; 
        [~,f,t,pwrD] = spectrogram(tSeries,window,olap,faxis,sRate); 
        t = t - t(1) ; % start with zero
        psd = [psd pwrD] ; time = [time t+tCount] ; % grow spectrogram PSD and Time
        tCount = tCount + t(end) + mode(diff(t)); % internal count for chunks
        disp(strcat(num2str(k/length(chunks)*100) , ' % processed')) % disp progress
        fLabel = 'Frequency (kHz)' ;
    end % END (k) produce spectrogram with provided "faxis"
end

% % Determine if "noPlot" was requested % % 
if nargin == 9; pltComm = varargin{5}; else pltComm = 'yesPlot' ; end
if strcmp(pltComm,'yesPlot') % PLOT !
    
    spec=10*log10(psd);  yVals = log10(f); % PLOT on logscale (better color contrast)    
    if f(1) == 0
        faxis = [f(2:end)' f(end)+f(end-1)];
        yVals = log10(faxis);
        %yVals(1) = [] ;  spec(1,:)=[]; % remove first bin if == 0
    end
    
	binSpike = histc(spikes,time) ;  
    psth = scaleValues(binSpike,max(yVals),min(yVals)) ;
    
    fH = figure;
    surf(time,yVals,spec,'edgecolor','none');view(0,90); colormap(jet); colorbar
    hold on; bar(time,psth,'facecolor','k','edgecolor','k') % plot psth
    xlabel('Time (sec)');ylabel(fLabel)
    ylim([min(yVals) max(yVals)]);xlim([time(1) time(end)])
    varargout{1} = fH ;
end

end % END (main function

% % % Subroutine % % %
function [new] = scaleValues(x,top,bot)
    if isinf(top) || isnan(top) || isinf(bot) || isnan(bot)
        disp('Check top/bottom for inf or nan')
    end
    if bot < 0  
        bot = bot + abs(bot); top = top + abs(bot) ;
        new = bot +(x - min(x)) .* (top-bot) ./ (max(x)-min(x) );
        new = new + abs(bot) ;
    end
end