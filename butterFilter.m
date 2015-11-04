% "timeSeries" - digital signal to process
% "frameRate" - sampling rate of timeSeries in Hz
% "order" - order of butter filter (higher is sharper/flatter passband)
% "highLow" - string determines if filter is high or low pass   
% MJRunfeldt, updated 2015_10_02

function [signalOut] = butterFilter(timeSeries,frameRate,cutoff,order,highLow)

    nyquist = frameRate/2;
    fMax = cutoff/nyquist;
    [b,a] = butter(order,fMax,highLow); % generate filter
    signalOut=filtfilt(b,a,timeSeries); % apply filter to timeseries

end