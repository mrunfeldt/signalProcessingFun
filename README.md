# signalProcessingFun
MATLAB functions for signal processing. Functions were originally written for analysis of neural activity in primary auditory cortex, but are also generalizable. Written by Melissa Jane Runfeldt, 2015_10_05

"spectro_chunkAndSetTime.m" - MATLAB function to generate a spectrogram from a timeseries input. Allows user to control temporal and spatial (frequency bins) resolution of spectrogram.
** "runScript_callSpecFun.m" - script to run example for use of above function. 
** "callwaveform.m" - sound pressure wave - example data for runScript 

"plot_pwrSpec_psth_overlay.m" - MATLAB function to generate spectrogram with enhanced user control. Also plots PSTH of neuronal spikes on top of spectrogram. Spikes are binned at same temporal resolution as spectrogram output.
** "runScript_plotSpecPSTH_overlay.m" - script to run example for use of above function. 
** "spikes.m" - neuron spike times - example data for runScript
