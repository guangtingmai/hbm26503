function [WaveOut b a] = butter_filtfilt(WaveIn, Fs, ftype, cut_off, order)
% Zero-phase butterworth filter
% Guangting Mai
% NIHR Nottingham BRC
%
% -- WaveIn should be [(channel) x (sample points)]
% -- Fs - sample frequency (Hz)
% -- ftype should be 'high', 'low' or 'band'
% -- cut_off should be a single frequency in Hz for 'high' and 'low'; a 
%    vector [low cutoff, high cutoff] for 'band'
% -- order - filter order

if ~(strcmp(ftype,'high')|strcmp(ftype,'low')|strcmp(ftype,'band'))
    error('filter type should be high, low or band');
end

if strcmp(ftype,'high')|strcmp(ftype,'low')&(length(cut_off)~=1)
    error('cutoff frequency should be a single frequency for lowpass and highpass filter');
end

if strcmp(ftype,'band')&(length(cut_off)~=2)
    error('cutoff frequencies should include two frequencies for bandpass filter');
end

% Set Nyquist frequency to sample rate/2
Nyquist = Fs/2;
% normalise cut off frequency in range 0 to 1 relative to Nyquist frequency
normalFc = cut_off/Nyquist;
% design butterworth filter
if strcmp(ftype,'high')|strcmp(ftype,'low')
    [b,a]=butter(order,normalFc,ftype); % lowpass/highpass
else
    [b,a]=butter(order,normalFc); % bandpass
end

% use filtfilt to lose phase distortion
% and apply filter coefficients to input signal 
[rows,cols] = size(WaveIn);
WaveIn = double(WaveIn);
for i = 1:rows
    WaveOut(i,:) = filtfilt(b,a,WaveIn(i,:));
end