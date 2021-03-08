% load data files
addpath('/Users/bryncrandles/Documents/Math-4F90-project/eeg_application/control_files/')
% folder of data files
datafolder = '/Users/bryncrandles/Documents/Math-4F90-project/eeg_application/';

% pattern of files of interest
filepattern1 = fullfile(datafolder, '*_ch24.mat'); 
filepattern2 = fullfile(datafolder, '*_ch65.mat');

% directory of data files
datafiles1 = dir(filepattern1); % channel 24 files for 401 and 441
datafiles2 = dir(filepattern2); % channel 65 files for 401 and 441

% variables
f = 10; 
fpass = 2;
t = 1:2000;
fs = 250; % sampling frequency of EEG data

for i = 2%:length(datafiles1)
    filename1 = datafiles1(i).name;
    fullfilename1 = fullfile(datafiles1(i).folder, filename1);
    filename2 = datafiles2(i).name;
    fullfilename2 = fullfile(datafiles2(i).folder, filename2);
    load(fullfilename1);
    load(fullfilename2);
    % complex demodulation
    x1 = x.*exp(-2*pi*1i*f.*t);
    y1 = y.*exp(-2*pi*1i*f.*t);
    % apply a 2 Hz low pass filter to isolate the alpha = (8, 12) frequency
    % band
    x2 = lowpass(x1, fpass, fs);
    y2 = lowpass(y1, fpass, fs);
    xphase = mod(angle(x2), 2*pi);
    yphase = mod(angle(y2), 2*pi);
    phasediff = xphase - yphase;
    figure()
    plot(phasediff)
end