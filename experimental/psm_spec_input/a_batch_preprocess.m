%   a_batch_preprocess - preprocess data, extract passband and save
%
%   Note that this script uses a local function, calc_spec_omega(), defined
%   at the end of the file.
%
%   2021-09-16  MHS - Initial version
%

%% Reset
close all
clearvars

%% Define path to data files
% Files are assumed to be named so that they are sorted in the same order
% that they were originally generated (i.e. the first file in a
% alphabetically sorted list contains the first part of the scan, etc.).

% Update this path if re-processing on your local computer
dataPath = 'E:\Synaptus\RadarAntarcticaData\';  

%% Create folder for saving output (if not already in place)
savePath = [dataPath 'Preprocessed\'];
if not(exist(savePath,'dir'))
    mkdir(savePath)
end

%% Constants
fs = 500e6;             % Sampling frequency
fLow = 170e6;           % Lower edge, passband
fHigh = 230e6;          % Upper edge, passband
xStep_orig = 0.1920;    % Original spatial sampling interval
downsamp_factor = 2;    % Optional downsampling factor. Set to 1 to disable downsampling. 

%% Get all file names
matFiles = dir([dataPath '*.mat']);
matFiles = {matFiles.name};
nFiles = length(matFiles);

%% Load data and calculate spectra
% for ii = 1:nFiles
for ii = 2
    % Status update
    disp(['Loading file ' matFiles{ii} ' (' num2str(ii) ' of ' num2str(nFiles) ')'])
    
    % Load data
    load([dataPath matFiles{ii}],'data');
    
    % Downsample to reduce data size
    if downsamp_factor ~= 1
        disp(['Downsampling, keeping every ' num2str(downsamp_factor) ' columns'])
        data = downsample(data.',downsamp_factor).';
    end
    xStep = xStep_orig * downsamp_factor;
    
    % FFT transform along time dimension, extract passband
    disp('Fourier transform along time dimension, extracting passband')
    [Pox, omega] = calc_spec_omega(data,fs,fLow,fHigh);
    
    % Save to file
    [~,origFileName,~] = fileparts(matFiles{ii});
    saveFileName = [savePath origFileName '_preprocessed.mat'];
    disp('Saving preprocessed data to')
    disp(['    ' saveFileName])
    save(saveFileName,'Pox','omega','xStep')
end

function [Pox,omega] = calc_spec_omega(ptx,fs,fLow,fHigh)
% calc_spec_omega - calculate (bandpass) omega spectrum for 2D dataset
%
%   Usage:
%   [Pox,omega] = calc_spec_omega(ptx,fs,fLow,fHigh)
%
%   Input parameters:
%   ptx     -    MxN matrix, dimensions (t,x)
%   fs      -    Sampling frequency [Hz]
%   fLow    -    Lower limit of transducer frequency band [Hz]
%   fHigh   -    Higher limit of transducer frequency band [Hz]
%
%   Output parameters:
%   Pox     -   Spectrum matrix, dimensions (omega,x)
%   omega   -   vector of omega values in output
%
%   Notes:
%   The function calculates the Fourier spectrum along the first dimension
%   of ptx, corresponding to the time axis. The spectrum is calculated at
%   nextpow2(M) points, where M is the number of time samples (the fft()
%   function is most efficient for input sizes that are powers of 2).
%   Optionally, this number can be modified using the 'tFftMult' parameter.
%   The FFT size is then given by tFftMult * nextpow2(M).
%
%   After calculating the "full" FFT, the passband (given by fLow and
%   fHigh) is extracted. Note that only positive omega (one-sided spectrum)
%   are used, i.e. the spectrum is extracted from the "first half" of the
%   output from fft().
%

%% Calculate FFT size, omega vector and omega bassband
nFFTt = 2^nextpow2(size(ptx,1));                                % FFT size, time axis
omega = ((0:(nFFTt-1)) - floor(nFFTt/2))'*((2*pi*fs)/nFFTt);    % Freq. vector
omega = ifftshift(omega);                                       % Shift as fft
omegaBandIndex = (omega >= 2*pi*fLow) & (omega <= 2*pi*fHigh);  % Index to passband

%% Fourier transform along time dimension
Pox = fft(ptx,nFFTt,1);

%% Extract passband
Pox = Pox(omegaBandIndex,:);
omega = omega(omegaBandIndex);
end