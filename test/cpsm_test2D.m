close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----        SYNAPTUS TOOLBOX TEST: cpsm_test2D             -----')
disp('-----    Cylindrical phase shift migration, 2D dataset      -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

%Add core and misc path
addpath(fullfile(toolboxPath,'core'),fullfile(toolboxPath,'misc'));

%% Load data, extract 2D slice and convert to double precision
disp('Loading data')
load(fullfile(toolboxPath,'datasets','CylScan3D_LeadSpheres.mat'),'ptpz','fs','phiStep','zStep',...
    'tDelay','r0','fLow','fHigh','cc');
zSliceInd = 30;                         % Index for 2D slice from 3D dataset
ptp = double(ptpz(:,:,zSliceInd));      % Extract slice, convert to double
clear ptpz;                             % Clear 3D data to save memory

%% Set processing parameters
rStep = (cc/2)/(fHigh-fLow);                % Range resolution
rStart = r0 + (tDelay*(cc/2));              % Start of focused image
rEnd = rStart + (size(ptp,1)/fs)*(cc/2);    % End of focused image

transFunc = 'haun';     % Transfer function, options 'haun', 'gardner', 'exact'

%% Focus
disp('Processing data')
tic
[im,phiIm,rIm] = cpsm(ptp,fs,tDelay,cc,fLow,fHigh,phiStep,r0,rStep,...
    rStart,rEnd,'transFunc',transFunc);
toc

%% Plot raw data and focused image
figure
polRawDataImage(abs(hilbert(ptp)),tDelay+(0:(size(ptp,1)-1))/fs,phiIm,...
    [],[],r0,cc)
figure
polFocusedImage(abs(im),phiIm,rIm,[],[])

%% Add blank line (nicer formatting for test text output).
disp(' ')