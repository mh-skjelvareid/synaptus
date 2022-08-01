close all
clearvars

%% Show test header (useful when run as part of test suite)
disp('-----------------------------------------------------------------')
disp('-----           SYNAPTUS TOOLBOX TEST: psm_test2D           -----')
disp('-----           Phase shift migration, 2D dataset           -----')
disp('-----------------------------------------------------------------')
disp(' ')

%% Add path to necessary functions
toolboxPath=fileparts(fileparts(mfilename('fullpath'))); %Get the toolbox path

% Add core and misc path
addpath(fullfile(toolboxPath,'core'),fullfile(toolboxPath,'misc'));

%% Load data
disp('Loading data');
load(fullfile(toolboxPath,'datasets','LineScan2D_PinsPlexiAluSDH.mat'),...
    'fs','ptx','tDelay','fLow','fHigh','thick','xStep','cc');

%% Parameters
[nT,nX] = size(ptx);                % Size of B-scan matrix
tPlot = (0:(nT-1))/fs + tDelay;     % Time vector for y axis
xPlot = (0:(nX-1))*xStep;           % X vector for x axis

%% Focus using PSM
disp('Processing data')
tic
[im,xIm,zIm] = psm(ptx,fs,tDelay,cc,thick,fLow,fHigh,xStep);
toc

%% Plot raw data
figure
imagesc(xPlot*1e3,tPlot*1e6,logImage(abs(hilbert(ptx))))
colorbar
caxis([-30 0])
xlabel('X [mm]')
ylabel('Time [\mu s]')
title('Raw data')

%% Plot focused image
figure
for ii = 1:3
    subplot(3,1,ii)
    imagesc(xIm*1e3,zIm{ii}*1e3,logImage(im{ii}))
    caxis([-30 0])
    ylabel('Z [mm]')
    title(['Focused image, layer ' num2str(ii)])
end
xlabel('X [mm]')

%% Add blank line (nicer formatting for test text output).
disp(' ')
